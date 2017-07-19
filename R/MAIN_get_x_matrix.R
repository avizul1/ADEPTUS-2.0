setwd('Adeptus2/data/expression_data/')

# Load the platform data
load('gpl_mappings_to_entrez.RData')

########### Functions ##############
# Convert a matrix x with probe ids as names into a new
# matrix of genes. 
# h is the list that maps each gene to its probe ids
# f is the function to apply, default is mean
probe2genes_conv<-function(x,h,f=mean, genes=NULL,...){
	if (is.null(genes)){genes = names(h)}
	y = sapply(genes,applyFunctionOnGroup,h=h,f=f,x=x,...)
	y = t(as.matrix(data.frame(y,check.names=F)))
	return (y)
}
# Apply a function f on values in the vector x (or a matrix)
# that are under the names of h[[g]]
# For example if h maps genes to probes then we can apply
# the mean value of the probes by setting f=mean
# We assume that all mapped probes are necessarily represented as 
# the names in x.
applyFunctionOnGroup<-function(g,h,x,f,intsct=F,...){
	curr_names = h[[g]]
	if (is.null(dim(x))){return (f(x[curr_names]))}
	if (length(curr_names)==1){return (sapply(x[curr_names,],f))}
	return (apply(x[curr_names,],2,f,...))
}

# clean a mapping list by the row names in m
clean_mapping<-function(m,h){
	newh = lapply(h,function(x,y)intersect(x,y),y=rownames(m))
	names(newh) = names(h)
	newh = newh[sapply(newh,length)>0]
	return(newh)
}

#### Tests
#load('GEO/GPL5175_expression_profiles.RData')
#curr_mapping = gpl_mappings_entrez2probes[["GPL5175"]]
#curr_mapping = clean_mapping(gsms_mat,curr_mapping)
#curr_mapping1 = curr_mapping[1:5]
#curr_probes = unique(unname(unlist(curr_mapping1)))
#currx = gsms_mat[curr_probes,1:4]
#probe2genes_conv(currx,curr_mapping1)
# Comparison of the results:
#> probe2genes_conv(currx,curr_mapping1)
#          GSM750476 GSM750477 GSM750478 GSM750479
#100133161      1.69      1.79     2.130     2.680
#140849         6.98      7.39     7.640     7.920
#402483         1.69      1.79     2.130     2.680
#344905         7.18      6.77     7.045     6.755
#254173         6.68      6.38     6.270     6.350
#> currx
#        GSM750476 GSM750477 GSM750478 GSM750479
#2998453      1.69      1.79      2.13      2.68
#3894047      6.98      7.39      7.64      7.92
#2711139      3.46      3.30      3.83      3.25
#4037638     10.90     10.24     10.26     10.26
#2315554      6.68      6.38      6.27      6.35
#> curr_mapping1
#$`100133161`
#[1] "2998453"
#
#$`140849`
#[1] "3894047"
#
#$`402483`
#[1] "2998453"
#
#$`344905`
#[1] "2711139" "4037638"
#
#$`254173`
#[1] "2315554"
###########
# A function that collects data on the NA cells in a matrix
collect_na_stats<-function(x){
	na_x = is.na(x)
	stats = list()
	stats[['percent NAs']] = sum(c(na_x))/length(na_x)
	stats[['num rows with NAs']] = sum(apply(na_x,1,any))
	stats[['num cols with NAs']] = sum(apply(na_x,2,any))
	stats[['row NA counts']] = rowSums(na_x)
	stats[['col NA counts']] = colSums(na_x)
	return(stats)
}
# A function that transforms a probe level data into a gene level matrix
# OPTIONAL: set out_path to save the new gene matrix to an RData file
# ... additional parameters for probe2genes_conv
transform_matrix_into_genes<-function(rdata_path,gpl,gpl_mappings_entrez2probes,out_path=NULL,...){
	gsms_mat = get(load(rdata_path))
	mode(gsms_mat) = 'numeric'
	print('gsm matrix was loaded')
	gsms_mat_na_stats = collect_na_stats(gsms_mat)
	curr_mapping = gpl_mappings_entrez2probes[[gpl]]
	curr_mapping = clean_mapping(gsms_mat,curr_mapping)
	print('difference before and after cleaning the gpl mapping list:')
	print(paste(length(gpl_mappings_entrez2probes[[gpl]]),length(curr_mapping)))
	print('merging genes by averaging probes')
	entrez_mat = probe2genes_conv(gsms_mat,curr_mapping,na.rm=T)
	print('done')
	print(dim(entrez_mat))
	entrez_mat_na_stats = collect_na_stats(entrez_mat)
	if (!is.null(out_path)){
		print('saving the new matrix into a file')
		save(entrez_mat,file=out_path)
	}
	return(list(gsms_mat_na_stats=gsms_mat_na_stats,entrez_mat=entrez_mat,entrez_mat_na_stats=entrez_mat_na_stats))
}

# Start with the GEO datasets
curr_files = list.files('GEO')
curr_files = curr_files[grepl('expression_profiles.RData',curr_files)]
gpls = sapply(curr_files,function(x)strsplit(x,split='_')[[1]][1])
curr_files = paste('GEO/',curr_files,sep='')
out_files = paste('GEO/',gpls,'_entrez_merged_profiles.RData',sep='')
# Add GSE60862
curr_files = c(curr_files,"GSE60862/GPL5175_expression_profiles.RData")
gpls = c(gpls,"GPL5175")
out_files = c(out_files,paste('GSE60862/','GPL5175','_entrez_merged_profiles.RData',sep=''))

# In windows: use a loop
for(i in 1:length(curr_files)){
	rdata_file = curr_files[i]
	print(rdata_file)
	out_file = out_files[i]
	gpl = gpls[i]
	new_mat = transform_matrix_into_genes(rdata_file,gpl,gpl_mappings_entrez2probes,out_file)
}

# In unix: use mclapply
library(parallel)
objs = mclapply(1:length(curr_files),
	function(x,a,b,c,d){transform_matrix_into_genes(a[x],b[x],c,d[x])},
	a = curr_files,b=gpls,c=gpl_mappings_entrez2probes,d=out_files,
	mc.cores=length(curr_files))

## Check some stats: go over the files
#gsm_stats = c()
#entrez_stats = c()
#gene_sets = c();sample_sets=c()
#for(f in out_files){
#	print(f)
#	f_x = get(load(f))
#	stats = collect_na_stats(f_x)
#	print(stats[1:3])
#	entrez_stats = rbind(entrez_stats,unlist(stats[1:3]))
#	gene_sets[[f]] = rownames(f_x)
#	sample_sets[[f]] = colnames(f_x)
#}
###### Analysis of the gene sets of the actual matrices #####
# Based on the analysis below we decided to exclude: GPL6102
# as the downloaded profiles had only ~15,000 probes unlike the 
# GPL data that indicated >30k probes. Using this platform leads
# to less than 8k genes in the intersection between the platforms
# intersect between the gene sets
get_list_intersect<-function(l){
	inter = l[[1]]
	for(j in 1:length(l)){inter = intersect(inter,l[[j]])}
	return(inter)
}
#length(get_list_intersect(gene_sets))
#platform2num_samples = c()
#platform2num_genes = c()
## try removing single gpls
#for(gpl in unique(gpls)){
#	print(gpl)
#	to_rem = grepl(names(gene_sets),pattern=gpl)
#	print(sum(to_rem))
#	print(length(get_list_intersect(gene_sets[!to_rem])))
#	platform2num_samples[gpl] = length(unlist(sample_sets[to_rem]))
#	platform2num_genes[gpl] = length(gene_sets[to_rem][[1]])
#}
## A problem was detected in "GPL6102": the number of genes in the entrez
## matrix is too small
#gsm_x = get(load("GEO/GPL6102_expression_profiles.RData"))
#curr_mapping = gpl_mappings_entrez2probes[["GPL6102"]]
#length(unique(unlist(curr_mapping)))
#colnames(gsm_x)
# Manual inspection of the profiles in the GEO website and the series matrix
# confirms that the profiles indeed have <16k probes.
########################################
# QA: case test
#e_x = get(load("GEO/GPL6947_entrez_merged_profiles.RData"))
#gsm_x = get(load("GEO/GPL6947_expression_profiles.RData"))
#hist(colSums(is.na(e_x)))
#sort(colSums(is.na(e_x)))
## looak at "GSM983584"
#which(is.na(e_x[,"GSM983584"]))
#curr_mapping = gpl_mappings_entrez2probes[["GPL6947"]]
#gsm_x[curr_mapping[["653321"]],"GSM983584"]
## do we have cells in e_x that are na whose corresponding data in
## gsm_x is not all NAs?
## If such cases exist then the following should print some output
#rows_with_NAs = which(apply(is.na(e_x),1,any))
#for (r in rows_with_NAs){
#	curr_e = rownames(e_x)[r]
#	curr_probes = curr_mapping[[curr_e]]
#	curr_NA_samples = which(is.na(e_x[r,]))
#	if(!all(is.na(gsm_x[curr_probes,curr_NA_samples]))){
#		print(curr_e)
#	}
#}
#######################################

# remove GPL6102
to_rem = grepl(curr_files,pattern='GPL6102')
out_files = out_files[!to_rem]
gpls = gpls[!to_rem]
curr_files = curr_files[!to_rem]
gene_sets = gene_sets[!to_rem]
array_gene_set = get_list_intersect(gene_sets)
sample_sets = sample_sets[!to_rem]
length(unique(unlist(sample_sets)))

# load the gtex and tcga datasets
non_geo_datasets = list()
non_geo_datasets[['gtex']] = get(load('gtex/raw_expression_profiles.RData'))
non_geo_datasets[['tcga']] = get(load('tcga/raw_expression_profiles.RData'))
sapply(non_geo_datasets,function(x,y)length(intersect(rownames(x),y)),y=array_gene_set)
non_geo_gene_sets = lapply(non_geo_datasets,rownames)
non_geo_gene_sets_inter = get_list_intersect(non_geo_gene_sets)

# As a secondary analysis: consider removing GPL96
adeptus2_primary_gene_set = intersect(non_geo_gene_sets_inter,array_gene_set)

##### Merge into one big matrix based on rankings ########

# This method receives as input a gene expression profile x
# and a set of feature names to consider fs.
# The method transforms the profile of the considered genes
# to ranks and then to exponential score of the ranks.
# Difference from Adeptus1 - to ease interpretation we now report
# max(w)-w. This way we keep the features of the function (i.e., the difference
# between highly expressed genes is higher that the difference between lowly expressed 
# genes), but the highly expressed genes have the top scores.
getRankedBasedProfile<-function(x,fs=NULL){
	if (!is.null(fs)){x = x[fs]}
	# get the ranks in decreasing order
	N = length(x)
	rs = N-rank(x,ties.method = "average")+1
	w_rs = rs*exp(-rs/N)
	names(w_rs) = names(rs)
	w_rs = -w_rs
	w_rs = w_rs+abs(min(w_rs))
	return (w_rs)
}

# QA
#x = runif(1000)
#par(mfrow=c(1,2))
#plot(x,getRankedBasedProfile(x))
#boxplot(getRankedBasedProfile(x))

x_arrays = c()
for(f in out_files){
	print(f)
	e_x = get(load(f))
	e_x_ranks = t(apply(e_x,2,getRankedBasedProfile,fs=adeptus2_primary_gene_set))
	e_x_ranks = e_x_ranks[,adeptus2_primary_gene_set]
	x_arrays = rbind(x_arrays,e_x_ranks)
	print(dim(x_arrays))
}
for(f in names(non_geo_datasets)){
	e_x = non_geo_datasets[[f]]
	e_x_ranks = t(apply(e_x,2,getRankedBasedProfile,fs=adeptus2_primary_gene_set))
	e_x_ranks = e_x_ranks[,adeptus2_primary_gene_set]
	x_arrays = rbind(x_arrays,e_x_ranks)
	print(dim(x_arrays))
}
x = x_arrays
rm(x_arrays);gc()
save(x,file='x_all_array_platforms.RData')























