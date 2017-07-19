# This script is based on the old MAIN_ files
# that created the x and y matrices separately.
# However, after some analysis we discovered that it is best if
# we first create the x matrix and afterwards the y matrix.
# The reason is that when x is specifies some platforms must be excluded
# due to low gene coverage. 
# This in turn will change the sample set for the analysis, which
# can lead to a different set of labels for the analysis.

#########################################################
# Set the working directory
try({setwd('Adeptus2/data/expression_data/')})
##########################################################

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIN_NUMBER_OF_GENES = 10000 # Determines which platforms will be excluded
ADD_TCGA = T
ADD_GTEX = T
PLATFORMS_TO_EXCLUDE = c() # The user can specify GEO platforms to be excluded
MAX_PERCENT_NAS_IN_PROFILE = 10 # Samples that have more than this percentage of NAs are removed
MIN_NUM_DATASETS_IN_LABEL = 3 # Determines which labels to include
MIN_NUM_SAMPLES_IN_LABEL = 100 # Determines which labels to include
UNIX = T
NUM_CORES = "Auto" # For unix only
TRANSFORM_TO_ENTREZ = F # False means that we assume that the entrez matrices are already present
TRANSFORM_TO_RANKS = F # False assumes that the output file already contains the normalized ranks
LOCAL_RLIBS_DIR = "Adeptus2/data/rlibs" # In Unix we sometimes have a local dir with R libraries
# In each expression data dir we have a sample annotation file, these specify important columns
# Other assumptions on these files: first columns are sample ids
ANNOT_FILE_TISSUE_COLUMN_IND = 6 
ANNOT_FILE_DOID_COLUMN_IND = 9
# When counting datasets per label, ignore datasets with less than this number of positive samples:
MIN_NUM_SAMPLES_IN_DATASET_FOR_COUNT_IN_LABEL=10
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_supervised_data_file = '../Supervised_data_atleast_10000_genes_max_percet_na_10.RData'
# Set the library paths
.libPaths( c( .libPaths(), LOCAL_RLIBS_DIR ) )
#biocLite('AnnotationDbi',lib='../rlibs',lib.loc='../rlibs',dependencies=T)
#biocLite('BiocGenerics',lib='../rlibs',lib.loc='../rlibs',dependencies=T)

##########################################################
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
#currx[4,] = NA
#> currx
#        GSM750476 GSM750477 GSM750478 GSM750479
#2998453      1.69      1.79      2.13      2.68
#3894047      6.98      7.39      7.64      7.92
#2711139      3.46      3.30      3.83      3.25
#4037638        NA        NA        NA        NA
#2315554      6.68      6.38      6.27      6.35
#> probe2genes_conv(currx,curr_mapping1,na.rm=T)
#          GSM750476 GSM750477 GSM750478 GSM750479
#100133161      1.69      1.79      2.13      2.68
#140849         6.98      7.39      7.64      7.92
#402483         1.69      1.79      2.13      2.68
#344905         3.46      3.30      3.83      3.25
#254173         6.68      6.38      6.27      6.35

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
##########################################################

# Load the platform data
load('gpl_mappings_to_entrez.RData')
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


############## Transform probes to entrez ####################
# In windows: use a loop
if(!UNIX && TRANSFORM_TO_ENTREZ){
	for(i in 1:length(curr_files)){
		rdata_file = curr_files[i]
		out_file = out_files[i]
		print(rdata_file)
		gpl = gpls[i]
		new_mat = transform_matrix_into_genes(rdata_file,gpl,gpl_mappings_entrez2probes,out_file)
	}
}else if (TRANSFORM_TO_ENTREZ) {
	# In unix: use mclapply
	library(parallel)
	if (NUM_CORES == "Auto"){
		NUM_CORES = length(curr_files)
	}
	objs = mclapply(1:length(curr_files),
		function(x,a,b,c,d){transform_matrix_into_genes(a[x],b[x],c,d[x])},
		a = curr_files,b=gpls,c=gpl_mappings_entrez2probes,d=out_files,
		mc.cores=NUM_CORES)
}

####################################################################
# QA: check GPL5175
#curr_raw_files = curr_files[grepl("GPL5175",curr_files)]
#gsms1 = get(load(curr_raw_files[1]))
#gsms2 = get(load(curr_raw_files[2]))
#dim(gsms1)
#dim(gsms2)
#v1 = get_gsm_profile(colnames(gsms1)[1])
#v2 = get_gsm_profile(colnames(gsms2)[1])
####################################################################
# Preprocessing before moving to rank-based scores
get_list_intersect<-function(l){
	inter = l[[1]]
	for(j in 1:length(l)){inter = intersect(inter,l[[j]])}
	return(inter)
}
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
xxx = runif(1000)
plot(xxx,getRankedBasedProfile(xxx),ylab="Rank transformed profile",xlab="Original profile")

if (TRANSFORM_TO_RANKS){
	# Analyze the platforms by the gene coverage
	platform_matrix_to_geneset = list()
	platform_matrix_to_samples = list()
	NA_count_per_sample = list()
	for(i in 1:length(out_files)){
		out_file = out_files[i]
		print(out_file)
		gpl = gpls[i]
		curr_x = get(load(out_file))
		num_NAs_in_cols = colSums(is.na(curr_x))
		print(table(is.na(curr_x))/length(curr_x))
		platform_matrix_to_geneset[[out_file]] = rownames(curr_x)
		NA_count_per_sample[[out_file]] = num_NAs_in_cols
		platform_matrix_to_samples[[out_file]] = colnames(curr_x)
	}
	perc_NA_count_per_sample = c()
	for(i in 1:length(out_files)){
		out_file = out_files[i]
		perc_NA_count_per_sample = c(perc_NA_count_per_sample,
			NA_count_per_sample[[out_file]]/length(platform_matrix_to_geneset[[out_file]]))
	}
	#length(perc_NA_count_per_sample)
	#print(sort(perc_NA_count_per_sample,decreasing=T)[1:10])
	gsms_to_exclude = names(perc_NA_count_per_sample)[perc_NA_count_per_sample>=MAX_PERCENT_NAS_IN_PROFILE/100]
	## Check direct download of "GSM571726"
	#v = get_gsm_profile("GSM571726")
	#table(is.na(v))/length(v)
	#perc_NA_count_per_sample["GSM571726"]

	removal_order = names(sort(sapply(platform_matrix_to_geneset,length)))
	rem_ind = 0
	curr_coverage = length(get_list_intersect(platform_matrix_to_geneset))
	print(paste(0,curr_coverage))
	while(rem_ind<length(removal_order)-1){
		curr_removals = removal_order[1:(rem_ind+1)]
		to_rem = is.element(names(platform_matrix_to_geneset),set=curr_removals)
		curr_coverage = length(get_list_intersect(platform_matrix_to_geneset[!to_rem]))
		print(paste(rem_ind+1,curr_coverage))
		rem_ind = rem_ind+1
		if(curr_coverage > MIN_NUMBER_OF_GENES){break}
	}

	files_to_exclude = c()
	if(rem_ind>0){
		files_to_exclude = curr_removals
	}
	for(f in out_files){
		for(gpl in PLATFORMS_TO_EXCLUDE){
			if(grepl(pattern=gpl,f)){
				files_to_exclude = union(files_to_exclude,f)
			}
		}
	}
	#print(length(unlist(platform_matrix_to_samples[files_to_exclude])))

	out_files = setdiff(out_files,files_to_exclude)
	array_gene_set = get_list_intersect(platform_matrix_to_geneset[out_files])
	adeptus2_gene_set = array_gene_set

	# load the gtex and tcga datasets
	non_geo_datasets = list()
	if (ADD_GTEX){
		non_geo_datasets[['gtex']] = get(load('gtex/raw_expression_profiles.RData'))
	}
	if(ADD_TCGA){
		non_geo_datasets[['tcga']] = get(load('tcga/raw_expression_profiles.RData'))
	}	
	if (length(non_geo_datasets)>0){
		non_geo_gene_sets = lapply(non_geo_datasets,rownames)
		non_geo_gene_sets_inter = get_list_intersect(non_geo_gene_sets)
		# As a secondary analysis: consider removing GPL96
		adeptus2_gene_set = intersect(non_geo_gene_sets_inter,array_gene_set)
	}

	# Transform to ranks using the preprocessing results above
	x_arrays = c()
	for(f in out_files){
		print(f)
		e_x = get(load(f))
		e_x = e_x[,setdiff(colnames(e_x),gsms_to_exclude)]
		e_x_ranks = t(apply(e_x,2,getRankedBasedProfile,fs=adeptus2_gene_set))
		e_x_ranks = e_x_ranks[,adeptus2_gene_set]
		x_arrays = rbind(x_arrays,e_x_ranks)
		print(dim(x_arrays))
	}
	for(f in names(non_geo_datasets)){
		e_x = non_geo_datasets[[f]]
		e_x_ranks = t(apply(e_x,2,getRankedBasedProfile,fs=adeptus2_gene_set))
		e_x_ranks = e_x_ranks[,adeptus2_gene_set]
		x_arrays = rbind(x_arrays,e_x_ranks)
		print(dim(x_arrays))
	}
	x = x_arrays
	rm(x_arrays);gc()
	save(x,file=output_supervised_data_file)
}
else{
	load(output_supervised_data_file)
}

####################################################################
####################################################################
####################################################################
# Get the Y matrix for the current generated x
# Load required data for this part
try({library(DO.db)})
try({library(IRanges,lib.loc=LOCAL_RLIBS_DIR)})
try({library(DO.db,lib.loc=LOCAL_RLIBS_DIR)})
dos = as.list(DOTERM)
ancs = as.list(DOANCESTOR)
tissue_slim_mapping_file = '../sample_labels/new_tissue_slim_mapping.txt'
tissue_slim_mapping = as.matrix(read.delim(tissue_slim_mapping_file,sep='\t',row.names=1))

dirs = c('GEO','GSE60862','gtex','tcga')
term2name<-function(n){
	if (is.null(DOTERM[[n]])){term=n}else{term = Term(DOTERM[[n]])}
	return (term)
}

sample2terms = list();sample2doids = list()
sample2study = c();sample2tissue_slim = c()
for(dir in dirs){
	curr_annots = as.matrix(read.delim(paste(dir,"/sample_annotation.txt",sep='')))
	
	if(dir=="gtex"){
		xgtex_samples = rownames(x)[grepl("gtex",rownames(x),ignore.case=T)]
		curr_annots[,1] = gsub(curr_annots[,1],pattern='-',replace='\\.')
	}
	
	to_keep = is.element(curr_annots[,1],set=rownames(x))
	curr_annots = curr_annots[to_keep,]
	# correct the tissue names to lowercase
	curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND] = tolower(curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND])
	curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND] = tissue_slim_mapping[curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND],1]
	sample2tissue_slim[curr_annots[,1]] = curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND]
		
	# add the tissue slim of each sample
	curr_tissues = curr_annots[,ANNOT_FILE_TISSUE_COLUMN_IND]
	# add all control samples
	curr_control_samples = curr_annots[curr_annots[,ANNOT_FILE_DOID_COLUMN_IND]=="",1]
	sample2terms[curr_control_samples] = paste('control',curr_tissues[curr_annots[,ANNOT_FILE_DOID_COLUMN_IND]==""],sep=';')
	# add the study of each sample
	sample2study[curr_annots[,1]] = curr_annots[,2]
	curr_dos = unique(curr_annots[,c(ANNOT_FILE_DOID_COLUMN_IND,ANNOT_FILE_TISSUE_COLUMN_IND)])
	if(length(curr_dos)==0){next}
	
	# analyze the current non control samples
	curr_dos = curr_dos[curr_dos[,1]!="",]
	if(is.null(dim(curr_dos))){curr_dos=matrix(curr_dos,nrow=1)}
	for(i in 1:nrow(curr_dos)){
		tissue = curr_dos[i,2]
		doids = curr_dos[i,1]
		# get all current samples with the same phenotype and from the same tissue
		curr_samples = curr_annots[curr_annots[,9]==doids & curr_annots[,6]==tissue,1]
		# split the doid field into all available do terms
		doids = strsplit(doids,split=";")[[1]]
		doids = gsub(doids,pattern = ' ', replace= '')
		# add all ancestors
		doids = unique(union(doids,unique((unname(unlist(ancs[doids]))))))
		# the name of a class: doid;disease_name;tissue_slim
		full_class_names = paste(doids,sapply(doids,term2name),sep=';')
		full_class_names = unique(paste(full_class_names,tissue,sep=';'))
		for(s in curr_samples){
			sample2terms[[s]] = full_class_names
			sample2doids[[s]] = doids
		}
	}
}

print(table(sample2study))
print(table(sample2tissue_slim))
print(length(sample2terms))

# Here we count for each possible term
# the number of datasets that have at least 10 samples
all_classes = unique(unlist(sample2terms))
all_studies = unique(sample2study)
class2num_datasets = rep(0,length(all_classes));names(class2num_datasets) = all_classes
all_doids = unique(sapply(all_classes,function(x)strsplit(x,split=';')[[1]][1]))
all_doids = all_doids[grepl("DOID",all_doids)]
doid2num_datasets = rep(0,length(all_doids));names(doid2num_datasets) = all_doids
for(study in all_studies){
	curr_samps = names(which(sample2study==study))
	curr_classes_counts = table(unlist(sample2terms[curr_samps]))
	curr_classes_counts = curr_classes_counts[curr_classes_counts>=MIN_NUM_SAMPLES_IN_DATASET_FOR_COUNT_IN_LABEL]
	if(length(curr_classes_counts)==0){next}
	class2num_datasets[names(curr_classes_counts)] = class2num_datasets[names(curr_classes_counts)]+1
	curr_classes = names(curr_classes_counts)
	curr_doids = unique(sapply(curr_classes,function(x)strsplit(x,split=';')[[1]][1]))
	curr_doids = intersect(all_doids,curr_doids)
	doid2num_datasets[curr_doids] = doid2num_datasets[curr_doids]+1
}
terms_for_analysis = names(which(class2num_datasets>=MIN_NUM_DATASETS_IN_LABEL))
doids_for_analysis = names(doid2num_datasets)[doid2num_datasets>=MIN_NUM_DATASETS_IN_LABEL ]
y = matrix(0,nrow=length(sample2study),ncol=length(terms_for_analysis)+length(doids_for_analysis))
rownames(y) = names(sample2study)
colnames(y) = c(terms_for_analysis,doids_for_analysis)
for(i in 1:length(sample2terms)){
	curr_samp = names(sample2terms)[i]
	curr_terms = sample2terms[[i]]
	curr_terms_doids = sapply(curr_terms,function(x)strsplit(x,split=';')[[1]][1])
	curr_terms = intersect(curr_terms,colnames(y))
	curr_terms_doids = intersect(curr_terms_doids,colnames(y))
	y[curr_samp,curr_terms]=1
	y[curr_samp,curr_terms_doids]=1
	if(i %% 500==0){print(i)}
}
y = y[,colSums(y)>=MIN_NUM_SAMPLES_IN_LABEL]


# Remove labels that have the same sample sets: keep the ancestor only
#source("http://bioconductor.org/biocLite.R")
#biocLite('DO.db',lib='rlibs')
#library(DO.db,lib.loc='rlibs')
#try(library(DO.db))
doancs = as.list(DOANCESTOR)
labels_to_r = c()
# We compare a pair of columns. 
# If they share the same positive sample sets then
# we have three cases:
# 1. Both labels have a tissue information
#	Solution: remove the ancestor 
# 2. One has a tissue and the other does not
#	Solution: remove the one without the tissue
# 3. Both do not have the tissue
#	Solution: remove the ancestor
# We check these cases as follows:
# step1: if both DOs are the same - take the one
#	 with the tissue
# step2: otherwise, remove the ancestor
for(i in 1:ncol(y)){
	for(j in 1:ncol(y)){
		if(i==j){next}
		curr_agreement = sum(y[,i]==y[,j])
		if(curr_agreement <= nrow(y)-5){next}
		print(colnames(y)[c(i,j)])
		label1 = strsplit(colnames(y)[i],split=';')[[1]]
		label2 = strsplit(colnames(y)[j],split=';')[[1]]
		do1 = label1[1];do2=label2[1]
		has_tissue1 = length(label1)>=3
		has_tissue2 = length(label2)>=3
		if (do1==do2){
			print(paste(i,j,paste(colnames(y)[i],colnames(y)[j],sep = "   ;   "),setp=','))
			if (!has_tissue1){labels_to_r = c(labels_to_r,i)}
			if (!has_tissue2){labels_to_r = c(labels_to_r,j)}
		}
		# Cases 1 + 3
		# do1 is an ancestor of do2: add do1 to the excluded
		else if (is.element(do1,set=doancs[[do2]])){
			if (has_tissue1 && has_tissue2){print(paste(i,j,paste(label1[2],label2[2],sep = "   <-   "),setp=','))}
			else{print(paste(i,j,paste(colnames(y)[i],colnames(y)[j],sep = "   <-   "),setp=','))}
			labels_to_r = c(labels_to_r,i)
		}
	}
	labels_to_r = unique(labels_to_r)
}
y = y[,-labels_to_r]
save(x,y,sample2study,sample2terms,sample2tissue_slim,file=output_supervised_data_file)

# Auxiliary method to add edges and labels for the analysis
# For a node D we add a dummy node D_star with all of D samples that
# are not in its children. We also add the the corresponding label to a
# returned matrix
add_star_dummy_nodes<-function(edges,Y,minSize=10){
	all_nodes = intersect(unique(c(edges)),colnames(Y))
	m = c();newnames=c()
	for (node in all_nodes){
		children = intersect(colnames(Y),get_children(node,edges))
		if(length(children)==0){next}
		S1 = rownames(Y)[Y[,node]==1]
		if(length(children)==1){
			S1_aux = rownames(Y)[Y[,children]>0]
		}
		else{
			S1_aux = rownames(Y)[rowSums(Y[,children])>0]
		}
		S1 = setdiff(S1,S1_aux)
		if(length(S1)<minSize){next}
		v = rep(0,nrow(Y))
		names(v) = rownames(Y)
		v[S1] = 1
		m = cbind(m,v)
		newnames = c(newnames,paste(node,"_star",sep=""))
		edges = rbind(edges,c(node,paste(node,"_star",sep="")))
	}
	rownames(m) = rownames(Y)
	colnames(m) = newnames
	return(list(edges=edges,m=m))
}

#Auxiliary methods for analysis of edge lists
# of the form: parent,son
get_parents<-function(node,edges){
	return(unique(edges[edges[,2]==node,1]))
}
get_children<-function(node,edges){
	return(unique(edges[edges[,1]==node,2]))
}
get_siblings<-function(node,parent,edges){
	return(setdiff(get_children(parent,edges),node))
}
source("../../rscripts/FUNCTIONS_multilabel_classification.R")
Y = y
Y = Y[,!grepl("_star$",colnames(Y),perl=T)]
library(bnlearn)
DOnet = getBayesianNetwork(colnames(Y))
mapping = DOnet[[2]]
asgraph = get_graph_for_cytoscape(arcs(DOnet[[1]]),mapping)
obj = add_star_dummy_nodes(asgraph,Y)
Y = cbind(Y,obj$m)

# exclude labels with too few samples or datasets
to_rem = c()
for(j in 1:ncol(Y)){
	S = names(which(Y[,j]==1))
	DS = sample2study[S]
	if(length(S) < MIN_NUM_SAMPLES_IN_LABEL || length(DS)< MIN_NUM_DATASETS_IN_LABEL){
		to_rem = c(to_rem,j)
	}
}
Y = Y[,-to_rem]
y = Y
grepl("_star$",colnames(Y),perl=T)
save(x,y,sample2study,sample2terms,sample2tissue_slim,file=output_supervised_data_file)

### Load the data and get some statistics ###
load(output_supervised_data_file)
samples = intersect(rownames(x),rownames(y))
table(sample2study[samples])
length(table(sample2study[samples]))
length(samples)
gtex_samples = grepl("gtex",samples,ignore.case=T)
tcga_samples = grepl("tcga",samples,ignore.case=T)
array_samples = (!gtex_samples) & (!tcga_samples)
table(array_samples)
annots_inf = get(load('../adeptus2_pheno_data.RData'))
unique(sample2tissue_slim[samples])
blood_samples = grepl("blood",sample2tissue_slim[samples])

array_samp = sample(which(array_samples&blood_samples))[1:400]
rnaseq_samp = sample(which(!array_samples&blood_samples))[1:400]
corrs = cor(t(x[c(array_samp,rnaseq_samp),]))
get_sq_mat_vals<-function(x){x[lower.tri(x)]}
par(mfrow=c(3,1))
hist(get_sq_mat_vals(corrs[1:400,1:400]),xlim=c(-0.5,1),main="Correlation between microarray profiles",xlab="Pearson rho")
hist(get_sq_mat_vals(corrs[401:800,1:400]),xlim=c(-0.5,1),main="Correlation between microarray and RNA-seq profiles",xlab="Pearson rho")
hist(get_sq_mat_vals(corrs[401:800,401:800]),xlim=c(-0.5,1),main="Correlation between RNA-seq profiles",xlab="Pearson rho")




















