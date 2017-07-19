# Set the working directory to the dir with the sample_annotation.txt file that
# you want to analyze
working_dir = 'Adeptus2/data/expression_data'
setwd(working_dir)

# Required libraries: GEOquery and the Human database
library('GEOquery')
library('org.Hs.eg.db')
library('GEOquery',lib.loc='rlibs')
library('org.Hs.eg.db',lib.loc='rlibs')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_GSM_ANNOT_FILE = 'GEO/sample_annotation.txt'
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# our flow gets as input an annotation table
# the first three columns are: GSM ids, GSE or GDS ids, and the GPL ids
# we use these data to get all gsm and gpl data for the rest of the project
gsm_labels = read.delim(MAIN_GSM_ANNOT_FILE)
all_datasets = unique(as.matrix(gsm_labels[,2:3]))

# Remove duplications 
# - done once, commented later
#sample_counts = table(gsm_labels[,1])
#uniques = names(sample_counts)[sample_counts==1]
#inds = which(is.element(gsm_labels[,1],set=uniques))
#not_uniques = names(sample_counts)[sample_counts>1]
#for(samp in not_uniques){
#	currind = which(gsm_labels[,1]==samp)[1]
#	inds = c(inds,currind)
#}
#inds = sort(inds)
#gsm_labels = gsm_labels[inds,]
#write.table(gsm_labels,file=MAIN_GSM_ANNOT_FILE,sep="\t",quote=F,row.names=F,col.names=T)
#

# Step 1: download all GPL data
gpl_objs = list()
all_gpls = as.character(unique(gsm_labels[,"Platform"]))
for (gpl in all_gpls){
	gpl_obj = getGEO(gpl)
	gpl_objs[[gpl]] = Table(gpl_obj) 	
}
sapply(gpl_objs,colnames)
has_entrez = sapply(gpl_objs,function(x){any(grepl(colnames(x),pattern='entrez',ignore.case=T))})
has_gb = sapply(gpl_objs,function(x){any(grepl(colnames(x),pattern='GB_',ignore.case=F))})
has_gene = sapply(gpl_objs,function(x){any(colnames(x)=="GENE")})
table(has_entrez|has_gene)
sapply(gpl_objs,colnames)[!(has_entrez|has_gene)]

# Step 2: for each row in the gpl data table try to map the probe
# to the relevant entrez gene ids.
# Here we use the bioconductor org.Hs.eg.db package to map the following
# ids into entrez genes: genebank, gene symbol, gene name, ensembl gene, and unigene.
# The relevant columns in the gpl matrix are those that contain the information
# of the ids above or a direct mapping to entrez genes.
# For each row we go over all relevant columns, map them to entrez ids
# and keep the union.
converters = list()
# note that this genebank mapping contains all refseq ids
converters[["GB_ACC"]] = as.list(org.Hs.egACCNUM2EG)
converters[["GB_LIST"]] = converters[["GB_ACC"]]
converters[["symbol"]] = as.list(org.Hs.egSYMBOL2EG)
converters[["gene"]] = converters[["symbol"]]
converters[["ensembl"]] = as.list(org.Hs.egENSEMBL2EG)
converters[["unigene"]] = as.list(org.Hs.egUNIGENE2EG)
# check thy self
# strsplit(c("A","A,B","A /// B"),split="( /// )|,",perl=T)
split_gpl_table_entry <-function(x){strsplit(x,split="( /// )|,|;",perl=T)}
# split_gpl_table_entry(c("A","A,B","A /// B"))
# merges into l1 and returns it as the output
merge_two_id_lists<-function(l1,l2){
	lnew = lapply(1:length(l1),function(x,y,z)union(y[[x]],z[[x]]),y=l1,z=l2)
	return(lnew)
}
map_gpl_rows_to_entrez<-function(x,converters){
	row2entrez=lapply(1:nrow(x),function(x)c())
	# look for entrez column
	is_entrez = which(grepl(colnames(x),pattern="entrez",ignore.case=T))
	if(length(is_entrez)>0){
		curr_entrez = as.character(x[,is_entrez[1]])
		row2entrez = split_gpl_table_entry(curr_entrez)
	}
	for(conv_fields in names(converters)){
		curr_cols = which(grepl(colnames(x),pattern=conv_fields,ignore.case=T))
		print(conv_fields)
		if(length(curr_cols)>0){
			for (j in curr_cols){
				currv = as.character(x[,j])
				currl = split_gpl_table_entry(currv)
				currl_vals = unique(unlist(currl))
				curr_conv = converters[[conv_fields]][currl_vals]
				curr_conv = curr_conv[!is.na(names(curr_conv))]
				currl = lapply(currl,function(x,y)unique(unlist(y[x])),y=curr_conv)
				row2entrez = merge_two_id_lists (row2entrez,currl)
			}
		}
	}
	return(row2entrez)
}

# correct the version issue in the GB_ACC columns
for (gpl in all_gpls){
	curr_col = which(colnames(gpl_objs[[gpl]])=="GB_ACC")
	if(length(curr_col)>0){
		gpl_objs[[gpl]][,"GB_OLD_ACC"] = gpl_objs[[gpl]][,"GB_ACC"]
		gpl_objs[[gpl]][,"GB_ACC"] = gsub(gpl_objs[[gpl]][,"GB_ACC"],pattern="\\.\\d+$",perl=T,replace="")
	}
}
gpl_mappings_to_entrez = lapply(gpl_objs,map_gpl_rows_to_entrez,converters=converters)

# Step 2.2: borrow data across the illumina platforms using the genebank data
# We observed that in the illumina platforms some rows are mapped only using the 
# genbank ids. On the other hand, these are not consistent across the platforms.
# We therefore obtain a unified mapping of genebank ids to entrez
# and then use it to fill gaps.
# Practically it added a few hundrend annotations in each GPL
manufacts = c()
for (gpl in all_gpls){
	gpl_obj = getGEO(gpl)
	manufacts[gpl] = Meta(gpl_obj)$manufacturer
}
illu_gpls = names(manufacts[grepl("Illumin",manufacts)])
all_gbs_to_entrez = list()
for(gpl in illu_gpls){
	curr_gbs = gpl_objs[[gpl]][,"GB_ACC"]
	curr_ez = gpl_mappings_to_entrez[[gpl]]
	for (j in 1:length(curr_gbs)){
		all_gbs_to_entrez[[curr_gbs[j]]] = union(all_gbs_to_entrez[[curr_gbs[j]]],curr_ez[[j]])
	}
}
table(sapply(all_gbs_to_entrez,length))
count_corrections1 = 0
count_corrections2 = 0
for(gpl in illu_gpls){
	curr_gbs = gpl_objs[[gpl]][,"GB_ACC"]
	curr_ez = gpl_mappings_to_entrez[[gpl]]
	for (j in 1:length(curr_gbs)){
		new_ezs = union(curr_ez[[j]],all_gbs_to_entrez[[curr_gbs[j]]])
		new_ezs = new_ezs[!is.na(new_ezs)]
		if(length(new_ezs)>length(curr_ez[[j]])){
			count_corrections1 = count_corrections1+1
			if(length(curr_ez[[j]])==0){count_corrections2 = count_corrections2+1}
			print(paste(c("#####",curr_ez[[j]]),collapse=';'))
			print(all_gbs_to_entrez[[curr_gbs[j]]])
			curr_ez[[j]] = new_ezs
		}
	}
	gpl_mappings_to_entrez[[gpl]] = curr_ez
}


# Step 3: use the correct rows from the GPL table as the probe ids.
# This is done by downloading a sample GSM profile for each gpl.
# We then go over the gpl table and seek the correct column that fits the ids
# in the names of the rows in the downloaded profile
gsms_and_platforms = as.matrix(gsm_labels[,c(1,3)])
for(gpl in all_gpls){
	curr_gsms = gsms_and_platforms[gsms_and_platforms[,2]==gpl,1]
	sample_gsm = curr_gsms[1]
	gsm_obj = Table(getGEO(sample_gsm))
	v = as.numeric(gsm_obj[,"VALUE"])
	names(v) = gsm_obj[,1]
	curr_table = gpl_objs[[gpl]]
	num_appearing = c()
	num_missing = c()
	for(j in 1:ncol(curr_table)){
		num_appearing[j] = length(intersect(names(v),curr_table[,j]))
		num_missing[j] = length(setdiff(names(v),curr_table[,j]))
	}
	curr_col = which(num_appearing==max(num_appearing))[1]
	print(num_missing)
	#if(length(curr_cols)>1){
	#	problematic_rows = which(!apply(curr_table[,curr_cols],1,function(x)all(x==x[1])))
	#	curr_table[problematic_rows,curr_cols]
	#}
	names(gpl_mappings_to_entrez[[gpl]]) = curr_table[,curr_col]
}

# Step 4: we now have a mapping of probes to entrez genes,
# get the reverse mapping
reverse_mapping_list<-function(l){
	newl = list()
	for(currname in names(l)){
		currvals = l[[currname]]
		for(v in currvals){
			newl[[v]] = union(newl[[v]],currname)
		}
	}
	return(newl)
}
gpl_mappings_entrez2probes = list()
for(gpl in all_gpls){
	gpl_mappings_entrez2probes[[gpl]] = reverse_mapping_list(gpl_mappings_to_entrez[[gpl]])
}
sapply(gpl_mappings_entrez2probes,length)
save(gpl_mappings_to_entrez,gpl_mappings_entrez2probes,file='gpl_mappings_to_entrez.RData')

################ QA: compare to the old mapping from Adeptus1 ##############
#load('gpl_mappings_to_entrez.RData')
#library(hash);library(hash,lib.loc='rlibs/')
#adeptus1_gpls = get(load(url('http://acgt.cs.tau.ac.il/adeptus/data/platform_entrez2probes_full.RData')))
#comparison_stats = c()
#counter=0
#for(gpl in names(adeptus1_gpls)){
#	adeptus1_entrez2probes = adeptus1_gpls[[gpl]]
#	if(!is.element(gpl,set=names(gpl_mappings_entrez2probes))){next}
#	adeptus2_entrez2probes = gpl_mappings_entrez2probes[[gpl]]
#	v = c()
#	v["numgenes_old"] = length(adeptus1_entrez2probes)
#	v["numgenes_new"] = length(adeptus2_entrez2probes)
#	in_both = intersect(names(adeptus1_entrez2probes),names(adeptus2_entrez2probes))
#	v["intersect_length"] = length(in_both)
#	same_annots = sapply(in_both,function(x,y,z)length(intersect(y[[x]],z[[x]]))==length(y[[x]]),y=adeptus1_entrez2probes,z=adeptus2_entrez2probes)
#	print(table(same_annots))
#	v["genes_with_same_annots"] = sum(same_annots)
#	comparison_stats = rbind(comparison_stats,v)
#	counter=counter+1
#	rownames(comparison_stats)[counter] = gpl
# 	}
#############################################################################
################# Test which platforms to remove ############################
# This code was used to check which platforms to remove before the analysis
# Here, we wanted to exclude platforms with a low number of samples, whose 
# removal will result in a larger gene set for the analysis.
# The current distributed annotation files do not contain the platforms
# we removed, therefore this code was commented out.
#gpls_to_rem = c("GPL5356")
#gpl2numsamples = sort(table(gsm_labels[,3]))
#all_gpls = unique(gsm_labels[,3])
#get_names_intersect<-function(l){
#	inter = names(l[[1]])
#	for(j in 2:length(l)){
#		inter = intersect(inter,names(l[[j]]))
#	}
#	return(inter)
#}
#without_removal = length(get_names_intersect(gpl_mappings_entrez2probes))
#for(j in 1:length(gpl2numsamples)){
#	curr_to_rem = names(gpl2numsamples)[1:j]
#	remaining_gpls = setdiff(all_gpls,curr_to_rem)
#	print(length(get_names_intersect(gpl_mappings_entrez2probes[remaining_gpls])))
#}
#gpls_to_rem = names(gpl2numsamples)[1:5]
#gsm_lines_to_remove = is.element(gsm_labels[,3],set=gpls_to_rem)
#newgsm_table = gsm_labels[!gsm_lines_to_remove,]
#write.table(newgsm_table,file='all_gsm_labels_cleaned.txt',sep="\t",col.names=T,row.names=F,quote=F)
#############################################################################
################ A framework to get the GSM matrix of a platform ############
get_gsm_profile<-function(gsm){
	try({
		gsm_obj = Table(getGEO(gsm))
		v = as.numeric(gsm_obj[,"VALUE"])
		names(v) = gsm_obj[,1]
		return(v)
	})
	return (NULL)
}
get_gsm_matrix<-function(gsms,mc.cores=1){
		gsms_list = mclapply(gsms,get_gsm_profile,mc.cores=mc.cores)
		names(gsms_list) = gsms
		# remove failed downloads
		to_rem = sapply(gsms_list,is.null)
		failed_gsms = names(which(to_rem))
		gsms_list = gsms_list[!to_rem]
		curr_rows = names(gsms_list[[1]])
		# tests: are the rows in the gsms in the same order
		same_order = sapply(gsms_list,function(x,y){all(names(x)==y)},y=curr_rows)
		print ("table of whether the downloaded gsm profiles have the same row order:")
		print(table(same_order))
		# if there are order "erros" this code also reorders
		gsms_mat = sapply(gsms_list,function(x,y)x[y],y=curr_rows)
		rownames(gsms_mat) = curr_rows
		colnames(gsms_mat) = names(gsms_list)
		return(list(gsms_mat=gsms_mat,failed_gsms=failed_gsms))
}
# This function assumes that all gses are of the same GPL
get_gsm_matrix_from_gses<-function(gses,gpl){
	gsms_mat = c()
	for (gse in gses){
		try({
			curr_mat = getGEO(gse,GSEMatrix=T,getGPL=F)
			curr_mat = exprs(curr_mat[[which(sapply(curr_mat,annotation)==gpl)[1]]])
			if(length(gsms_mat)==0){
				gsms_mat = curr_mat
			}
			else{
				gsms_mat = cbind(gsms_mat,curr_mat[rownames(gsms_mat),])
			}
		})
	}
	return(gsms_mat)
}

################ Tests: compare to obtaining the GSE matrices ############
#testgse = 'GSE46602'
#testgse = 'GSE32225'
#testgse = 'GSE32665'
#currgsms = newgsm_table[newgsm_table[,2]==testgse,1]
#gse_mat = getGEO(testgse,GSEMatrix=T,getGPL=F)
#gse_mat = exprs(gse_mat[[1]])
#gsms_mat = get_gsm_matrix(currgsms,mc.cores=10)
#length(intersect(colnames(gsms_mat),colnames(gse_mat)))==ncol(gsms_mat)
#length(intersect(colnames(gsms_mat),colnames(gse_mat)))==ncol(gse_mat)
#length(intersect(rownames(gsms_mat),rownames(gse_mat)))==nrow(gsms_mat)
#length(intersect(rownames(gsms_mat),rownames(gse_mat)))==nrow(gse_mat)
#gse_mat = gse_mat[rownames(gsms_mat),colnames(gsms_mat)]
# The following correlations should all be 1
#for(j in 1:ncol(gsms_mat)){
#	print(cor(gsms_mat[,j],gse_mat[,j],method='spearman'))
#}
#############################################################################
load('gpl_mappings_to_entrez.RData')
setwd('Adeptus2/data/expression_data/GEO')
gsm_labels = read.delim('sample_annotation.txt')
all_gpls = as.character(unique(gsm_labels[,3]))
for(gpl in all_gpls){
	# This internal loop is needed because the GSM download process might 
	# fail when a large list of GSMs is given.
	download_completed = F
	num_iterations = 0
	while(!download_completed && num_iterations < 5){
		print (gpl)
		curr_file = paste(gpl,'_expression_profiles.RData',sep='')
		curr_gsms = unique(as.character(gsm_labels[gsm_labels[,3]==gpl,1]))
		previous_gsm_x = c()
		# If we already have a matrix for this GPL, check if we
		# have samples with failed download
		if(is.element(curr_file,set=list.files('.'))){
			previous_gsm_x = get(load(curr_file))
			# re-download samples whose profiles are all NAs
			curr_gsms = setdiff(curr_gsms,colnames(previous_gsm_x))
			curr_gsms = colnames(previous_gsm_x)[apply(is.na(previous_gsm_x),2,all)]
			gc()
			print (paste("Previous matrix was founded, number of missing GSMs:",length(curr_gsms)))
			if(length(curr_gsms)==0){download_completed = T; next}
			previous_gsm_x = previous_gsm_x[,setdiff(colnames(previous_gsm_x),curr_gsms)]
		}
		print(paste(gpl,length(curr_gsms)))
		gsms_profiles = get_gsm_matrix(curr_gsms,mc.cores=5)
		gsms_mat = gsms_profiles[[1]]
		failed_gsms = gsms_profiles[[2]]
		if(length(previous_gsm_x)>0){gsms_mat = cbind(gsms_mat,previous_gsm_x)}
		save(gsms_mat,file=paste(gpl,'_expression_profiles.RData',sep=''))
		if(length(failed_gsms)>0){save(failed_gsms,file=paste(gpl,'_failed_gsms.RData',sep=''))}
		print(dim(gsms_mat))
		rm(gsms_mat);rm(previous_gsm_x);gc()
		num_iterations = num_iterations+1
	}
}

setwd('Adeptus2/data/expression_data/GSE60862')
gsm_labels = read.delim('sample_annotation.txt')
all_gpls = as.character(unique(gsm_labels[,3]))
for(gpl in all_gpls){
	# This internal loop is needed because the GSM download process might 
	# fail when a large list of GSMs is given.
	download_completed = F
	num_iterations = 0
	while(!download_completed && num_iterations < 5){
		print (gpl)
		curr_file = paste(gpl,'_expression_profiles.RData',sep='')
		curr_gsms = unique(as.character(gsm_labels[gsm_labels[,3]==gpl,1]))
		previous_gsm_x = c()
		# If we already have a matrix for this GPL, check if we
		# have samples with failed download
		if(is.element(curr_file,set=list.files('.'))){
			previous_gsm_x = get(load(curr_file))
			# re-download samples whose profiles are all NAs
			curr_gsms = setdiff(curr_gsms,colnames(previous_gsm_x))
			curr_gsms = colnames(previous_gsm_x)[apply(is.na(previous_gsm_x),2,all)]
			gc()
			print (paste("Previous matrix was founded, number of missing GSMs:",length(curr_gsms)))
			if(length(curr_gsms)==0){download_completed = T; next}
			previous_gsm_x = previous_gsm_x[,setdiff(colnames(previous_gsm_x),curr_gsms)]
		}
		print(paste(gpl,length(curr_gsms)))
		gsms_profiles = get_gsm_matrix(curr_gsms,mc.cores=5)
		gsms_mat = gsms_profiles[[1]]
		failed_gsms = gsms_profiles[[2]]
		if(length(previous_gsm_x)>0){gsms_mat = cbind(gsms_mat,previous_gsm_x)}
		save(gsms_mat,file=paste(gpl,'_expression_profiles.RData',sep=''))
		if(length(failed_gsms)>0){save(failed_gsms,file=paste(gpl,'_failed_gsms.RData',sep=''))}
		print(dim(gsms_mat))
		rm(gsms_mat);rm(previous_gsm_x);gc()
		num_iterations = num_iterations+1
	}
}

# OPTIONAL: for the failed cases try to get the gse matrix
#for(gpl in all_gpls){
#	# if file exists in records then move on
#	curr_expr_file = paste(gpl,'_expression_profiles.RData',sep='')
#	try({
#		curr_failed_gsms_file = paste(gpl,'_failed_gsms.RData',sep='')
#		load(curr_failed_gsms_file)
#		#print(paste(gpl,length(failed_gsms)))
#		if(length(failed_gsms)==0){next}
#		curr_gses = unique(as.character(gsm_labels[is.element(gsm_labels[,1],set=failed_gsms),2]))
#		print(unique(curr_gses))
#		gse_mat = get_gsm_matrix_from_gses(curr_gses,gpl)
#		load(curr_expr_file)
#		gse_mat = gse_mat[,failed_gsms]
#	})
#}












