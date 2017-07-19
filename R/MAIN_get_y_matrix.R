setwd('Z:/davidama/Adeptus2/data/expression_data/')
setwd('C:/Users/davidama/Desktop/Adeptus2/data/expression_data/')

library(DO.db)
dos = as.list(DOTERM)
ancs = as.list(DOANCESTOR)

tissue_slim_mapping_file = '../sample_labels/new_tissue_slim_mapping.txt'
tissue_slim_mapping = as.matrix(read.delim(tissue_slim_mapping_file,sep='\t',row.names=1))

dirs = c('GEO','GSE60862','gtex','tcga')
term2name<-function(n){
	if (is.null(DOTERM[[n]])){term=n}else{term = Term(DOTERM[[n]])}
	return (term)
}

###################
MIN_NUM_DATASETS = 3
MIN_NUM_SAMPLES = 100
###################

sample2terms = list()
sample2study = c()
sample2tissue_slim = c()
tissue_column_ind = 6
doid_column_ind = 9
for(dir in dirs){
	curr_annots = as.matrix(read.delim(paste(dir,"/sample_annotation.txt",sep='')))
	# correct the tissue names to lowercase
	curr_annots[,tissue_column_ind] = tolower(curr_annots[,tissue_column_ind])
	curr_annots[,tissue_column_ind] = tissue_slim_mapping[curr_annots[,tissue_column_ind],1]
	sample2tissue_slim[curr_annots[,1]] = curr_annots[,tissue_column_ind]
		
	# add the tissue slim of each sample
	curr_tissues = curr_annots[,tissue_column_ind]

	# add all control samples
	curr_control_samples = curr_annots[curr_annots[,doid_column_ind]=="",1]
	sample2terms[curr_control_samples] = paste('control',curr_tissues[curr_annots[,doid_column_ind]==""],sep=';')
	
	# add the study of each sample
	sample2study[curr_annots[,1]] = curr_annots[,2]
	curr_dos = unique(curr_annots[,c(doid_column_ind,tissue_column_ind)])
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
		doids = union(doids,unique((unname(unlist(ancs[doids])))))
		# the name of a class: doid;disease_name;tissue_slim
		full_class_names = paste(doids,sapply(doids,term2name),sep=';')
		full_class_names = paste(full_class_names,tissue,sep=';')
		for(s in curr_samples){
			sample2terms[[s]] = full_class_names
		}
	}
}
table(sample2study)
table(sample2tissue_slim)
length(sample2terms)

all_classes = unique(unlist(sample2terms))
all_studies = unique(sample2study)
class2num_datasets = rep(0,length(all_classes));names(class2num_datasets) = all_classes
class2num_samples = rep(0,length(all_classes));names(class2num_samples) = all_classes
all_doids = unique(sapply(all_classes,function(x)strsplit(x,split=';')[[1]][1]))
all_doids = all_doids[grepl("DOID",all_doids)]
doid2num_datasets = rep(0,length(all_doids));names(doid2num_datasets) = all_doids
min_size=10
for(study in all_studies){
	curr_samps = names(which(sample2study==study))
	curr_classes_counts = table(unlist(sample2terms[curr_samps]))
	curr_classes_counts = curr_classes_counts[curr_classes_counts>=min_size]
	if(length(curr_classes_counts)==0){next}
	class2num_datasets[names(curr_classes_counts)] = class2num_datasets[names(curr_classes_counts)]+1
	class2num_samples[names(curr_classes_counts)] = class2num_samples[names(curr_classes_counts)] + curr_classes_counts
	curr_classes = names(curr_classes_counts)
	curr_doids = unique(sapply(curr_classes,function(x)strsplit(x,split=';')[[1]][1]))
	curr_doids = intersect(all_doids,curr_doids)
	doid2num_datasets[curr_doids] = doid2num_datasets[curr_doids]+1
}
sum(class2num_datasets>=MIN_NUM_DATASETS & class2num_samples>=MIN_NUM_SAMPLES)

terms_for_analysis = names(which(class2num_datasets>=MIN_NUM_DATASETS & class2num_samples>=MIN_NUM_SAMPLES))
doids_for_analysis = names(doid2num_datasets)[doid2num_datasets>=MIN_NUM_DATASETS]
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
save(y,sample2study,sample2tissue_slim,file='adeptus2_pheno_data.RData')
sort(apply(y,2,sum))
y = y[,colSums(y)>=MIN_NUM_SAMPLES]


# Remove labels that have the same sample sets: keep the ancestor only
#source("http://bioconductor.org/biocLite.R")
#biocLite('DO.db',lib='rlibs')
#library(DO.db,lib.loc='rlibs')
try(library(DO.db))
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
		#print(paste(i,j))
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
save(y,sample2study,sample2tissue_slim,file='adeptus2_pheno_data.RData')








