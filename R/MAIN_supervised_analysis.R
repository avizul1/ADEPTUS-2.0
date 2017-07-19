# CMD example:
# nohup R-3.1.1 --no-save --no-restore --args 100 2 < MAIN_supervised_analysis.R > test_method2_100cv_5000vargenes_200.Rout &

setwd('Adeptus2/data/')

args = commandArgs(T)
if (length(args)!=2){
        print ("Usage:<num folds><method>")
        q("no")
}
num_folds = as.numeric(args[1])
method_num = as.numeric(args[2])

##########################################################################
## installation in unix
#source("http://bioconductor.org/biocLite.R")
#cran_libs = c('hash','e1071','pROC','randomForest','bnlearn','ranger','LiblineaR')
#bioconductor_libs = c('S4Vectors','CMA','preprocessCore','DO.db','limma',"RBGL","BiocGenerics",'gRbase','gRain','IRanges','AnnotationDbi')
#for(b in bioconductor_libs){	
#	try({biocLite(b,lib = paste(getwd(),'/rlibs',sep=''), 
#		lib.loc = paste(getwd(),'/rlibs',sep=''), destdir = paste(getwd(),'/rlibs',sep=''),dependencies=T)})
#}
#for (l in cran_libs){
#	try({install.packages(l,lib = paste(getwd(),'/rlibs',sep=''), 
#		lib.loc = paste(getwd(),'/rlibs',sep=''), destdir = paste(getwd(),'/rlibs',sep=''),dependencies=T)})
#}
#all_libs = c(bioconductor_libs,cran_libs)
#for(j in 1:length(all_libs)){
#	library(all_libs[j],lib.loc = paste(getwd(),'/rlibs',sep=''),character.only=T)
#}
##########################################################################
library('parallel')
cran_libs = c('hash','e1071','pROC','randomForest','bnlearn','ranger','LiblineaR')
bioconductor_libs = c('S4Vectors','CMA','preprocessCore','DO.db','limma',"RBGL","BiocGenerics",'gRbase','gRain')
all_libs = c(bioconductor_libs,cran_libs)
for(j in 1:length(all_libs)){
	try({library(all_libs[j],lib.loc = paste(getwd(),'/rlibs',sep=''),character.only=T)})
}
#.libPaths( c( .libPaths(), paste(getwd(),"/rlibs",sep='')) )
#for (p in (c(bioconductor_libs,cran_libs))){try({library(p,character.only=T)})}
##########################################################################
print ("##################################################")
print ("Info of this running session:")
print(sessionInfo())

##########################################################################
NUM_FEATURES = 200
NUM_VAR_GENES = 5000
USE_TOP_VAR_GENES = T
MAX_NUM_THREADS = 40
SHUFFLE_SAMPLE_LABELS = F
CLASSIFICATION_RESULTS_FILE = 'classification_results_40k_samples.RData'
DATA_FILE = 'Supervised_data_atleast_10000_genes_max_percet_na_10.RData'
USE_GTEX = T
##########################################################################

if(!any(grepl(pattern = CLASSIFICATION_RESULTS_FILE,list.files()))){
	classification_results = list()
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

load(DATA_FILE)
all_samples = intersect(rownames(x),rownames(y))
if(!USE_GTEX){
	gtex_samples = grepl(all_samples,pattern='gtex',ignore.case=T)
	all_samples = all_samples[!gtex_samples]
}
x = x[all_samples,];y=y[all_samples,]
sample2study = sample2study[all_samples]

## QA: check that all labels have at least 100 samples from at least 4 studies
#class2num_studies = c()
#class2num_samples = c()
#for(j in 1:ncol(y)){
#	curr_samples = y[,j]!=0
#	class2num_studies[j] = length(unique(sample2study[curr_samples]))
#}
#sort(class2num_studies)
#sort(colSums(y))[1:5]

source("../rscripts/FUNCTIONS_multilabel_classification.R")

# for tests
#samp_rows = sample(1:nrow(x))[1:3000]
#samp_labels = sample(1:ncol(y))[1:20]
#x = x[samp_rows,1:500];y=y[samp_rows,samp_labels];gc()
#d = sample2study[rownames(x)]
#d_ldo_vec = getKfoldLeaveDatasetOutVector(d,length(unique(d))/10)

#svm_br_ldo_results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=length(unique(d_ldo_vec)),
#	baseClassifierArgs=list(classification_function=svm,numFeatures=100,probability=T),
#	classification_function = simple_binary_relevance,baseClassifierPredArgs=list(probability=T),
#	baseClassifier=featureSelectionClassifier)

#print(svm_br_ldo_results$class2ROC)
#q("no")
#########################

# For faster results: take the top 1000 variance genes
sds = apply(x,2,sd)
gene_var_thr = sort(sds,decreasing=USE_TOP_VAR_GENES)[NUM_VAR_GENES]
if (USE_TOP_VAR_GENES){
	cols_for_analysis = sds >= gene_var_thr
}
if (!USE_TOP_VAR_GENES){
	cols_for_analysis = sds <= gene_var_thr
}
x = x[,cols_for_analysis];gc()
print ("##################################################")
print("size of x before the cv process starts:")
print(dim(x))
############################################################

## Shuffle the sample labels - QA
if (SHUFFLE_SAMPLE_LABELS){
	samps = rownames(x)
	old_rownames = rownames(x)
	samps = sample(samps)
	x = x[samps,]
	rownames(x) = old_rownames
	print("rows where shuffled in x")
	# for fast results:
	# take the top 20 variance y vectors
	# select 10 features, not 100
	# use 10 cv, not 20
	#y_vars = apply(y,2,var)
	#y_thr = sort(y_vars,decreasing=T)[20]
	#y = y[,y_vars>=y_thr]
	#num_folds = 10
}
############################

getKfoldLeaveDatasetOutVector<-function(d,K,tosample=T,seed=10){
	Nd = length(d);	Nds = length(unique(d))
	data2fold = c(sapply(1:K,rep,(Nds/K + 1)))[1:Nds]
	if(!is.null(seed)){set.seed(seed)}
	if (tosample){data2fold = sample (data2fold)}
	names(data2fold)<-unique(d)
	newd = data2fold[d]
	names(newd)=names(d)
	return (newd)
}

d = sample2study[rownames(x)]
num_folds = min(num_folds,length(unique(d)))
d_ldo_vec = getKfoldLeaveDatasetOutVector(d,num_folds)
print ("##################################################")
print("sizes of the folds:")
print(table(d_ldo_vec))
num_cores = min(MAX_NUM_THREADS,num_folds)

print ("##################################################")
print ("Number of labels:")
print(ncol(y))

print ("##################################################")
print ("Starting the cross validation")

# Standard SVM implementation
if (method_num==1){
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=svm,numFeatures=NUM_FEATURES,probability=T),
		classification_function = simple_binary_relevance,baseClassifierPredArgs=list(probability=T),
		baseClassifier=featureSelectionClassifier)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[['svm_default_100']] = results
	classification_results[['svm_default_100']][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

# Fast SVM (libLineaR) + sampling
if (method_num==2){
	method_name = paste('liblinear_svm_subsampling_2_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=100,max_num_samples = 500),
		classification_function = simple_binary_relevance,baseClassifierPredArgs=list(),
		baseClassifier=simple_sampling_based_learner)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

# Fast SVM (libLineaR) + sampling + BC, no CV
if (method_num==3){
	method_name = paste('liblinear_svm_subsampling_BCnoCV_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=5),
		classification_function = multitaskChainClassifier,baseClassifierPredArgs=list(),
		baseClassifier=simple_sampling_based_learner,
		secondClassifierPredArgs = list(),secondClassifier = bayesianCorrectionClassifier_new,
		secondClassifierArgs=list(),isSecondClassifierBinaryRelevance=F,cvfolds=1
	)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

if (method_num==4){
	method_name = paste('liblinear_svm_subsampling_BC2CV_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=10),
		classification_function = multitaskChainClassifier,baseClassifierPredArgs=list(),
		baseClassifier=simple_sampling_based_learner,
		secondClassifierPredArgs = list(),secondClassifier = bayesianCorrectionClassifier_new,
		secondClassifierArgs=list(),isSecondClassifierBinaryRelevance=F,cvfolds=2
	)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

if (method_num==5){
	method_name = paste('liblinear_svm_subsampling_viamulticlass_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		classification_function = multitaskViaMulticlassClassifier,
		baseClassifierArgs=list(baseClassifier=liblinear_binary_svm_wrapper,baseClassifierArgs=list(),baseClassifierPredArgs=list()),useFS=F)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

if (method_num==6){
	method_name = paste('liblinear_svm_subsampling_viamulticlassRF_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		classification_function = multitaskViaMulticlassClassifier,
		baseClassifierArgs=list(baseClassifier=liblinear_binary_svm_wrapper,baseClassifierArgs=list(ntree=5,nodesize=50)),useFS=F)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

if (method_num==7){
	method_name = paste('chain_noCV_liblinear_svm_subsampling_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=10),
		classification_function = multitaskChainClassifier,baseClassifierPredArgs=list(),
		baseClassifier=simple_sampling_based_learner,
		secondClassifierPredArgs = list(),secondClassifier = simple_sampling_based_learner,
		secondClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=50,d=sample2study,reps=10),
		isSecondClassifierBinaryRelevance=T,cvfolds=1
	)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

if (method_num==8){
	method_name = paste('chain_2CV_liblinear_svm_subsampling_',NUM_FEATURES,'_',num_folds,"cv_",NUM_VAR_GENES,"topVarGenes",sep="")
	currtime1 = Sys.time()
	results = getMultitaskLeaveDatasetOutClassificationPerformance_multicore(t(x),y,d_ldo_vec,numCores=num_cores,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=10),
		classification_function = multitaskChainClassifier,baseClassifierPredArgs=list(),
		baseClassifier=simple_sampling_based_learner,
		secondClassifierPredArgs = list(),secondClassifier = simple_sampling_based_learner,
		secondClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=50,d=sample2study,reps=10),
		isSecondClassifierBinaryRelevance=T,cvfolds=2
	)
	currtime2 = Sys.time()
	results[['running time']] = list(start=currtime1,end=currtime2)
	load(CLASSIFICATION_RESULTS_FILE)
	classification_results[[method_name]] = results
	classification_results[[method_name]][["method_number"]] = method_num
	save(classification_results,file=CLASSIFICATION_RESULTS_FILE)
}

# Tests from BC_qa
#baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,numFeatures=NUM_FEATURES,d=sample2study,reps=2)
#baseClassifiers = apply(y,2,getBaseClassifier,x=x,func=simple_sampling_based_learner,baseClassifierArgs)
#basePredictions = sapply(baseClassifiers,getPredictionVector,x=x)
#BC = bayesianCorrectionClassifier_new(basePredictions,y)
#inds=1:20000
#pred = predict(BC,basePredictions[inds,])
#table((pred>0.5) == (y[inds,]>0.5))
#table((basePredictions[inds,]>0.5) == (y[inds,]>0.5))

print ("##################################################")
print ("cross validation scores")
# Measures
P = results$predictions
Y = results$y
D = sample2study[rownames(Y)]
print(paste("rownames of P and Y are the same:",all(rownames(P)==rownames(Y))))

auc_scores = getMultitaskPerformancePositivesVsNegatives(P,Y,D)
auc_scores = sapply(auc_scores,function(x)x)

print ("########### PB- and PN- AUC scores ####################")
print ("classes with PN-ROC>0.75 and PB-ROC>0.75")
v = auc_scores[,1]>0.75 & auc_scores[,3]>0.75
print(table(v))

# Study-based analysis
get_meta_analysis_eval<-function(p,y,d,minP=10,minN=10){
	qvals = c();all_pvals=c()
	for (j in 1:ncol(y)){
		P = rownames(y)[which(y[,j]==1)]
		DS = unique(d[P])
		N = rownames(y)[which(y[,j]==0 & is.element(d[rownames(y)],set=DS))]
		zs = c();ws=c();ps=c()
		for (ds in DS){
			#print (ds)
			P1 = P[d[P]==ds]
			N1 = N[d[N]==ds]
			#print (paste(length(P1),length(N1)))
			if (length(P1)<minP || length(N1) < minN){next}
			pval = wilcox.test(p[P1,j],p[N1,j],alternative="greater",exact=F)$p.value
			z = qnorm(1-pval)
			zs[ds]=z
			ws[ds] = length(P1)+length(N1)
			ps[ds] = pval
		}
		if (length(zs) < 1){next}
		all_pvals[[colnames(y)[j]]] = ps
		zs[zs<0 & is.infinite(zs)] = -10
		zs[zs>0 & is.infinite(zs)] = 10
		qvals[colnames(y)[j]] = getStouferPval(zs)
	}
	qvals = p.adjust(qvals,method='fdr')
	sum(sapply(all_pvals,function(x)sum(x<0.1)/length(x)) > 0.5)
	return (list(qvals=qvals,all_pvals=all_pvals))
}

get_meta_analysis_other_stats<-function(p,y,d,minP=10,minN=10,donames=F){
	stats = c()
	for (j in 1:ncol(y)){
		P = rownames(y)[which(y[,j]==1)]
		DS = unique(d[P])
		N = rownames(y)[which(y[,j]==0 & is.element(d[rownames(y)],set=DS))]
		num_dss = 0
		num_pss = 0
		for (ds in DS){
			P1 = P[d[P]==ds]
			N1 = N[d[N]==ds]
			#print (paste(length(P1),length(N1)))
			if (length(P1)<minP || length(N1) < minN){next}
			num_dss = num_dss+1
			num_pss = num_pss+length(P1)
		}
		stats = rbind(stats,c(num_dss,num_pss,length(P)))
	}
	colnames(stats) = c("Num datasets for q analysis", "Num positives covered","total positives")
	rownames(stats) = colnames(y)
	if (donames){rownames(stats) = sapply(colnames(y),term2name)}
	return (stats)
}

meta_sig_tests = get_meta_analysis_eval(P,Y,D)
meta_analysis_stats = get_meta_analysis_other_stats(P,Y,D)
print ("########### meta analysis results ####################")
print ("classes with a significant q-value, 0.001")
print (table(meta_sig_tests$qvals<0.001))
print ("number of datasets with a p-value <0.01 in at least half of the labels")
print (table(sapply(meta_sig_tests$all_pvals,function(x)sum(x<0.01)/length(x))>0.5))


















