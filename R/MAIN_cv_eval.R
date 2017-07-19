# nohup R --no-save --no-restore < MAIN_cv_eval.R > test_cv_eval.Rout &
try(setwd('/Adeptus2/data/'))
try(setwd('/home/gaga/davidama/Adeptus2/data/'))
try(setwd('C:/Users/davidama/Desktop/Adeptus2/data/'))

##########################################################################
library('parallel')
cran_libs = c('hash','e1071','pROC','randomForest','bnlearn','ranger','LiblineaR')
bioconductor_libs = c('CMA','preprocessCore','DO.db','limma')
#for (p in cran_libs){try({install.packages(p,lib='rlibs')})}
#source("http://bioconductor.org/biocLite.R")
#for (p in bioconductor_libs){try({biocLite(p,lib='rlibs')})}
.libPaths( c( .libPaths(), paste(getwd(),"/rlibs",sep='')) )
for (p in (c(bioconductor_libs,cran_libs))){try({library(p,character.only=T)})}

library(PRROC); library(ROCR)
##########################################################################

##########################################################################
CLASSIFICATION_RESULTS_FILE = 'classification_results_40k_samples.RData'
DATA_FILE = 'Supervised_data_atleast_10000_genes_max_percet_na_10.RData'
##########################################################################

source("../rscripts/FUNCTIONS_multilabel_classification.R")

# Implementation for the rep-analysis
##########################################################################
# Study-based analysis
# This function receives as input a vector y of labels, an assignment of
# samples to their datasets, and a vector of scores p.
# We go over all datasets and collect the following statistics in a matrix
#	A row for each dataset
#	The p-value, z-score, and the number of samples
get_study_based_evaluation<-function(p,y,d,minP=10,minN=10){
	P = names(y)[y==1]
	DS = unique(d[P])
	N = names(y)[y==0 & is.element(d[names(y)],set=DS)]
	zs = c();ps=c();num_pos=c();num_neg=c();rocs=c();auprs=c();ps2=c();ws = c()
	for (ds in DS){
		P1 = P[d[P]==ds];N1 = N[d[N]==ds]
		if (length(P1)<minP || length(N1) < minN){next}
		pval = wilcox.test(p[P1],p[N1],alternative="greater",exact=F)$p.value
		pval2 = wilcox.test(p[P1],p[N1],alternative="two.sided",exact=F)$p.value
		z = qnorm(1-pval)
		zs[ds]=z;ps[ds] = pval;ps2[ds]=pval2
		num_pos[ds] = length(P1)
		num_neg[ds] = length(N1)
		ws[ds] = length(P1) + length(N1)
		curr_scores = c(p[P1],p[N1])
		curr_lab = c(rep(1,length(P1)),rep(0,length(N1)))
		inds = sample(1:length(curr_scores))
		curr_scores = curr_scores[inds];curr_lab=curr_lab[inds]
		rocs[ds] = calcAupr(curr_scores,curr_lab,roc=T)
		auprs[ds] = calcAupr(curr_scores,curr_lab,roc=F)
	}
	if (length(zs) < 1){return(NULL)}
	m = cbind(ps,ps2,zs,rocs,auprs,num_pos,num_neg,ws)
	return(m)
}
get_study_weights<-function(y,d,minP=10,minN=10){
	P = names(y)[y==1]
	DS = unique(d[P])
	N = names(y)[y==0 & is.element(d[names(y)],set=DS)]
	ws = c()
	for (ds in DS){
		P1 = P[d[P]==ds];N1 = N[d[N]==ds]
		if (length(P1)<minP || length(N1) < minN){next}
		ws[ds] = length(P1) + length(N1)
	}
	return(ws)
}
get_pval_zscore<-function(p){
	if (p<0.5){return(qnorm(p))}
	return (qnorm(p,lower.tail=T))
}
get_matrices_study_based_evaluation<-function(p,y,d,...){
	l = list()
	for (j in 1:ncol(y)){
		currname = colnames(y)[j]
		l[[currname]] = get_study_based_evaluation(P[,j],Y[,j],D,...)
	}
	return (l)
}

# Implementation for the meta-analysis
##########################################################################
# 1. Auxiliary methods for analysis of edge lists
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
#get_siblings("DOID:4;disease;brain_star","DOID:4;disease;brain",asgraph)
#Y[,"DOID:4;disease;brain_star"]
get_all_siblings<-function(node,edges){
	sibs = c()
	pars = get_parents(node,edges)
	for(p in pars){
		sibs = union(sibs,get_siblings(node,p,edges))
	}
	return(unique(sibs))
}

# 2. Auxiliary method to add edges and labels for the analysis
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
		S1_aux = rownames(Y)[Y[,children]>0]
		S1 = setdiff(S1,S1_aux)
		if(length(S1)<minSize){next}
		v = rep(0,nrow(Y));names(v) = rownames(Y)
		v[S1] = 1
		m = cbind(m,v)
		newnames = c(newnames,paste(node,"_star",sep=""))
		edges = rbind(edges,c(node,paste(node,"_star",sep="")))
	}
	rownames(m) = rownames(Y);colnames(m) = newnames
	return(list(edges=edges,m=m))
}

# Tests
#DOnet = getBayesianNetwork(colnames(Y))
#mapping = DOnet[[2]]
#asgraph = get_graph_for_cytoscape(arcs(DOnet[[1]]),mapping)
#get_parents("DOID:1612;breast cancer;breast",asgraph)
#get_children("DOID:1612;breast cancer;breast",asgraph)
#get_siblings("DOID:1612;breast cancer;breast","breast",asgraph)
#get_siblings("DOID:7;disease of anatomical entity;brain","DOID:7",asgraph)

# Network-based analysis: edge-based analysis
# This function receives as input a matrix of labels Y, a network N, and a score vector p
# For each edge in the network N (in the format specified in the auxiliary functions above)
# We define two sets of samples:
#	S1: the samples of the son term
#	S2: the samples of the siblings union samples of the parent minus S1
#	We then compare the scores in p for the given groups
evaluate_an_edge_vs_scores<-function(parent,son,Y,N,p,filter_siblings_by="DOID",minSize=50){
	S1 = rownames(Y)[Y[,son]==1]
	S2 = c()
	siblings = intersect(colnames(Y),get_siblings(son,parent,N))
	siblings = setdiff(siblings,son)
	if(!is.null(filter_siblings_by)){
		siblings=siblings[grepl(siblings,pattern=filter_siblings_by)]
	}
	if(length(siblings)>1){
		S2_aux = rownames(Y)[rowSums(Y[,siblings])>0]
		S2 = union(S2,S2_aux)
	}
	if(length(siblings)==1){
		S2_aux = rownames(Y)[Y[,siblings]>0]
		S2 = union(S2,S2_aux)
	}
	S3 = rownames(Y)[Y[,parent]>0]
	S3 = setdiff(S3,union(S1,S2))
	print(length(S3))
	inters = intersect(S2,S1)
	S2 = setdiff(S2,inters)
	S2 = union(S2,S3)
	S1 = setdiff(S1,inters)
	curr_scores = c(p[S1],p[S2])
	curr_lab = c(rep(1,length(S1)),rep(0,length(S2)))
	inds = sample(1:length(curr_scores))
	curr_scores = curr_scores[inds];curr_lab=curr_lab[inds]
	roc = calcAupr(curr_scores,curr_lab,roc=T)
	aupr = calcAupr(curr_scores,curr_lab,roc=F)
	if(length(S1)<minSize||length(S2)<minSize){
		roc=NA;aupr=NA
	}
	return(c(roc=roc,aupr=aupr,num_S1 = length(S1),num_S2 = length(S2)))
}
evaluate_an_edge_vs_scores_TEST<-function(parent,son,Y,N,p,filter_siblings_by="DOID",minSize=50){
	S1 = rownames(Y)[Y[,son]==1]
	S2 = rownames(Y)[Y[,parent]>0]
	S2 = setdiff(S2,S1)
	curr_scores = c(p[S1],p[S2])
	curr_lab = c(rep(1,length(S1)),rep(0,length(S2)))
	inds = sample(1:length(curr_scores))
	curr_scores = curr_scores[inds];curr_lab=curr_lab[inds]
	roc = calcAupr(curr_scores,curr_lab,roc=T)
	aupr = calcAupr(curr_scores,curr_lab,roc=F)
	if(length(S1)<minSize||length(S2)<minSize){
		roc=NA;aupr=NA
	}
	return(c(roc=roc,aupr=aupr,num_S1 = length(S1),num_S2 = length(S2)))
}

get_network_edge_based_evaluation<-function(Y,N,P){
	scores = c()
	for(i in 1:nrow(N)){
		son = N[i,2];parent=N[i,1]
		if(!is.element(son,set=colnames(Y)) || !is.element(parent,set=colnames(Y))){
			scores = rbind(scores,c(N[i,],NA,NA))
			next
		}
		curr_scores = evaluate_an_edge_vs_scores(parent,son,Y,N,P[,son])
		scores = rbind(scores,c(N[i,],curr_scores))
	}
	return(scores)
}
get_network_edge_based_evaluation_for_gene<-function(G,N,Y){
	scores = c()
	for(i in 1:nrow(N)){
		son = N[i,2];parent=N[i,1]
		if(!is.element(son,set=colnames(Y)) || !is.element(parent,set=colnames(Y))){
			scores = rbind(scores,c(N[i,],NA,NA))
			next
		}
		curr_scores = evaluate_an_edge_vs_scores(parent,son,Y,N,G)
		scores = rbind(scores,c(N[i,],curr_scores))
	}
	return(as.numeric(scores[,3]))
}
get_network_from_Y<-function(Y){
	n = ncol(Y)
	YY =  t(Y)%*%Y
	nY =  colSums(Y)
	m = matrix(0,nrow=n,ncol=n)
	colnames(m) = colnames(Y)
	rownames(m) = colnames(Y)
	# all parent - son relations
	inds = c()
	for(i in 1:n){
		for(j in 1:n){
			if(i==j || YY[i,j]==0){next}
			if (YY[i,j] < nY[j]){next}
			# now all samples in j are also in i
			inds = rbind(inds,c(i,j))
			# mark that i is an ancestor of j
			m[i,j]=1
		}
	}
	edges = cbind(colnames(Y)[inds[,1]],colnames(Y)[inds[,2]])

	to_rem = c()
	for(l in 1:nrow(inds)){
		i = inds[l,1];j=inds[l,2]
		# if there is another k such that j is a parent of k and i remove the i,j node
		for(k in 1:n){
			if (m[k,j]>0 && m[i,k]>0){
				to_rem = c(to_rem,l)
				break
			}
		}
	}
	inds = inds[-to_rem,]
	edges = cbind(colnames(Y)[inds[,1]],colnames(Y)[inds[,2]])
	return(edges)
}

### Other auxiliary methods
# Stouffer's method from 
# http://stats.stackexchange.com/questions/127328/stouffers-method-in-log-space-a-problem-with-my-r-implementation
# here modified - we apply flooring og the p-values
Stouffer.test <- function(p, w) { # p is a vector of p-values
	# p: is a vector of p-values
    # w: Vector of weights
	p = pmax(p,1e-5)
	p = pmin(p,0.99999)
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 2*(1-pnorm(abs(Z)))
  return(c(Z = Z, p.value = p.val))
}
weighted_avg<-function(x,w){
	return(sum(x*w)/sum(w))
}
# A report for a single gene
# v,y,d must be defined on the same sample set, same names
get_profile_by_datasets<-function(v,y,d,col1="blue",col2="red",minSize=10){
	curr_datasets = unique(d[y==1])
	in_datasets = sapply(d,is.element,set=curr_datasets)
	summary_groups = list()
	summary_groups[[paste("P:",sum(y==1),sep="")]] = v[y==1]
	summary_groups[[paste("N:",sum(y==0 & in_datasets),sep="")]] = v[y==0 & in_datasets]
	summary_groups[[paste("BGC:",sum(y==0 & !in_datasets),sep="")]] = v[y==0 & !in_datasets]
	y = y[in_datasets]
	v = v[in_datasets]
	d = d[in_datasets]
	l = list();cols=c()
	ps = c();ps_one_sided = c()
	for(ds in curr_datasets){
		curr_pos = y==1 & d==ds
		curr_neg = y==0 & d==ds
		if(sum(curr_pos)<minSize && sum(curr_neg)<minSize){next}
		if(sum(curr_pos)>=minSize){
			l[[paste(ds,"P",sum(curr_pos))]] = v[curr_pos]
			cols = c(cols,col1)
		}
		if(sum(curr_neg)>=minSize){
			l[[paste(ds,"N",sum(curr_neg))]] = v[curr_neg]
			cols = c(cols,col2)
		}
		if(sum(curr_pos)>=minSize && sum(curr_neg)>=minSize){
			ps[ds] = wilcox.test(v[curr_pos],v[curr_neg])$p.value
			ps_one_sided[ds] = wilcox.test(v[curr_pos],v[curr_neg],alternative="greater")$p.value
		}

	}
	return(list(values=l,cols=cols,ps=ps,summary_groups=summary_groups,ps_one_sided=ps_one_sided))
}
get_scores_vector_based_report<-function(x,Y,currterm,D,toplot=T){
	res = get_profile_by_datasets(x,Y[,currterm],D)
	# ROC scores
	positives = sample(rownames(Y)[Y[,currterm]==1])
	datasets = unique(D[positives])
	in_datasets = is.element(D[rownames(Y)],set=datasets)
	negatives = sample(rownames(Y)[Y[,currterm]==0 & in_datasets])
	bgcs = sample(rownames(Y)[Y[,currterm]==0 & !in_datasets])
	samps = sample(c(positives,bgcs))
	roc1<-calcAupr(x[samps],Y[samps,currterm],roc = TRUE)
	study_pvals = res[["ps"]]
	study_smp = Stouffer.test(res[["ps_one_sided"]])["p.value"]
	rep_score = sum(study_pvals<0.05)/length(study_pvals)
	if(toplot){
		par(mfrow=c(2,2),mar=c(5,5,5,5))
		boxplot(res[[4]],main=paste("Positive-BGC ROC:", format(roc1,digits=3)))
		boxplot(res[[1]],col=res[[2]],las=2,
			main=paste(
				"Meta-analysis p:", format(study_smp,digits=4),"\n",
				"REP score:", format(rep_score,digits=3)),
			horizontal=T)
	}
	# use prroc to get auc values
	fg = x[positives];bg = x[negatives]
	roc<-roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
	pr<-pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
	pr2 = calcAupr(c(fg,bg),c(rep(1,length(fg)),rep(0,length(bg))))
	if(toplot){
		# use rocr to plot (added because prroc initializes the figure)
		samples = sample(c(positives,negatives))
		pred <- prediction( x[samples], Y[samples,currterm])
		perf <- performance(pred,"tpr","fpr")
		plot(perf,main=paste("PN-ROC curve, AUC:",format(roc[2],digits=3)))
		perf1 <- performance(pred, "prec", "rec")
		plot(perf1,main=paste("PN-AUPR curve, AUC:",format(pr[2],digits=3)),ylim=c(0,1),xlim=c(0,1))
		QAval = abs(calcAupr(x[samples], Y[samples,currterm]) - pr[2]$auc) < 0.01
		#plot(pr)
	}
	res[["pnroc"]] = roc;res[["pnaupr"]] = pr;res[["pnaupr2"]] = pr2
	return(res)
}
get_gene_study_based_report<-function(x,curr_gene,Y,currterm,D,toplot=T){
	res = get_profile_by_datasets(x[,curr_gene],Y[,currterm],D)
	# ROC scores
	positives = sample(rownames(Y)[Y[,currterm]==1])
	datasets = unique(D[positives])
	in_datasets = is.element(D[rownames(Y)],set=datasets)
	negatives = sample(rownames(Y)[Y[,currterm]==0 & in_datasets])
	bgcs = sample(rownames(Y)[Y[,currterm]==0 & !in_datasets])
	samps = sample(c(positives,bgcs))
	roc1<-calcAupr(x[samps,curr_gene],Y[samps,currterm],roc = TRUE)
	study_pvals = res[["ps"]]
	study_smp = Stouffer.test(res[["ps_one_sided"]])["p.value"]
	rep_score = sum(study_pvals<0.05)/length(study_pvals)
	if(toplot){
		par(mfrow=c(2,2),mar=c(5,5,5,5))
		boxplot(res[[4]],main=paste("Positive-BGC ROC:", format(roc1,digits=3)))
		boxplot(res[[1]],col=res[[2]],las=2,
			main=paste(
				"Meta-analysis p:", format(study_smp,digits=4),"\n",
				"REP score:", format(rep_score,digits=3)),
			horizontal=T)
	}
	# use prroc to get auc values
	fg = x[positives,curr_gene]
	bg = x[negatives,curr_gene]
	roc<-roc.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
	pr<-pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
	if(toplot){
		# use rocr to plot (added because prroc initializes the figure)
		samples = sample(c(positives,negatives))
		pred <- prediction( x[samples,curr_gene], Y[samples,currterm])
		perf <- performance(pred,"tpr","fpr")
		plot(perf,main=paste("PN-ROC curve, AUC:",format(roc[2],digits=3)))
		perf1 <- performance(pred, "prec", "rec")
		plot(perf1,main=paste("PN-AUPR curve, AUC:",format(pr[2],digits=3)))
		QAval = abs(calcAupr(x[samples,curr_gene], Y[samples,currterm]) - pr[2]$auc) < 0.01
	}
	res[["pnroc"]] = roc
	res[["pnaupr"]] = pr
	return(res)
}
############ Format term names #############
format_label_names<-function(x){
	is_star = grepl(x,pattern = "_star$")
	x = gsub(x,pattern = "_star$",replace="")
	for(i in 1:length(x)){
		do1 = x[i]
		if(!is.null(DOTERM[[do1]])){
			currname = Term(DOTERM[[do1]])
			do1 = paste(do1,currname,sep=';')
		}
		x[i] = do1
	}
	is_doid = grepl(x,pattern = "^DOID")
	for(i in which(is_doid)){
		do1 = strsplit(x[i],split=";")[[1]][-1]
		do1 = paste(do1,collapse=';')
		x[i] = do1
	}
	has_tissue = grepl(x,pattern=';')
	x[is_star & !has_tissue] = paste(x[is_star & !has_tissue],'*',sep="")
	x[is_star & has_tissue] = gsub(x[is_star & has_tissue],pattern=';',replace='*;')
	return(x)
}
#format_label_names(colnames(y))[201]
#colnames(y)[201]
#############################################
get_wilcoxon_pval<-function(x,g,group1_val=NULL,...){
	if(is.null(group1_val)){group1_val = g[1]}
	x1 = x[g==group1_val]
	x2 = x[g!=group1_val]
	return(wilcox.test(x1,x2,...)$p.value)
}
# This function returns a matrix for each dataset: genes vs. p-values
get_study_based_analysis_for_a_column<-function(j,Y,x,D,minSize=10,...){
	y = Y[,j]
	pos_samples = names(which(y==1))
	curr_datasets = unique(D[pos_samples])
	res = list()
	for(ds in curr_datasets){
		curr_samples = names(which(D==ds))
		curry = as.character(y[curr_samples])
		if(length(unique(curry)) <2  || min(table(curry))<10){next}
		currx = x[curr_samples,]
		two_tail = apply(currx,2,get_wilcoxon_pval,g=curry,group1_val="1",alternative="two.sided")
		one_tail = apply(currx,2,get_wilcoxon_pval,g=curry,group1_val="1",alternative="greater")
		res[[ds]] = cbind(two_tail,one_tail)
	}
	return(res)
}
##########################################################################
load(CLASSIFICATION_RESULTS_FILE)
load(DATA_FILE)
load('sample2study.RData')
library('DO.db');library(bnlearn)
names(classification_results)
label_names = colnames(y)

COMPARE_CLASSIFIERS = T
classifier2selected_diseases = list()
classifier_pnrocs = c();classifier_pbrocs = c()
classifier_smqs = c();classifier_reps = c()
classifier_edgerocs = c()
if(COMPARE_CLASSIFIERS){
	method_names = c()
	for(ind in 1:length(classification_results)){
		currname = names(classification_results)[ind]
		results = classification_results[[ind]]
		if(ind==1){label_names = intersect(label_names,colnames(results$predictions))}
		if(length(intersect(label_names,colnames(results$predictions)))<length(label_names)){next}
		P = results$predictions[,label_names]
		Y = y[rownames(results$y),label_names]
		D = sample2study[rownames(Y)]
		Dtissue = sample2tissue_slim[rownames(Y)]
		CV = results$cv_batch_info
		DOnet = getBayesianNetwork(colnames(Y))
		mapping = DOnet[[2]]
		asgraph = get_graph_for_cytoscape(arcs(DOnet[[1]]),mapping)

		# Add star nodes?
		if (!any(grepl(colnames(Y),pattern = "_star$",perl=T))){
			obj = add_star_dummy_nodes(asgraph,Y)
			star_nodes = colnames(obj$m)
			star_nodes = gsub(star_nodes,pattern = "_star",replace="")
			Y = cbind(Y,obj$m)
			newP = P[,star_nodes]
			colnames(newP) = colnames(obj$m)
			P = cbind(P,newP)
			asgraph = obj$edges
		}

		# 1. Study-based
		stud_based_analysis = get_matrices_study_based_evaluation(P,Y,D)
		smqs1 = p.adjust(sapply(stud_based_analysis,function(x)unname(Stouffer.test(x[,"ps"])[2])))
		smqs2 = p.adjust(sapply(stud_based_analysis,function(x)unname(Stouffer.test(x[,"ps"],x[,"ws"])[2])))
		node2ps = sapply(stud_based_analysis,function(x)x[,"ps"])
		node2num_pos = sapply(stud_based_analysis,function(x)sum(x[,"num_pos"]))
		label2num_pos = colSums(Y)
		study_based_analysis_prop = node2num_pos/ label2num_pos[names(node2num_pos)]
		study_based_reps = sapply(node2ps,function(x)sum(x<=0.05)/length(x))
		study_based_terms = names(which((study_based_reps>=0.5 & smqs1<=0.1) | study_based_analysis_prop < 0.2))
 		print(all(names(study_based_analysis_prop)==names(study_based_reps))) # QA
		plot(-log(smqs2),study_based_reps)
		classifier_smqs = cbind(classifier_smqs,smqs2)
		classifier_reps = cbind(classifier_reps,study_based_reps)
		
		# 2. Edge-based
		net_edge_based_analysis = get_network_edge_based_evaluation(Y,asgraph,P)
		sons = unique(net_edge_based_analysis[,2])
		son2rocs = sapply(sons,function(x,y){v = as.numeric(y[y[,2]==x,3]);names(v)=y[y[,2]==x,1];return(v)},y=net_edge_based_analysis)
		son2rocs = lapply(son2rocs,function(x)x[!is.na(x)])
		edge_based_terms = names(which(unlist(lapply(son2rocs,median)>0.7) | unlist(lapply(son2rocs,length)==0)))
		edge_control_terms = edge_based_terms[grepl("control",edge_based_terms)]
		edge_based_terms_concise = union(names(which(unlist(lapply(son2rocs,median)>0.7))),edge_control_terms)
		classifier_edgerocs[[currname]] = son2rocs
		table(sapply(son2rocs,length))
		
		# 3. Node-based
		auc_scores = getMultitaskPerformancePositivesVsNegatives(P,Y,D)
		auc_scores = sapply(auc_scores,function(x)x)
		node_based_terms = names(which(auc_scores[,1]>0.7 & auc_scores[,3]>0.7))
		classifier_pnrocs = cbind(classifier_pnrocs,auc_scores[,1])
		classifier_pbrocs = cbind(classifier_pbrocs,auc_scores[,3])
		table(auc_scores[,1]>0.7)
		
		# Summary
		library(gplots);library('VennDiagram')
		venn(list(study_based=study_based_terms,edge_based=edge_based_terms,node_based=node_based_terms))
		#venn.diagram(list(study_based=study_based_terms,edge_based=edge_based_terms,node_based=node_based_terms))
		int1 = intersect(study_based_terms,edge_based_terms_concise)
		int2 = intersect(int1,node_based_terms)
		classifier2selected_diseases[[currname]] = int2
		method_names = c(method_names,currname)

		# QA and specific examples
		setdiff(edge_based_terms,edge_based_terms_concise)
		###
		smqs2["DOID:162;cancer;immune system or blood"]
		son2rocs[["DOID:162;cancer;immune system or blood"]]
		auc_scores["DOID:162;cancer;immune system or blood",]
		study_based_reps[["DOID:162;cancer;immune system or blood"]]
		###
		smqs2["DOID:1289;neurodegenerative disease;immune system or blood"]
		son2rocs[["DOID:1289;neurodegenerative disease;immune system or blood"]]
		auc_scores["DOID:1289;neurodegenerative disease;immune system or blood",]
		study_based_reps[["DOID:1289;neurodegenerative disease;immune system or blood"]]
		###
		smqs2["DOID:14566;disease of cellular proliferation;skin"]
		son2rocs[["DOID:14566;disease of cellular proliferation;skin"]]
		auc_scores["DOID:14566;disease of cellular proliferation;skin",]
		study_based_reps[["DOID:14566;disease of cellular proliferation;skin"]]
		smqs2["DOID:14566;disease of cellular proliferation;skin"]
		son2rocs[["DOID:1909;melanoma;skin"]]
		son2rocs[["DOID:1909"]]
		auc_scores["DOID:1909;melanoma;skin" ,]
		study_based_reps[["DOID:1909;melanoma;skin" ]]
		node2ps[["DOID:1909;melanoma;skin"]]
		all_mel_samples = rownames(Y)[Y[,"DOID:1909"]==1]
		skin_mel_samples = rownames(Y)[Y[,"DOID:1909;melanoma;skin"]==1]
		nonskin_mel_samples = setdiff(all_mel_samples,skin_mel_samples)
		sample2study[nonskin_mel_samples]
		table(sample2study[all_mel_samples])
		skcm_samples = names(which(sample2study=="SKCM"))
		table(Y[skcm_samples,"DOID:1909"])
		uvm_samples = names(which(sample2study=="UVM"))
		table(Y[uvm_samples,"DOID:1909"])
		asgraph[grepl("DOID:1909",asgraph[,1])|grepl("DOID:1909",asgraph[,2]),]
		test_graph = net_edge_based_analysis[grepl("DOID:1909",asgraph[,1])|grepl("DOID:1909",asgraph[,2]),][4,]
		test_parent = test_graph[1]
		test_son = test_graph[2]
		evaluate_an_edge_vs_scores(test_parent,test_son,Y,asgraph,P[,test_son])
		break
	}
}

# QA: look at disease specific reports
get_scores_vector_based_report(P[,"DOID:1909"],Y,"DOID:1909",sample2study[rownames(Y)],toplot=T)
get_scores_vector_based_report(P[,"DOID:1909;melanoma;skin"],Y,"DOID:1909;melanoma;skin",sample2study[rownames(Y)],toplot=T)
get_scores_vector_based_report(P[,"DOID:1289;neurodegenerative disease;immune system or blood"],Y,
	"DOID:1289;neurodegenerative disease;immune system or blood",sample2study[rownames(Y)],toplot=T)

apply(classifier_pnrocs,2,median,na.rm=T)
apply(classifier_pbrocs,2,median,na.rm=T)
sapply(classifier2selected_diseases,length)
apply(classifier_reps,2,sum,na.rm=T)
classifier_scores_matrices = list()
classifier_scores_matrices[["classifier_pnrocs"]] = classifier_pnrocs
classifier_scores_matrices[["classifier_pbrocs"]] = classifier_pbrocs
classifier_scores_matrices[["classifier2selected_diseases"]] = classifier2selected_diseases
classifier_scores_matrices[["classifier_edgerocs"]] = classifier_edgerocs
classifier_scores_matrices[["classifier_reps"]] = classifier_reps
classifier_scores_matrices[["classifier_smqs"]] = classifier_smqs
save(classifier_scores_matrices,file="classification_performance_scores_summary.RData")

# load the data and analyze
load("classification_performance_scores_summary.RData")
ind=1
pnrocs = classifier_scores_matrices[["classifier_pnrocs"]][,ind]
pbrocs = classifier_scores_matrices[["classifier_pbrocs"]][,ind]
edgerocs = classifier_scores_matrices[["classifier_edgerocs"]][[ind]]
edgerocs = unlist(lapply(edgerocs,median,na.rm=T))
edgerocs = edgerocs[!is.na(edgerocs)]
repscores = classifier_scores_matrices[["classifier_reps"]][,ind]
plot(edgerocs,pnrocs[names(edgerocs)]);abline(0,1)
wilcox.test(edgerocs,pnrocs[names(edgerocs)],paired=T)
plot(edgerocs,pbrocs[names(edgerocs)]);abline(0,1)
wilcox.test(edgerocs,pbrocs[names(edgerocs)],paired=T)
wilcox.test(edgerocs,repscores[names(edgerocs)],paired=T)
plot(edgerocs,repscores[names(edgerocs)])


#auc_scores_tissues = getMultitaskPerformancePositivesVsNegatives(P,Y,Dtissue)
#auc_scores_tissues = sapply(auc_scores_tissues,function(x)x)
#global_aupr = calcAupr(c(P),c(Y),roc=F)
#global_roc = calcAupr(c(P),c(Y),roc=T)
# for each term get the number of positives, negatives and bgcs
#term_stats = c()
#for(term in rownames(auc_scores)){
#	positives = rownames(Y)[Y[,term]==1]
#	datasets = unique(D[positives])
#	negatives = names(D[is.element(D,set=datasets)])
#	negatives = setdiff(negatives,positives)
#	bgcs = setdiff(rownames(Y),union(positives,negatives))
#	term_stats = rbind(term_stats,c(length(positives),length(negatives),length(datasets),length(bgcs)))	
#}
#colnames(term_stats) = c("#positives","#negatives","#datasets","#bgcs")
#auc_scores = cbind(auc_scores,term_stats)
#m = auc_scores[!apply(is.nan(auc_scores),1,any),]
#term_stats_tissues = c()
#for(term in rownames(auc_scores_tissues)){
#	positives = rownames(Y)[Y[,term]==1]
#	datasets = unique(Dtissue[positives])
#	negatives = names(Dtissue[is.element(Dtissue,set=datasets)])
#	negatives = setdiff(negatives,positives)
#	bgcs = setdiff(rownames(Y),union(positives,negatives))
#	term_stats_tissues = rbind(term_stats_tissues,c(length(positives),length(negatives),length(datasets),length(bgcs)))	
#}
#colnames(term_stats_tissues) = c("#positives","#negatives","#datasets","#bgcs")
#rownames(term_stats_tissues) = rownames(auc_scores_tissues)
#auc_scores_tissues = cbind(auc_scores_tissues,term_stats_tissues)
#m = auc_scores_tissues[!apply(is.nan(auc_scores_tissues),1,any),]

#######################################################################################################################################################################
# Gene-based analysis
Y = y
D = sample2study[rownames(Y)]
x = x[rownames(Y),]
gc()

# If you want to exclude the gtex samples then run this code:
all_samples = intersect(rownames(x),rownames(y))
gtex_samples = grepl(all_samples,pattern='gtex',ignore.case=T)
all_samples = all_samples[!gtex_samples]
x = x[all_samples,];y=y[all_samples,];Y=Y[all_samples,]
sample2study = sample2study[all_samples]
D = D[all_samples]
gc()

# Step1: Node-based analysis
pn_roc_scores = c()
for(j in 1:ncol(Y)){
	print(j)
	curry = Y[,j]
	curr_datasets = unique(D[curry==1])
	in_datasets_inds = sample(which(sapply(D,is.element,set=curr_datasets)))
	curry = curry[in_datasets_inds];currx=x[in_datasets_inds,]
	if(length(unique(curry))<2){next}
	curr_scores = apply(currx,2,calcAupr,gs=curry,roc=T)
	pn_roc_scores = cbind(pn_roc_scores,curr_scores)
	colnames(pn_roc_scores)[ncol(pn_roc_scores)] = colnames(Y)[j]
	if(any(is.na(c(pn_roc_scores)))){break}
	gc()
}
rownames(pn_roc_scores) = colnames(x)
save(pn_roc_scores,file = 'gene_pn_roc_scores_wo_gtex.RData')
pb_roc_scores = c()
for(j in 1:ncol(Y)){
	print(j)
	curry = Y[,j]
	curr_datasets = unique(D[curry==1])
	in_datasets = sapply(D,is.element,set=curr_datasets)
	curr_inds = (curry==1 | !in_datasets)
	print(table(curr_inds))
	curr_inds = sample(which(curr_inds))
	curry = curry[curr_inds];currx=x[curr_inds,]
	if(length(unique(curry))<2){next}
	curr_scores = apply(currx,2,calcAupr,gs=curry,roc=T)
	pb_roc_scores = cbind(pb_roc_scores,curr_scores)
	colnames(pb_roc_scores)[ncol(pb_roc_scores)] = colnames(Y)[j]
	gc()
}
rownames(pb_roc_scores) = colnames(x)
save(pb_roc_scores,file = 'gene_pb_roc_scores_wo_gtex.RData')

# Compare the ROCs with and without the gtex samples
yyy = get(load("classification_performance_scores_summary.RData"))
currterms = yyy$classifier2selected_diseases[[1]]
x1 = get(load('gene_pb_roc_scores.RData'))
x2 = get(load('gene_pb_roc_scores_wo_gtex.RData'))
col_inter = intersect(colnames(x1),colnames(x2))
x1 = x1[,col_inter];x2=x2[,col_inter];gc()
corrs = cor(x1,x2)
hist(diag(corrs))
hist(diag(corrs[currterms,currterms]))
rand_ind = sample(1:ncol(x1))[1]
plot(y=x1[,rand_ind],x2[,rand_ind],main=colnames(x1)[rand_ind],ylab="With GTeX",xlab="Without GTex");abline(0,1)


# Step 2: we use the results of step 1 to filter genes
# For each disease we create matrices of genes vs. their scores
num_threads = 40
inds = 1:ncol(Y)
gene_dataset_p_matrices = list()
for(i in inds){gene_dataset_p_matrices[[colnames(Y)[i]]] = NULL}
while(length(inds)>0 && !is.null(inds) && !is.na(inds)){
	curres = mclapply(inds,get_study_based_analysis_for_a_column,mc.cores=num_threads,Y=Y,x=x,D=D)
	names(curres) = colnames(Y)[inds]
	for (nn in names(curres)){
		gene_dataset_p_matrices[[nn]] = curres[[nn]]
	}
	inds = which(grepl(sapply(gene_dataset_p_matrices,class),pattern='error'))
	save(gene_dataset_p_matrices,file='gene_dataset_p_matrices.RData')
	curres = NULL
	gc()
}

# Edge-based analysis
DOnet = getBayesianNetwork(colnames(Y))
mapping = DOnet[[2]]
asgraph = get_graph_for_cytoscape(arcs(DOnet[[1]]),mapping)
net_edge_based_analysis = apply(x,2,get_network_edge_based_evaluation_for_gene,N=asgraph,Y=Y)
save(net_edge_based_analysis,file="gene_net_edge_based_analysis.RData")
sons = unique(asgraph[,2])
# network: edges are pairs of parent,son
# y: a matrix of scores, each row corresponds to an edge
get_son_edge_scores<-function(son,y,network){
	m = t(as.matrix(y[network[,2]==son,]))
	if(nrow(m)==1){m=t(m)}
	print(dim(m))
	colnames(m) = network[network[,2]==son,1]
	return(m)
}
son2rocs = list()
for(son in sons){
	m = get_son_edge_scores(son,y=net_edge_based_analysis,network=asgraph)
	to_rem = apply(is.na(m),2,all)
	m = m[,!to_rem]
	son2rocs[[son]] = as.matrix(m)
}
son2rocs = son2rocs[sapply(son2rocs,length)>0]
save(son2rocs,file='gene_edge_based_son2rocs.RData')
sapply(son2rocs,dim)
table(sapply(son2rocs,length))

##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################

# Summarize
# (1) For PN and PB ROCS use a threshold on the two-tailed statistic
# (2) For the study-based analysis we use a simple REP score and a sign test
# (3) For the edge-based analysis we use a threshold on the ROC that should apply to all siblings
load('gene_pb_roc_scores.RData')
load('gene_pn_roc_scores.RData')
load('gene_dataset_p_matrices.RData')
load('gene_study_roc_matrices.RData')
load('gene_edge_based_son2rocs.RData')
load("classification_performance_scores_summary.RData")
classifier2selected_diseases = classifier_scores_matrices[["classifier2selected_diseases"]]
# Edge-based analysis
DOnet = getBayesianNetwork(colnames(Y))
mapping = DOnet[[2]]
asgraph = get_graph_for_cytoscape(arcs(DOnet[[1]]),mapping)


class2two_sided_pvals = c()
for(i in 1:length(gene_dataset_p_matrices)){
	print(i)
	l = gene_dataset_p_matrices[[i]]
	if(length(l)==0){next}
	m = c()
	for(j in 1:length(l)){
		m = cbind(m,l[[j]][,1])
	}
	colnames(m) = names(l)
	print(dim(m))
	class2two_sided_pvals[[names(gene_dataset_p_matrices)[i]]] = m
}

# QA for observed NAs - should be cases without positives
#terms_with_nas = names(which(sapply(class2two_sided_pvals,function(x)sum(is.na(x)))>0))
#currterm = sample(terms_with_nas)[1]
#mat = class2two_sided_pvals[[currterm]]
#table(apply(is.na(mat),1,sum))
#table(apply(is.na(mat),2,sum))
#currgenes = names(which(apply(is.na(mat),1,sum)>0))
#currgene = sample(currgenes)[1]
#mat[currgene,]
#get_gene_study_based_report(x,currgene,Y,currterm,D,toplot=T)
#curr_profile = get_profile_by_datasets(x[,currgene],Y[,currterm],D)
#dataset_values = curr_profile$values
#sapply(dataset_values,length)
#curr_profile$ps
#v1 = dataset_values[["KICH P 66"]]
#v2 = dataset_values[["KICH N 25"]]
#table(x[,currgene]==0)
#plot(sort(x[,currgene]))
# QA for NAs in the PN-ROC
mat = pn_roc_scores
table(is.na(mat))
table(apply(is.na(mat),1,sum))
which(apply(is.na(mat),2,sum)>0)

# Correct NAs and NaNs in the matrix
# Put 1 p-values
for(j in 1:length(class2two_sided_pvals)){
	class2two_sided_pvals[[j]][is.na(class2two_sided_pvals[[j]])] = 1
	class2two_sided_pvals[[j]][is.nan(class2two_sided_pvals[[j]])] = 1
}

class2one_sided_pvals = c()
for(i in 1:length(gene_dataset_p_matrices)){
	print(i)
	l = gene_dataset_p_matrices[[i]]
	if(length(l)==0){next}
	m = c()
	for(j in 1:length(l)){
		m = cbind(m,l[[j]][,2])
	}
	colnames(m) = names(l)
	print(dim(m))
	class2one_sided_pvals[[names(gene_dataset_p_matrices)[i]]] = m
}
############################################################################################

p_thr = 0.05
rep_thr = 0.5
roc_thr = 0.7
q_thr = 0.01
gene_sets = list()
disease2summary_scores = list()
exclude_gtex=T
currterms = classifier2selected_diseases[[1]]
#currterms = colnames(Y)
for(currterm in currterms){
	currm = c()
	print(currterm)
	if(grepl(currterm,"_star")){break}
	# Validation type 1: PN, PB
	pn_roc = pn_roc_scores[,currterm]
	pn_roc = pmax(pn_roc,1-pn_roc)
	pb_roc = pb_roc_scores[,currterm]
	pb_roc = pmax(pb_roc,1-pb_roc)
	roc_selected_genes = names(pn_roc)[pn_roc>=roc_thr & pb_roc>=roc_thr]
	currm = cbind(pn_roc,pb_roc)

	# Validation type 2: study-based
	curr_set = roc_selected_genes
	pvals = class2two_sided_pvals[[currterm]]
	pvals1 = class2one_sided_pvals[[currterm]]
	if(length(pvals1) >0 && is.null(dim(pvals1))){pvals1 = as.matrix(pvals1,ncol=1)}
	if(length(pvals) >0 && is.null(dim(pvals))){pvals = as.matrix(pvals,ncol=1)}
	if(exclude_gtex){
		pvals = pvals[,!grepl(colnames(pvals),pattern="GTEx")]
		pvals1 = pvals1[,!grepl(colnames(pvals1),pattern="GTEx")]
	}
	if(length(pvals1) >0 && is.null(dim(pvals1))){pvals1 = as.matrix(pvals1,ncol=1)}
	if(length(pvals) >0 && is.null(dim(pvals))){pvals = as.matrix(pvals,ncol=1)}
	if(length(pvals) >0 && !is.null(dim(pvals))){
		smqs = apply(pvals1,1,Stouffer.test)["p.value",][rownames(currm)]
		smqs = p.adjust(smqs,method="BY")
		pval_scores = (apply(pvals<=p_thr,1,sum)/ncol(pvals))[rownames(currm)]
		pval_selected_genes = names(pval_scores)[pval_scores>=rep_thr &smqs<=0.01]
		curr_set = intersect(roc_selected_genes,pval_selected_genes)
		currm = cbind(currm,cbind(pval_scores,smqs))
	}
	getSYMBOL(curr_set,data='org.Hs.eg')
	# Validation type 3: edge-based
	edge_based_scores = son2rocs[[currterm]]
	edge_based_scores = pmax(edge_based_scores,1-edge_based_scores)
	if(!is.null(edge_based_scores) && length(edge_based_scores)>0){
		child_specific_genes = apply(edge_based_scores>=roc_thr,1,all)
		child_specific_genes = names(child_specific_genes)[child_specific_genes]
		curr_set = intersect(curr_set,child_specific_genes)
		edge_based_scores = as.matrix(edge_based_scores[rownames(currm),])
		currm = cbind(currm,apply(edge_based_scores,1,min))
	}
	gene_sets[[currterm]] = curr_set 
	disease2summary_scores[[currterm]] = currm
}
save(gene_sets,disease2summary_scores,file="selected_genes_adeptus2.RData")
gene_sets[["DOID:1909"]]

# correlation among the scores
score_corrs = list()
for (nn in names(disease2summary_scores)){
	scores = disease2summary_scores[[nn]]
	na_inds = apply(is.na(scores),1,any)
	if(all(na_inds)){next}
	score_corrs[[nn]] =  cor(scores[!na_inds,],method='spearman')
}
score_corrs[["DOID:1909"]]
save(gene_sets,score_corrs,disease2summary_scores,file="selected_genes_adeptus2.RData")

# QA and reports for some specific genes
selected_diseases = classifier_scores_matrices[["classifier2selected_diseases"]][[1]]
currterm = sample(selected_diseases)[1]
DOTERM[[currterm]]
x1 = disease2summary_scores[[currterm]][,"smqs"]
x2 = pn_roc_scores[names(x1),currterm]
x3 = disease2summary_scores[[currterm]][,"pval_scores"]
x4 = disease2summary_scores[[currterm]][,5]
chisq.test(table(x1<0.01,pmax(x2,1-x2)>0.7))$p.value
chisq.test(table(x3>0.5,pmax(x2,1-x2)>0.7))$p.value
chisq.test(table(x3,x1<0.01))$p.value
mat = cbind(x1,x3,x4)
library(scatterplot3d)
#x2 = pmax(x2,1-x2)
cols = rep("black",length(x1))
cols[x3>=0.5 & pmax(x2,1-x2)>0.7 & x4>0.7] = "blue"
main = currterm
if (!is.null(DOTERM[[currterm]])){
	main = paste(currterm,Term(DOTERM[[currterm]]))
}
scatterplot3d(x4,x2,x3,color=cols,
	xlab = "Disease specificity",ylab="Signal intensity",zlab="Replicability",
	main = main
)
pairs(cbind(x4,x2,x3),col=cols)

# case 1: bad ROC, good meta-analysis
curr_gene = names(which(abs(x2-0.5) <= 0.01 & x1 <0.001)[1])
# case 2: good neg ROC, good meta-analysis
curr_gene = names(which(0.5-x2 >= 0.35 & x1 <0.001)[1])
# case 3: good positive ROC, good meta-analysis
curr_gene = names(which(x2-0.5 >= 0.4 & x1 <0.001)[1])
# case 4: 
curr_gene = names(which(abs(x2-0.5) <= 0.01 & x1 >0.5)[1])
get_gene_study_based_report(x,curr_gene,Y,currterm,D,toplot=T)

# Look at the jaccard among the sets of the selected diseases
n = length(selected_diseases)
Js = matrix(0,n,n)
colnames(Js) = sapply(selected_diseases,format_label_names)
rownames(Js) = colnames(Js)
for(i in 1:n){
	set1 = gene_sets[[selected_diseases[i]]]
	if(length(set1)==0 || is.null(set1)){next}
	for(j in 1:n){
		set2=gene_sets[[selected_diseases[j]]]
		if(length(set2)==0 || is.null(set2)){next}
		jacc = length(intersect(set1,set2))/length(union(set1,set2))
		Js[i,j] = jacc
	}
}
quantile(c(Js[lower.tri(Js)]))
library(corrplot)
JJs = Js;rownames(JJs)=NULL;colnames(JJs)=NULL
corrplot(JJs,cl.lim=c(0,1))

###############################################################################
# Create the files for the webtool and text files
###############################################################################

# For learning the classifiers
load(DATA_FILE)
Y = y
D = sample2study[rownames(Y)]
X = x[rownames(Y),]
gc()
#################################

library(org.Hs.eg.db);library(annotate)
# entrez to symbol
library('org.Hs.eg.db')
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
entrez2symbol <- unlist(as.list(x[mapped_genes]))
# reverse
x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
symbol2entrez <- unlist(as.list(x[mapped_genes]))

# set a dir to contain the output text files:
newsite_loc = 'analysis_output_text_files/'
ind = 1
dir.create(newsite_loc)
load("classification_performance_scores_summary.RData")
selected_terms = classifier_scores_matrices[["classifier2selected_diseases"]][[ind]]

# Modify X and Y and run the classifiers, comment out if not needed
colnames(X) = entrez2symbol[colnames(X)]
Y = Y[,selected_terms]
colnames(Y) = sapply(colnames(Y),format_label_names)
NUM_FEATURES = 200
NUM_VAR_GENES = 5000
USE_TOP_VAR_GENES = T
sds = apply(X,2,sd)
gene_var_thr = sort(sds,decreasing=USE_TOP_VAR_GENES)[NUM_VAR_GENES]
if (USE_TOP_VAR_GENES){cols_for_analysis = sds >= gene_var_thr}
if (!USE_TOP_VAR_GENES){cols_for_analysis = sds <= gene_var_thr}
X = X[,cols_for_analysis];gc()
svm_classifiers = simple_binary_relevance(X,Y,
		baseClassifierArgs=list(classification_function=liblinear_binary_svm_wrapper,
			numFeatures=NUM_FEATURES,d=sample2study,reps=50,max_num_samples = 500),
		baseClassifierPredArgs=list(),baseClassifier=simple_sampling_based_learner)
classifier_genes = colnames(X)
all_genes = names(sds)
label_names = colnames(Y)
save(svm_classifiers,classifier_genes,label_names,entrez2symbol,symbol2entrez,all_genes,file='svm_classifiers.RData')
#########################################################################

load("selected_genes_adeptus2.RData")
# write the results of the selcted terms to a text file
for(dterm in selected_terms){
	term = format_label_names(dterm)
	genes = gene_sets[[dterm]]
	if(length(genes) == 0 ){next}
	gene_symbs = getSYMBOL(genes,data='org.Hs.eg')
	gene_symbs2 = entrez2symbol[genes]
	print(table(gene_symbs == gene_symbs2))
	if(length(genes)==0 || is.null(genes)){next}
	v = rep(term,length(genes))
	m = cbind(genes,gene_symbs,v)
	to_append = j>1
	write.table(m[,c(1,3)],file = paste(newsite_loc,"selected_genes_adeptus2.txt",sep=''),sep="\t",quote=F,row.names=F,col.names=F,append=to_append)
	write.table(m[,c(2,3)],file = paste(newsite_loc,"selected_genes_adeptus2_symbols.txt",sep=''),sep="\t",quote=F,row.names=F,col.names=F,append=to_append)
}
########## Print the classification results output for cytoscape ########
selected_terms_network = asgraph[is.element(asgraph[,1],set=selected_terms) & is.element(asgraph[,2],set=selected_terms),]
selected_terms_formatted = format_label_names(selected_terms)
selected_terms_network_formatted = cbind(format_label_names(selected_terms_network[,1]),format_label_names(selected_terms_network[,2]))
node_scores = c()
for(score_type in names(classifier_scores_matrices)){
	if(score_type == "classifier2selected_diseases"){next}
	currterms = rownames(classifier_scores_matrices[[score_type]])
	if(is.null(currterms)){
		currterms = names(classifier_scores_matrices[[score_type]][[ind]])
	}
	v = rep(NA,length(selected_terms))
	names(v) = selected_terms
	if(score_type=="classifier_edgerocs"){
		v[intersect(currterms,selected_terms)] = sapply(classifier_scores_matrices[[score_type]][[ind]][intersect(currterms,selected_terms)],min)
		v[is.infinite(v)] = NA
	}
	else{
		v[intersect(currterms,selected_terms)] = classifier_scores_matrices[[score_type]][intersect(currterms,selected_terms),ind]
	}
	node_scores = cbind(node_scores,v)
	colnames(node_scores)[ncol(node_scores)] = score_type
}
rownames(node_scores) = selected_terms_formatted
colnames(node_scores) = gsub(colnames(node_scores),pattern="classifier_",replace="")
dirname = paste(newsite_loc,"BR_classifiers_perf_output/",sep='')
dir.create(dirname)
selected_terms_network_formatted = rbind(selected_terms_network_formatted,
	cbind(selected_terms_formatted,selected_terms_formatted)
)
selected_terms_network_formatted = unique(selected_terms_network_formatted)
write.table(selected_terms_network_formatted,file=paste(newsite_loc,"selected_terms_network_formatted.txt",sep=""),
	sep="\t",quote=F,row.names=F,col.names=F)
write.table(node_scores,file=paste(newsite_loc,"selected_terms_node_scores.txt",sep=""),
	sep="\t",quote=F,row.names=T,col.names=T)
selected_terms_tissues = unique(sapply(selected_terms_network_formatted,function(x)try(strsplit(x,split=';')[[1]][2])))
selected_terms_tissues = selected_terms_tissues[!is.na(selected_terms_tissues)]
selected_terms_tissues = unique(gsub(selected_terms_tissues,pattern='\\*',replace=''))
for(tissue in selected_terms_tissues){
	curr_edges = grepl(tissue,selected_terms_network_formatted[,1]) & grepl(tissue,selected_terms_network_formatted[,2])
	curr_graph = selected_terms_network_formatted[curr_edges,]
	write.table(curr_graph,file=paste(dirname,tissue,"_selected_terms_network_formatted.txt",sep=""),
		sep="\t",quote=F,row.names=F,col.names=F)
}
all_terms = rownames(classifier_scores_matrices[[1]])
all_terms_network = asgraph
all_terms_formatted = format_label_names(all_terms)
all_terms_network_formatted = cbind(format_label_names(all_terms_network[,1]),format_label_names(all_terms_network[,2]))
colnames(all_terms_network_formatted) = c('parent','son')
write.table(all_terms_network_formatted,file=paste(newsite_loc,"all_terms_network_formatted.txt",sep=""),
	sep="\t",quote=F,row.names=F,col.names=F)
write.table(all_terms_network_formatted,file=paste(newsite_loc,"newtwork_edges.txt",sep=""),
	sep="\t",quote=F,row.names=F,col.names=F)
write.table(all_terms_formatted,file=paste(newsite_loc,"all_terms_formatted.txt",sep=""),
	sep="\t",quote=F,row.names=F,col.names=F)
write.table(selected_terms_formatted,file=paste(newsite_loc,"selected_terms_formatted.txt",sep=""),
	sep="\t",quote=F,row.names=F,col.names=F)

###########################################################
# Print the gene scores to txt files
# RELOAD THE DATA AFTER THIS CODE, IT CHANGES ALL LABEL AND GENE NAMES IN THE OBJECTS
all_labels = sapply(colnames(y),format_label_names)
# 1. ROC scores
rownames(pb_roc_scores) = entrez2symbol[rownames(pb_roc_scores)]
rownames(pn_roc_scores) = entrez2symbol[rownames(pn_roc_scores)]
colnames(pb_roc_scores) = sapply(colnames(pb_roc_scores),format_label_names)
colnames(pn_roc_scores) = sapply(colnames(pn_roc_scores),format_label_names)
# In a new version of the DB these tests may change (and must be attended if it happens)
table(is.na(rownames(pb_roc_scores))) # QA: all should be false
length(unique(rownames(pb_roc_scores))) == length(rownames(pb_roc_scores)) # QA: should be true
write.table(pb_roc_scores[,all_labels],file=paste(newsite_loc,'background_roc.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
pn_roc_scores_cols = names(which(apply(is.na(pn_roc_scores),2,sum)==0))
pn_roc_scores_cols = intersect(pn_roc_scores_cols,all_labels)
write.table(pn_roc_scores[,pn_roc_scores_cols],file=paste(newsite_loc,'negative_roc.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
# 2. Replicability scores
exclude_gtex = T
rep_scores = c()
p_thr = 0.05
all_genes = names(rownames(pb_roc_scores))
all_symbs = rownames(pb_roc_scores)
for(j in 1:ncol(Y)){
	currterm = colnames(Y)[j]
	pvals = class2two_sided_pvals[[currterm]]
	if(exclude_gtex){pvals = pvals[,!grepl(colnames(pvals),pattern="GTEx")]}
	if(length(pvals) >0 && is.null(dim(pvals))){pvals = as.matrix(pvals,ncol=1)}
	if(length(pvals) >0 && !is.null(dim(pvals))){
		pval_scores = (apply(pvals<=p_thr,1,sum)/ncol(pvals))
		rep_scores = cbind(rep_scores,pval_scores[all_genes])
		colnames(rep_scores)[ncol(rep_scores)] = currterm
	}
}
table(rep_scores < 0.5)/length(rep_scores)
colnames(rep_scores) = sapply(colnames(rep_scores),format_label_names)
rownames(rep_scores) = all_symbs
write.table(rep_scores,file=paste(newsite_loc,'replicability_scores.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
table(rep_scores < 1e-10)
# 3. specificity scores
load('gene_edge_based_son2rocs.RData')
specf_scores = c()
for(j in 1:ncol(Y)){
	currterm = colnames(Y)[j]
	edge_based_scores = son2rocs[[currterm]]
	edge_based_scores = pmax(edge_based_scores,1-edge_based_scores)
	if(!is.null(edge_based_scores) && length(edge_based_scores)>0){
		edgerocs = apply(edge_based_scores,1,min)
		specf_scores = cbind(specf_scores,edgerocs[all_genes])
		colnames(specf_scores)[ncol(specf_scores)] = currterm
	}
}
colnames(specf_scores) = sapply(colnames(specf_scores),format_label_names)
rownames(specf_scores) = all_symbs
write.table(specf_scores,file=paste(newsite_loc,'specificity_scores.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
# 4. SMQs based on one sided p-values
smqs = c()
all_genes = names(rownames(pb_roc_scores))
all_symbs = rownames(pb_roc_scores)
for(j in 1:ncol(y)){
	currterm = colnames(y)[j]
	pvals1 = class2one_sided_pvals[[currterm]][all_genes,]
	if(length(pvals1) >0 && is.null(dim(pvals1))){pvals1 = as.matrix(pvals1,ncol=1)}
	if(exclude_gtex){pvals1 = pvals1[,!grepl(colnames(pvals1),pattern="GTEx")]}
	if(length(pvals1) >0 && is.null(dim(pvals1))){pvals1 = as.matrix(pvals1,ncol=1)}
	if(length(pvals1) >0 && !is.null(dim(pvals1))){
		term_smqs = apply(pvals1,1,Stouffer.test)["p.value",]
		term_smqs = p.adjust(term_smqs,method='BY')
		smqs = cbind(smqs,term_smqs[all_genes])
		colnames(smqs)[ncol(smqs)] = currterm
	}
}
colnames(smqs) = sapply(colnames(smqs),format_label_names)
rownames(smqs) = all_symbs
write.table(smqs,file=paste(newsite_loc,'smqs_scores.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
# 5. Print the p-value matrices and info for showing gene reports
# We need to load the data to get some gene reports - needed for statistics
matrix_names = sapply(names(class2two_sided_pvals),format_label_names)
pval_data_dir = paste(newsite_loc,'pvalue_matrices/',sep='')
dir.create(pval_data_dir)
load(DATA_FILE)
load('sample2study.RData')
Y = y
D = sample2study[rownames(Y)]
x = x[rownames(Y),]
gc()
curr_gene = 1
label2metadata = c()
for (j in 1:length(matrix_names)){
	print(j)
	term = matrix_names[j]
	if(!grepl(term,pattern='\\*')){next}
	if (!is.element(names(matrix_names)[j],set=colnames(Y))){next}
	res = get_gene_study_based_report(x,curr_gene,Y,names(matrix_names)[j],D,toplot=F)
	dataset_info = t(sapply(names(res[[1]]),function(x)strsplit(x,split=" ")[[1]]))
	curr_datasets = unique(dataset_info[,1])
	dataset_info_summary = c()
	for(ds in curr_datasets){
		if (grepl(ds,pattern = "gtex",ignore.case=T)){next}
		if (sum(dataset_info[,1]==ds)<2){next}
		curr_info = dataset_info[dataset_info[,1]==ds,]
		num_samples = curr_info[,3];names(num_samples) = curr_info[,2]
		dataset_info_summary = rbind(dataset_info_summary,num_samples)
		rownames(dataset_info_summary)[nrow(dataset_info_summary)] = ds
	}
	if(length(dataset_info_summary)==0){next}
	dataset2name = c()
	for (ds in rownames(dataset_info_summary)){
		dataset2name[ds] = paste(c(ds,dataset_info_summary[ds,]),collapse=';')
	}
	group_sizes = sapply(res[[4]],length)
	mat = class2two_sided_pvals[[j]]
	if(is.null(dim(mat))){next}
	# needed because of the star terms: 
	# we need only datasets of the positives and not all of the doid
	curr_shared_datasets = intersect(colnames(mat),curr_datasets)
	mat = mat[,curr_shared_datasets]
	if(is.null(dim(mat)) && length(mat)>0){
		mat=as.matrix(mat,ncol=1)
		colnames(mat) = curr_shared_datasets
	}
	if(length(mat)==0){next}
	genes = entrez2symbol[rownames(mat)]
	mat_colnames = dataset2name[colnames(mat)]
	keep_cols = !grepl(colnames(mat),pattern='gtex',ignore.case=T)
	mat_colnames = mat_colnames[keep_cols]
	mat = mat[,keep_cols]
	if(is.null(dim(mat))){mat=as.matrix(mat,ncol=1)}
	rownames(mat) = genes
	colnames(mat) = mat_colnames
	label2metadata = rbind(label2metadata,
		c(term,length(curr_datasets),ncol(mat),group_sizes)
	)
	term = gsub(term,pattern = '/',replace='-')
	term = gsub(term,pattern = '\\*',replace='_star')
	write.table(mat,file=paste(pval_data_dir,term,'.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)
}
rownames(label2metadata) = label2metadata[,1]
label2metadata = label2metadata[,-1]
mode(label2metadata) = 'numeric'
colnames(label2metadata) = c("TotalDatasets","RepAnalysisDatasets","P","N","BGC")
write.table(label2metadata,file=paste(pval_data_dir,'label2metadata.txt',sep=''),sep="\t",quote=F,col.names=T,row.names=T)


#################################################################################


















