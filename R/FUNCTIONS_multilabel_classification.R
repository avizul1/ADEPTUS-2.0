
################## Implementation of wrappers for base classifiers ####################################
# These are used within some multilabel algorithms, most notably the naive binary relevance method.
# We implement below wrappers for binary and non binary data. These methods
# First run feature selection, record the selected clustering and use them to learn the
# classifier. 
# The non-binary data analysis is based on the CMA package in R.
# For binary data we used Fisher's exact test.
# The main reason for these wrappers: the need to save running time.
featureSelectionClassifier<-function(x,y,classification_function=randomForest,fsmethod="limma",numFeatures=100,verbose=F,...){
	selection = GeneSelection(x,y,method=fsmethod,trace=F)
	selected_inds = toplist(selection,k=numFeatures,show=F)$index
	if(verbose){print ("features selected")}
	classifier = classification_function(x[,selected_inds],y,...)
	if(verbose){print (paste("classifier learned",class(classifier)))}
	obj = list(inds=selected_inds,classifier=classifier)
	class(obj) = "featureSelectionClassifier"
	return (obj)
}

predict.featureSelectionClassifier<-function(obj,newx,...){
	#print ("in fs predict")
	args = c(list(obj$classifier,newx[,obj$inds]),...)
	return (do.call(predict,args=args))
}

get_hg_pval<-function(x,y,...){
	tab = table(x,y)
	if (length(tab) !=4){return (1)}
	return (fisher.test(tab,...)$p.value)
}

featureSelectionClassifier_binary_data<-function(x,y,classification_function=svm,numFeatures=250,verbose=F,...){
	f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'greater')
	names(f_pvals) = colnames(x)
	selected_features1 = names(sort(f_pvals)[1:numFeatures])
	f_pvals = apply(x,2,get_hg_pval,y=y,alternative = 'less')
	names(f_pvals) = colnames(x)
	selected_features2 = names(sort(f_pvals)[1:numFeatures])
	selected_features = union(selected_features1,selected_features2)
	if(verbose){print (paste(length(selected_features)," features selected"))}
	classifier = classification_function(as.matrix(x[,selected_features]),y,...)
	if(verbose){print (paste("classifier learned",class(classifier)))}
	obj = list(inds=selected_features,classifier=classifier)
	class(obj) = "featureSelectionClassifier"
	return (obj)
}
################### other tests for analysis of binary data ##########################
# These did not prove as more useful than Fisher's exact test and are not
# used by default.
get_chisq_pvalue<-function(x,y){return (chisq.test(table(x,y))$p.value)}
get_chisq_statistic<-function(x,y){return (chisq.test(table(x,y))$statistic)}
get_logit_statistic<-function(x,y){
	lfit = glm(y ~ x, family = "binomial")
	return (summary(lfit)$aic)
}
get_logit_pvalue<-function(x,y){
	lfit = glm(y ~ x, family = "binomial")
	return (summary(lfit)$coefficients[2,4])
}
get_gtest_pvalue<-function(x,y){return (g.test(table(x,y))$p.value)}
g.test = function(x,print=T) {
  row.sum = rowSums(x); col.sum = colSums(x);n = sum(x)
  # x.expected[i,j] = row.sum[i] * col.sum[j] / n
  x.expected = row.sum %o% col.sum / n
  g = 2*sum( x*log(x/x.expected) )
  degf = (nrow(x)-1)*(ncol(x)-1)
  p.g = 1 - pchisq(g,degf)
  return(invisible(list(expected=x.expected,test.stat=g,p.value=p.g)))
}
#############################################################

######### Wrapper for liblinear (FAST) ################
# The code below implements our API above for binary classification using the liblinear package.
# These are useful in cases in which it is hard to for default SVM implementations to converge.
# For SVM, we get probabilistic predictions using bootstrap sampling.

# Wrapper for SVM classification + tunning of V
liblinear_binary_svm_wrapper<-function(xTrain,yTrain,type=1){
	yTrain = as.character(yTrain)
	s=scale(xTrain,center=TRUE,scale=TRUE)
	co=heuristicC(s)
	m=LiblineaR(data=s,target=yTrain,type=type,cost=co,bias=TRUE,verbose=FALSE)
	m1=LiblineaR(data=s,target=yTrain,type=0,cost=co,bias=TRUE,verbose=FALSE)
	m1W=m$W
	m$Type=0 # enforce logistic regr prediction
	train_center = attr(s,"scaled:center")
	train_scale = attr(s,"scaled:scale")
	obj = list(m=m1,train_center=train_center,train_scale=train_scale)
	class(obj) = "liblinear_binary_svm_wrapper"
	return(obj)
}
predict.liblinear_binary_svm_wrapper<-function(obj,xTest){
	s2=scale(xTest,obj$train_center,obj$train_scale)
	p=predict(obj$m,s2,proba=T)
	return(p$probabilities)
}

# Logistic regression models + internal CV
liblinear_wrapper<-function(x,y,cross=2,scale=F,tryTypes=c(0,6),tryCosts=c(1000,1,0.001)){
	means = 0
	sds = 1
	if (scale){
		means = apply(x,2,mean)
		sds = apply(x,2,sd)
		x = apply(x,2,function(x){(x-mean(x))/sd(x)})
	}
	bestCost=NA;bestAcc=0;bestType=NA
	for(ty in tryTypes){
		for(co in tryCosts){
			acc=LiblineaR(data=x,target=y,type=ty,cost=co,bias=TRUE,cross=5,verbose=FALSE)
			#cat("Results for C=",co," : ",acc," accuracy.\n",sep="")
			if(acc>bestAcc){bestCost=co;bestAcc=acc;bestType=ty}
		}
	}
	final_model = LiblineaR(data=x,target=y,type=bestType,cost=bestCost,bias=TRUE)
	obj = list(model = final_model,type=bestType,cpst = bestCost,scale=scale,means=means,sds=sds)
	class(obj) = "liblinear_wrapper"
	return (obj)
}
predict.liblinear_wrapper<-function(obj,newx){
	if (obj$scale){
		for (j in 1:ncol(newx)){
			newx[,j] = (newx[,j] - obj$means[j])/obj$sds[j]
		}
	}
	return (predict(obj$model,newx,proba=T)$probabilities)
}

# Wrapper for fast svm regression
liblinear_svm_wrapper<-function(x,y,reps = 50, bootstrap_frac = 0.5, cross=2,scale=F,tryTypes=c(1,5),tryCosts=c(1000,1,0.001)){
	means = 0
	sds = 1
	if (scale){
		means = apply(x,2,mean)
		sds = apply(x,2,sd)
		x = apply(x,2,function(x){(x-mean(x))/sd(x)})
	}
	models = list()
	for (i in 1:reps){
		samps = sample(1:nrow(x),replace=T)[1:(nrow(x)*bootstrap_frac)]
		currx = x[samps,]
		curry = y[samps]
		models[[i]] = liblinear_wrapper(currx,curry,cross=cross,scale=F,tryTypes=tryTypes,tryCosts=tryCosts)
	}
	obj = list(models = models,bootstrap_frac=bootstrap_frac,scale=scale,means=means,sds=sds)
	class(obj) = "liblinear_svm_wrapper"
	return (obj)
}

predict.liblinear_svm_wrapper<-function(obj,newx){
	preds = matrix(nrow=nrow(newx),ncol=2)
	if (obj$scale){
		for (j in 1:ncol(newx)){
			newx[,j] = (newx[,j] - obj$means[j])/obj$sds[j]
		}
	}
	return (predict(obj$model,newx,proba=T)$probabilities)
}
#############################################################

# This algorithm implements a simple multitask correction 
# similar to Groves et al. 2012
# x is the data matrix (columns are features)
# y is a {0,1} matrix in which each column represent a different task
# When the base classifiers are given from an external analysis, we use them
# instead of running the training again. In these cases cvVec must be given
# We use internal cross validation to make the chain analysis 
# independent of the binary relevance predictions.
# If the number of folds is set to 1 then we use the same data for
# base classifier (the binary relevance models) and the chain rule learning.
multitaskChainClassifier<-function(x,y,
		baseClassifier=featureSelectionClassifier,secondClassifier=randomForest,
		baseClassifierArgs=list(),secondClassifierArgs=list(),
		baseClassifierPredArgs=list(type="prob"),secondClassifierPredArgs=list(type="prob"),
		isFirstClassifierSequential = T,isSecondClassifierBinaryRelevance = T,addOriginalfeaturesToSeconfClassifier=F,
			cvfolds=1,baseClassifiers = NULL,cvVec=NULL,d=NULL){
	
	if (is.null(cvVec)){cvVec = getStratifiedMultitaskCVDOStructure(y,cvfolds,d=d)}
	# for each fold we have a base classifier and a second classifiers
	trainBaseClassifiers = F
	if (is.null(baseClassifiers)){
		baseClassifiers = list()
		trainBaseClassifiers = T
	}
	basePredictions = list()
	secondClassifiers = NULL
	if (!is.null(secondClassifier)){secondClassifiers = list()}
	for (fold in 1:cvfolds){
		basex = x[-which(cvVec==fold),]; basey = y[-which(cvVec==fold),]
		if (length(basey)==0){basex = x;basey = y}
		# learn the base classifiers
		if (trainBaseClassifiers){
			baseClassifiers[[fold]] = apply(basey,2,getBaseClassifier,x=basex,func=baseClassifier,baseClassifierArgs)
		}
		# learn the second classifiers
		secondx = x[which(cvVec==fold),];secondy = y[which(cvVec==fold),]
		if (length(secondx)==0){secondx=basex;secondy=basey}
		if (!is.null(secondClassifier)){
			basePredictions[[fold]] = sapply(baseClassifiers[[fold]],getPredictionVector,x=secondx,baseClassifierPredArgs)
			rownames(basePredictions[[fold]]) = rownames(secondx)
			if (addOriginalfeaturesToSeconfClassifier){
				basePredictions[[fold]] = cbind(secondx,basePredictions[[fold]])
			}
			if (!isSecondClassifierBinaryRelevance){
				currargs = c(list(basePredictions[[fold]],secondy),secondClassifierArgs)
				secondClassifiers[[fold]] = do.call(secondClassifier,args=currargs)
			}
			else{secondClassifiers[[fold]] = apply(secondy,2,getBaseClassifier,x=basePredictions[[fold]],func=secondClassifier,secondClassifierArgs)}
		}
	}
	obj = list(
			baseClassifiers=baseClassifiers,basePredictions=basePredictions,
			secondClassifiers=secondClassifiers,
			cvfolds=cvfolds,cvVec = cvVec ,
			baseClassifierPredArgs=baseClassifierPredArgs,
			secondClassifierPredArgs=secondClassifierPredArgs,
			isSecondClassifierBinaryRelevance = isSecondClassifierBinaryRelevance,
			isFirstClassifierSequential = isFirstClassifierSequential,
			addOriginalfeaturesToSeconfClassifier = addOriginalfeaturesToSeconfClassifier
		    )
	
	class(obj) <- "multitaskChainClassifier"
	return (obj)
}
predict.multitaskChainClassifier<-function(obj,newx){
	if (is.null(dim(newx))){newx = matrix(newx,nrow=1)}
	numclasses = length(obj$baseClassifiers[[1]])
	numsamples = dim(newx)[1]
	cvpreds = matrix(rep(0,times=numclasses*numsamples),nrow=numsamples)
	rownames(cvpreds) = rownames(newx)
	colnames(cvpreds) = names(obj$baseClassifiers[[1]])
	for (fold in 1:obj$cvfolds){
		basePredictions = sapply(obj$baseClassifiers[[fold]],getPredictionVector,x=newx,obj$baseClassifierPredArgs)
		if (is.null(dim(basePredictions))){
			basePredictions = matrix(basePredictions,nrow=1)
			colnames(basePredictions) = names(obj$baseClassifiers[[fold]])
		}
		rownames(basePredictions)=rownames(newx)
		print (dim(basePredictions))
		if (dim(basePredictions)[1] != dim(newx)[1]){stop("base predictor failure: not all samples were used")} 
		if (dim(basePredictions)[2] != length(obj$baseClassifiers[[fold]])){stop("base predictor failure: not all predictors were used")}
		if (is.null(obj$secondClassifiers)){cvpreds = cvpreds + basePredictions; next}
		if (obj$addOriginalfeaturesToSeconfClassifier){basePredictions = cbind(newx,basePredictions)}
		if (obj$isSecondClassifierBinaryRelevance){cvpreds = cvpreds + sapply(obj$secondClassifiers[[fold]],getPredictionVector,x=basePredictions,obj$secondClassifierPredArgs)} 
		else{
			curr_pred =  predict(obj$secondClassifiers[[fold]],basePredictions)
			cvpreds = cvpreds + curr_pred
		}
	}
	return (cvpreds/obj$cvfolds)
}
#############################################################

# An implementation of simple BR analysis
simple_binary_relevance<-function(x,y,
		baseClassifier=featureSelectionClassifier,baseClassifierArgs=list(),baseClassifierPredArgs=list(type="prob")){
	
	# learn the base classifiers
	baseClassifiers = list()
	baseClassifiers[[1]] = list()
	for (j in 1:ncol(y)){
		print(colnames(y)[j])
		baseClassifiers[[1]][[colnames(y)[j]]] = getBaseClassifier(y[,j],x=x,func=baseClassifier,baseClassifierArgs)
		gc()
	}
	obj = list(
			secondClassifiers = NULL,
			baseClassifiers=baseClassifiers,cvfolds=1,
			baseClassifierPredArgs=baseClassifierPredArgs
		    )
	class(obj) <- "simple_binary_relevance"
	return (obj)
}

predict.simple_binary_relevance<-function(obj,newx){
	print (paste("in prediction",obj$cvfolds))
	if (is.null(dim(newx))){newx = matrix(newx,nrow=1)}
	numclasses = length(obj$baseClassifiers[[1]])
	numsamples = dim(newx)[1]
	cvpreds = matrix(rep(0,times=numclasses*numsamples),nrow=numsamples)
	rownames(cvpreds) = rownames(newx)
	colnames(cvpreds) = names(obj$baseClassifiers[[1]])
	basePredictions = sapply(obj$baseClassifiers[[1]],getPredictionVector,x=newx,obj$baseClassifierPredArgs)
	if (is.null(dim(basePredictions))){
		basePredictions = matrix(basePredictions,nrow=1)
		colnames(basePredictions) = names(obj$baseClassifiers[[1]])
	}
	rownames(basePredictions)=rownames(newx)
	print (dim(basePredictions))
	if (dim(basePredictions)[1] != dim(newx)[1]){stop("base predictor failure: not all samples were used")} 
	if (dim(basePredictions)[2] != length(obj$baseClassifiers[[1]])){stop("base predictor failure: not all predictors were used")}
	cvpreds = cvpreds + basePredictions
	return (cvpreds)
}
#############################################################

# This is a simple transformation method that treats multilabel 
# classification as a multiclass classification task.
# It uses a base multiclass classifier to obtain a predictive model.
# For prediction of a label v we sum over the predictions of all classes
# that contain samples of v.
multitaskViaMulticlassClassifier<-function(x,y,
		baseClassifier=randomForest,baseClassifierArgs=list(),baseClassifierPredArgs=list(type="prob"),useFS=F,numFeatures=50,d=NULL){

	print (useFS)
	features = colnames(x)
	if (useFS){
		features = c()
		for (j in 1:ncol(y)){
			currscores = calcRocOrAuprScores(x,y[,j],useROC=T,datasets=d)
			currfeatures = names(sort(currscores,decreasing=T)[1:numFeatures])
			features = unique(c(features,currfeatures))
		}
		print (paste('total number of features after feature selection: ',length(features)))
	}
	original_tasks = colnames(y)
	y_labels = apply(y,1,paste,collapse=',',sep='')
	y_labs_factor = as.factor(y_labels)
	print (paste('total number of classes: ',length(unique(y_labels))))
	args = c(list(x=x[,features],y=y_labs_factor),baseClassifierArgs)
      model = do.call(baseClassifier,args=args)
	print (class(model))
	obj = list(features=features,model=model,baseClassifierPredArgs=baseClassifierPredArgs,y=y,y_labels=y_labels,y_labs_factor=y_labs_factor,original_tasks=original_tasks)
	class(obj) <- "multitaskViaMulticlassClassifier"
	return (obj)
}

predict.multitaskViaMulticlassClassifier<-function(obj,newx,...){
	if (is.null(dim(newx))){newx = matrix(newx,nrow=1)}
	args = c(list(obj$model,newx[,obj$features]),obj$baseClassifierPredArgs,...)
	# predObj is a matrix: rows are samples 
      predObj = do.call(predict,args=args)
	if (class(obj$model)=='svm'){predObj = attr(predObj, "probabilities")}
	preds = t(apply(predObj,1,multitaskViaMulticlassClassifier_transformPrediction,xnames=colnames(predObj)))
	print (dim(preds))
	# fix in case we predict a single sample
	if (is.null(dim(preds))){preds = matrix(preds,nrow=1)}
	rownames(preds) = rownames(newx)
	colnames(preds) = obj$original_tasks
	return (preds)
}

multitaskViaMulticlassClassifier_transformPrediction<-function(x,xnames){
	names_matrix = sapply(xnames,function(x)strsplit(x,split=',')[[1]])
	mode(names_matrix) = 'numeric'
	return(names_matrix%*%x)
}


#############################################################

# This is a fast implementation of the Bayesian correction methos used by previous
# papers (e.g., Huang et al., 2010 PNAS).
# The method is similar to the classifier chains algorithm in that it uses a base
# classifier to obtain initial classification results.
# Unlike classifier chains, the different labels are combined using a known Bayesian network
# that represents the dependencies among the labels.
# This method makes sense when the network is known. For example for tasks involving ontologies 
# (e.g., gene or disease ontology)
#revArcs=F;nbins=5;minEntryCount=1
bayesianCorrectionClassifier_new<-function(basePredictions,y,revArcs=F,nbins=5,minEntryCount=1){
	doids = colnames(basePredictions)
	bn_obj = getBayesianNetwork(doids,show=F,revArcs=revArcs)
	initialBN = bn_obj[[1]]
	new_colnames = bn_obj[[2]]
	colnames(basePredictions) = new_colnames[colnames(basePredictions)]
	# discretization of the base predictions
	breaks = list()
	for(j in 1:ncol(basePredictions)){
		breaks[[j]] = find_breaks(basePredictions[,j],numbins=nbins,minEntryCount=minEntryCount)
	}
	basePredictions = apply(basePredictions,2,discPredVector,numbins=nbins,minEntryCount=minEntryCount)
	dataForBN_preds = data.frame(basePredictions,check.names = F,stringsAsFactors=T)
	dataForBN_real = y
	colnames(dataForBN_real) = new_colnames[colnames(dataForBN_real)]
	colnames(dataForBN_preds) = paste(colnames(dataForBN_preds),"_hat",sep="")
	# At this point our generated network might have added some
	# additional classes that are not represented in y
	# by definition, these must be classes that represent either a tissue
	# or a DOID (but not together)
	net_nodes = nodes(initialBN)
	unrepresented_nodes = setdiff(net_nodes,c(colnames(dataForBN_preds),colnames(dataForBN_real)))
	mode(dataForBN_real)='numeric'
	for(nodename1 in unrepresented_nodes){
		nodename = names(new_colnames)[new_colnames==nodename1]
		curr_cols_in_y = which(grepl(colnames(y),pattern=nodename))
		if(length(curr_cols_in_y)==1){
			currv = y[,curr_cols_in_y]
		}
		else{
			currv = as.numeric(apply(y[,curr_cols_in_y],1,sum)>0.5)
		}
		dataForBN_real = cbind(dataForBN_real,currv)
		colnames(dataForBN_real)[ncol(dataForBN_real)] = nodename1
	}
	mode(dataForBN_real)='character'
	dataForBN_real = data.frame(dataForBN_real,check.names = F,stringsAsFactors=T)
	# make sure that all prediction columns have at least one level
	for(j in 1:ncol(dataForBN_preds)){
		if(length(levels(dataForBN_preds[,j]))==1){
			levels(dataForBN_preds[,j]) = c(levels(dataForBN_preds[,j]),"dummy level (no variance in training data)")
		}
	}
	# make sure that all real class columns have "0" and "1" as labels
	for(j in 1:ncol(dataForBN_real)){
		if(length(levels(dataForBN_real[,j]))==1){
			levels(dataForBN_real[,j]) = c(levels(dataForBN_real[,j]),setdiff(c("0","1"),levels(dataForBN_real[,j])))
		}
	}
	dataForBN = cbind(dataForBN_preds,dataForBN_real)
	fittedBN = bn.fit(initialBN,dataForBN)
	grainNet = NULL
	try({grainNet = as.grain(fittedBN)})
	# test
	#  get_no_hat_nodes_prediction_vector(dataForBN_preds[20,],fittedBN=fittedBN,non_hat_nodes=setdiff(nodes(fittedBN),netNames))
	hat_cloumns = nodes(fittedBN)[grepl(nodes(fittedBN),pattern="_hat$")]
	training_data_levels = lapply(dataForBN_preds,levels)
	obj = list(fittedBN=fittedBN,hat_cloumns = hat_cloumns,training_data_levels=training_data_levels, new_colnames=new_colnames,
			nbins=nbins,minEntryCount=minEntryCount,breaks=breaks,tr_node_names = colnames(dataForBN),grainNet=grainNet)
	class(obj) <- "bayesianCorrectionClassifier_new"
	return (obj)
}

# numbins is the maximal number of bins in the data
# A simple discretization algorithm: 
# Split the interval of y into numbins equal inteval
# omit intervals whose count is too low
# merge the resulting intervals to cover the enrite
# 0,1 interval.
find_breaks<-function(y,numbins=5,minEntryCount=1){
	if(sd(y)==0){return(c(-0.0001,1.0001))}
	# get the initial breaks, this is from the source code of cut
	nb <- as.integer(numbins + 1) # one more than #{intervals}
	dx <- diff(rx <- range(y, na.rm = TRUE))
	if(dx == 0) {
		dx <- abs(rx[1L])
		initial_breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000,length.out = nb)
	} 
	else {
		initial_breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
		initial_breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + dx/1000)
	}
	initial_breaks = round(initial_breaks,digits=4)
	##
	biny = cut(y,breaks = initial_breaks,dig.lab=4)
	counts = table(biny)
	breaks = names(counts)
	breaks = gsub(breaks,pattern="\\(|\\]",replace="")
	breaks = strsplit(breaks,split=',')
	breaks = sapply(breaks,as.numeric)
	to_keep = counts>=minEntryCount
	breaks = breaks[,to_keep]
	counts = counts[to_keep]
	for(j in 2:ncol(breaks)){
		if(abs(breaks[1,j] - breaks[2,j-1])>1e-10){
			breaks[2,j-1] = breaks[1,j]
		}
	}
	breaks = sort(unique(c(breaks)))
	if (breaks[1] > 0 ){breaks[1] = -0.0001}
	if (breaks[length(breaks)]<1.0001){ breaks[length(breaks)]=1.0001}
	return(breaks)
}

# Important: the labels of x are kept as the upper bounds
# of the intervals.
discPredVector<-function(x,numbins=5,minEntryCount=1){
	x[is.na(x)] = 0;x[is.nan(x)] = 0
	x[is.infinite(x) & x>0] = max(x[!is.infinite(x)])
	x[is.infinite(x) & x<0] = min(x[!is.infinite(x)])
	breaks=find_breaks(x,numbins,minEntryCount)
	y = cut(x,breaks=breaks,breaks,ordered_result=T,dig.lab=4,
		labels = paste("<",breaks[-1]))
	if(length(levels(y))==1){
		levels(y) = c(levels(y),"dummy level (no variance in training data)")
	}
	return (y)
}

# v is a discrete vector of values of the _hat vectors
# this function returns the predictions for the non _hat nodes
# the predictions are the probabilities to observe "1"
# the last argument "non_hat_nodes" is also given to keep the correct
# order of the nodes
getEventForBayesianNet<-function(variables,values,returnString=T){
	variables = unlist(variables);values = unlist(values)
	return (paste("(", variables, "==","'",values, "'",")", sep = "",collapse = " & "))
}
get_no_hat_nodes_prediction_vector<-function(v,fittedBN,non_hat_nodes,method='cpquery-ls',...){
	preds = c()
	for(n in non_hat_nodes){
		preds[n] = 0
		try({
			# Currently not working (bnlearn version 4.0) - sometimes leads to an R failure
			if(method=='predict'){
				pred = predict(fittedBN,node = n,data=v,method='bayes-lw',prob=T,...)
				preds[n] = attr(pred,"prob")["1",1]
			}
			if(method == 'cpquery-ls'){
				evid <<- getEventForBayesianNet(names(v),v)
				evnt <<- getEventForBayesianNet(n,1)
				preds[n] = cpquery(fittedBN,event=eval(parse(text=evnt)),evidence=eval(parse(text = evid)),...)
			}
			# Currently not working (bnlearn version 4.0) - leads to an R failure
			if(method == 'cpquery-lw'){
				evnt <<- getEventForBayesianNet(n,1)
				preds[n] = cpquery(fittedBN,event=eval(parse(text=evnt)),evidence=as.list(v),method='lw',...)
			}
		})
	}
	return (preds)
}

get_grain_prediction<-function(node,evid,grain_net){
	pred = querygrain(grain_net,nodes=c(node,names(evid)),evidence=evid,type="conditional")[[1]]["1"]
	return(pred)
}

predict.bayesianCorrectionClassifier_new<-function(obj,basePredictions,...){
	if (is.null(dim(basePredictions))){
		currnames = names(basePredictions)
		basePredictions = matrix(basePredictions,nrow=1)
		colnames(basePredictions) = currnames
	}
	input_rownames = rownames(basePredictions);input_colnames = colnames(basePredictions)
	colnames(basePredictions) = obj$new_colnames[colnames(basePredictions)]
	hat_names = paste(colnames(basePredictions),'_hat',sep='')
	colnames(basePredictions) = hat_names
	non_hat_nodes = obj$new_colnames[input_colnames]
	# discretize based on the training data
	m = list()
	for(j in 1:ncol(basePredictions)){
		currCol = as.character(cut(basePredictions[,j],breaks = obj$breaks[[j]],labels = paste("<",obj$breaks[[j]][-1])))
		currCol = factor(currCol,levels = obj$training_data_levels[[j]])
		names(currCol) = rownames(basePredictions)
		m[[colnames(basePredictions)[j]]] = currCol
	}
	basePredictions = data.frame(m,check.names = F,stringsAsFactors=F)
	#check_bc_model_vs_data(obj,basePredictions)
	preds = c()
	if(is.null(obj$grainNet)){
		for(i in 1:nrow(basePredictions)){
			currp = get_no_hat_nodes_prediction_vector(basePredictions[i,],fittedBN=obj$fittedBN,non_hat_nodes=non_hat_nodes,...)
			preds = rbind(preds,currp)
			if (i%%100 == 0){print (i)}
		}
	}
	else{
		# The test data might be too large and the predict method
		# might require too much space and fail.
		# send at most 500 samples at a time
		chunks = c();last_chunk=0
		while (length(chunks) < nrow(basePredictions)){
			last_chunk = last_chunk+1
			currv = rep(last_chunk,500)
			chunks = c(chunks,currv)
		}
		chunks = chunks[1:nrow(basePredictions)]
		for(j in sort(unique(chunks))){
			inds = which(chunks==j)
			preds_list = predict(obj$grainNet,non_hat_nodes,colnames(basePredictions),basePredictions[inds,], type="distribution")
			currpreds = c()
			for(nodename in non_hat_nodes){
				currpreds = cbind(currpreds,preds_list[[1]][[nodename]][,"1"])
			}
			colnames(currpreds) = non_hat_nodes
			preds = rbind(preds,currpreds)
			gc()
		}
	}

	# make sure that the output is a matrix
	if (is.null(dim(preds))){
		currnames = names(preds);preds = matrix(preds ,nrow=1)
		colnames(preds) = currnames
	}
	rownames(preds) = input_rownames
	colnames(preds) = input_colnames
	print ('DONE with Bayesian network-based prediction')
	return (preds)
}
# This method receives as input a set of class names and returns the
# corresponding network structure.
# The name of each class is split by semicolon and the resulting array specifies:
# 	entry 1: the DOID or control (always there)
#	if the array length is >1 then the last entry specifies the tissue
#	all doids and tissues are taken as parent terms in the network
#	the class names as is are taken as the "initial" leaves
#	each initial leaf has a son with the same name and a "_hat" suffix - these nodes
#	are the nodes that specify the actual discrete predictions 
# We add four types of arcs:
#	(1) From a class name to its "_hat" son
#	(2) DO arcs within tissues
#	(3) An arc from a tissue to its: control class, top ancestors in its tissue DAG
#	(4) Arcs between "clean" do terms: terms that either appear as classes or are represented in
#		multiple tissues
# Dec 2016: a new type of arcs was added: "_star" node to its parent.
#	a "_star" node is dummy nodes that is added to a term D and represents samples
#	of D that are not in any of D's children.
#	If classes_names contains _star nodes then these are added with their arcs
#	and then removed as these are not required for downstream analyses.
#	We ASSUME that for each _star node another term in classes_names represents its parent.
getBayesianNetwork<-function(class_names,addControlArcs=F,show=F,revArcs=F,splitPattern=';'){
	class_names = gsub(class_names,pattern="Control",replace="control")
	class_names2node_names = paste("Node",1:length(class_names),sep='')
	names(class_names2node_names) = class_names
	node_counter = length(class_names)
	nodes = c();edges = c()
	# (1) add the class name nodes and their "hat" nodes
	for (d in class_names){
		nodes = c(nodes,paste("[",class_names2node_names[d],"]",sep=""))
		nodes = c(nodes,paste("[",class_names2node_names[d],"_hat]",sep=""))
		edges = c(edges,paste("[",class_names2node_names[d],"_hat|",class_names2node_names[d],"]",sep=""))
	}
	# (1.1) Deal with the "star" nodes
	star_nodes = class_names[grepl(class_names,pattern="_star$",perl=T)]
	if (length(star_nodes)>0){
		for (d in star_nodes){
			p = gsub(d,pattern = "_star$",perl=T,replace='')
			edges = c(edges,paste("[",class_names2node_names[d],"|",class_names2node_names[p],"]",sep=""))
		}
	}
	class_names = setdiff(class_names,star_nodes)
	# Split the current class names
	class2arr=list();tissues = c()
	for (d in class_names){
		arr = strsplit(d,split=splitPattern)[[1]]
		class2arr[[d]]=arr
	}
	# (2) Go over all doids, within each tissue add all arcs of the DO
	parents = as.list(DOPARENTS);ancs = as.list(DOANCESTOR)
	doids = unique(sapply(class2arr,function(x)x[1]))
	doids = setdiff(doids,"control")
	# Go over all doids, add all arcs from a parent term to 
	# a child term if they have the same tissue
	for (d in doids){
		# Get all class names of the current doid that have tissue
		nodes_of_d = setdiff(class_names[grepl(pattern=paste(d,splitPattern,sep=''),class_names)],d)
		for (d1 in nodes_of_d){
			# we get all ancestory classes - ancestor doid, same tissue
			curr_arr = class2arr[[d1]]
			d1_tissue = curr_arr[length(curr_arr)]
			curr_ps = ancs[[d]]
			curr_ps = intersect(curr_ps,doids)
			# keep curr ancestors that have a class and have the same tissue
			new_ps_set = c()
			for(p in curr_ps){
				nodes_of_p = setdiff(class_names[grepl(pattern=paste(p,splitPattern,sep=''),class_names)],p)
				for (pclasses in nodes_of_p){
					curr_p_arr = class2arr[[pclasses]]
					if(length(curr_p_arr)==1){next}
					same_tissue = d1_tissue == curr_p_arr[length(curr_p_arr)]
					if(length(nodes_of_p)>0 && same_tissue){
						new_ps_set = union(new_ps_set,p)
					}
				}
			}
			curr_ps = new_ps_set
			if(length(curr_ps)>1){curr_ps = clean_ancestor_list(curr_ps,ancs)}
			# here we have a compact set of ancestors: add the arcs
			if (is.null(curr_ps)|length(curr_ps)==0){next}
			for (p in curr_ps){
				# We now have a doid d and a parent in the network p
				# Add all edges from a term with d to p with the same tissue
				# Also, add a direct edge from d to p (without any tissue)
				nodes_of_d = setdiff(class_names[grepl(pattern=paste(d,splitPattern,sep=''),class_names)],d)
				nodes_of_p = setdiff(class_names[grepl(pattern=paste(p,splitPattern,sep=''),class_names)],p)
				for (d1 in nodes_of_d){
					curr_arr = class2arr[[d1]]
					d1_tissue = curr_arr[length(curr_arr)]
					for (p1 in nodes_of_p){
						curr_arr = class2arr[[p1]]
						p1_tissue = curr_arr[length(curr_arr)]
						if(p1_tissue!=d1_tissue){next}
						edges = c(edges,paste("[",class_names2node_names[d1],"|",class_names2node_names[p1],"]",sep=""))
					}
				}
			}
		}
	}
	# (3) Add tissue nodes as parrents
	# add links to the DO tree top ancestros and to the control
	tissues = c()
	for (d in class_names){
		arr = class2arr[[d]]
		if(length(arr)==1){next}
		tissue = arr[length(arr)]
		if(!is.element(tissue,set=names(class_names2node_names))){
			node_counter = node_counter+1
			class_names2node_names[tissue] = paste("Node",node_counter,sep='')
			nodes = c(nodes,paste("[",class_names2node_names[tissue],"]",sep=""))
			tissues = union(tissues,tissue)
		}
	}
	for(tissue in tissues){
		all_classes = grepl(class_names,pattern=paste(';',tissue,sep=''))
		all_classes = class_names[all_classes]
		control_node = all_classes[grepl(all_classes,pattern='control;')]
		if(length(control_node)>0){
			edges = c(edges,paste("[",class_names2node_names[control_node],"|",class_names2node_names[tissue],"]",sep=""))
			all_classes = setdiff(all_classes,control_node)
		}
		curr_doids = sapply(class2arr[all_classes],function(x)x[[1]])
		curr_ancs = names(curr_doids)
		if(length(curr_doids)>1){
			curr_ancs = names(clean_ancestor_list(curr_doids,ancs,remvAncestors=F))
		}
		for(curranc in curr_ancs){
			edges = c(edges,paste("[",class_names2node_names[curranc],"|",class_names2node_names[tissue],"]",sep=""))
		}
	}
	# (4) Add the DO backbone
	arr_lengths = sapply(class2arr,length)
	single_dos_as_classes = names(class2arr)[arr_lengths==1]
	dos_class2arr = class2arr[sapply(class2arr,function(x)grepl(x[1],pattern="DOID")&length(x)>1)]
	dos_class2arr_dos = sapply(dos_class2arr,function(x)x[[1]])
	multiple_class_dos = table(dos_class2arr_dos)
	multiple_class_dos = names(multiple_class_dos)[multiple_class_dos>1]
	doids = union(single_dos_as_classes,multiple_class_dos)
	# (4.1)
	# Go over all class names.
	# If they represent a DO from the doids set then add the doid
	# as a node (if it is not represented).
	# Add an arc from its specific tissue nodes to it
	for (d in class_names){
		arr = class2arr[[d]]
		doid = arr[1]
		if(!is.element(doid,set=doids)){next}
		if(!is.element(doid,set=names(class_names2node_names))){
			node_counter = node_counter+1
			class_names2node_names[doid] = paste("Node",node_counter,sep='')
			nodes = c(nodes,paste("[",class_names2node_names[doid],"]",sep=""))
		}
		# Add an arc from the current name to its doid
		if(doid != d){
			edges = c(edges,paste("[",class_names2node_names[d],"|",class_names2node_names[doid],"]",sep=""))
		}
	}
	# (4.2)
	# Add all arcs between the doids
	for (d in doids){
		curr_ps = parents[[d]]
		curr_ps = intersect(curr_ps,doids)
		# if there is no parent then add all arcs from their lowest present ancestors
		if (is.null(curr_ps)|length(curr_ps)==0){
			curr_ps = ancs[[d]]
			curr_ps = intersect(curr_ps,doids)
			if(length(curr_ps)>1){curr_ps = clean_ancestor_list(curr_ps,ancs)}
		}
		if (is.null(curr_ps)|length(curr_ps)==0){next}
		for (p in curr_ps){
			edges = c(edges,paste("[",class_names2node_names[d],"|",class_names2node_names[p],"]",sep=""))
		}
	}
	nodes = unique(nodes);edges=unique(edges)
	# Go over all doids that 
	str_rep = paste(paste(nodes,collapse=''),paste(edges,collapse=''),collapse='')
	node_names = gsub(nodes,pattern='\\[|\\]',replace='')
	res = empty.graph(node_names)
	modelstring(res) = str_rep
	if (show){plot(res)}
	return (list(node = res,class_names2node_names=class_names2node_names,tissues=tissues))
}

bn_get_arcs_as_matrix<-function(arcs,mapp){
	arcs = gsub(arcs,pattern='\\[|\\]',replace='')
	arcs = t(sapply(arcs,function(x)strsplit(x,split='\\|')[[1]]))
	arcs = cbind(arcs[,2],arcs[,1]);colnames(arcs)=c("from","to")
	return(arcs)
}
# Do not run
QA_bn_find_cycles<-function(edges,class_names2node_names){
	arcs = bn_get_arcs_as_matrix(edges)
	ig = igraph::graph_from_edgelist(arcs, directed = TRUE)
	is_dag(ig)
	has_cycle = setdiff(names(V(ig)),names(topo_sort(ig)))
	subg_arcs = arcs[is.element(arcs[,1],set=has_cycle) & is.element(arcs[,2],set=has_cycle),]
	subg_ig = igraph::graph_from_edgelist(subg_arcs, directed = TRUE)
	plot(subg_ig)
	is_dag(subg_ig)
	#currnodes = c("Node46","Node49","Node50")
	revmap = names(class_names2node_names);names(revmap)=class_names2node_names
	#d = revmap[currnodes][1]
	#d = class2arr[[d]][1]
	return(revmap[has_cycle])
}
bn_map_non_hat_edges_to_names<-function(edges,class_names2node_names){
	revmap = names(class_names2node_names);names(revmap)=class_names2node_names
	arcs = bn_get_arcs_as_matrix(edges)
	arcs = arcs[!grepl(arcs[,2],pattern="_hat"),]
	newarcs = cbind(revmap[arcs[,1]],revmap[arcs[,2]])
	return(newarcs)
}
get_graph_for_cytoscape<-function(arcs,mapping){
	revmap = names(mapping);names(revmap) = mapping
	arcs = arcs[!grepl(arcs[,2],pattern="_hat"),]
	newarcs = cbind(revmap[arcs[,1]],revmap[arcs[,2]])
	return(newarcs)
}

# Tests
#library(DO.db)
#class_names = c("DOID:1612;breast","DOID:10283;prostate","control;prostate")
#getBayesianNetwork(doids)

clean_ancestor_list<-function(s,doancs,remvAncestors=T){
	to_rem = c()
	n = length(s)
	for(i in 1:n){
		do1 = s[i]
		for(j in 1:n){
			do2 = s[j]
			if(i==j){next}
			#print(class(do2))
			if(!is.element(do1,set=c(doancs[[do2]]))){next}
			#print("$$$")
			# do1 is an ancestor of do2: remove according remvAncestors
			if (remvAncestors){to_rem = c(to_rem,i)}
			else{to_rem = c(to_rem,j)}
		}
	}
	to_rem = unique(to_rem)
	return(s[-to_rem])
}

###################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################



getBaseClassifier<-function(y,x,func,...){
	if (length(unique(y))<2){stop("cannot create base classifer: less than two classes!")}
	# print (sum(table(y,as.factor(y))) - sum(diag(table(y,as.factor(y)))))
	args = c(list(x=x,y=as.factor(y)),...)
	return (do.call(func,args=args))
}

getPredictionVector<-function(obj,x,...){
	args = c(list(obj,x),...)
	predObj = do.call(predict,args=args)
	predsMat = getPredProbabilities(predObj)
	# in case we do not use the probabilities (e.g., we use the actual predictions)
	if (is.null(dim(predsMat)) || is.null(predsMat)){return(predObj)}
	# Deal with the prediction matrix in case we want the actual probabilities
	if (mode(predsMat)!='numeric'){return(predsMat)}
	currClassInd = which(colnames(predsMat)=="1")
	currVec = predsMat[,currClassInd]
	names(currVec) = rownames(predsMat)
	#print (currVec)
	return (currVec)
}

getPredProbabilities <- function(pred_obj){
	probs = attr(pred_obj,"prob")
	if (!is.null(probs)){return(probs)}
	probs = attr(pred_obj,"probabilities")
	if (!is.null(probs)){return(probs)}
	if (is.element("matrix",class(pred_obj))){return (pred_obj)}
	if (length(pred_obj)>=2) {probs = pred_obj[[2]]}
	if (!is.null(probs)){return(probs)}
	return (NULL)
}
#############################################################

# Methods for getting a stratified partitioning of samples into folds.

# Method 1: a simple analysis for a single vector y
# y is categorical vector
# return value is a matrix in which each row is
# a sample, and each column is a fold
# RETij is true iff sample i is a test
# in fold j
getStratifiedCV<-function(y,folds){
	cvVec = rep(1,times = length(y))
	names(cvVec) = names(y)
	for (cl in unique(y)){
		curr_inds = which(y==cl)
		curr_inds_part = partitionToK(curr_inds,folds)
		cvVec[curr_inds]=curr_inds_part
	}
	return (cvVec)
}
partitionToK<-function(v,k){
	x = sample(v)
	assignment = rep(1:k,times = floor(length(v)/k)+1)
	assignment = assignment[1:length(x)]
	names(assignment)=x
	return (assignment)
}

# Method 2: a more involved analysis that considers the 
# partitioning of the samples by their studies.
# y is a matrix, each column is a task
# return value is a matrix in which each row is
# a sample, and each column is a fold
# RETij is true iff sample i is a test in fold j
getStratifiedMultitaskCVDOStructure<-function(y,folds,d=NULL){
	# Deal with class overlap, assign samples to their new class:
	# the set of all do terms of the sample
	newClasses = c()
	for (i in 1:dim(y)[1]){
		newClasses[i] = paste(y[i,],collapse=';')
	}
	# run standard stratified CV on the new vector
	if (is.null(d)|| folds <2){return (getStratifiedCV(newClasses,folds))}
	# create a cv vector when a mapping to datasets is given
	d = d[rownames(y)]
	# create a vector of bins for the different diseases
	datasets = unique(d)
	datasets = sample(datasets)
	# for each fold keep the number of samples
	# and the number of datasets for each disease
	fold_samples = matrix(nrow=folds,ncol=ncol(y))
	fold_samples[,]=0
	fold_datasets = matrix(nrow=folds,ncol=ncol(y))
	fold_datasets[,] = 0
	cvVec = rep(1,nrow(y))
	names(cvVec) = rownames(y)
	for (currd in datasets){
		# current sample inds
		curr_inds = which(d==currd)
		# get y's submatrix that corresponds to these samples
		curr_sub_y = y[curr_inds,]
		# the classes for which there are positive samples
		curr_contributions = which(apply(curr_sub_y,2,sum)>0)
		# add current datasets to the next fold 
		# select the fold with minimal number of datasets
		# get the fold_datasets submatrix - the columns of
		# the classes in which the current datasets has positive samples
		curr_contr_datasets = matrix(fold_datasets[,curr_contributions],nrow=folds)
		min_cell = min(curr_contr_datasets)
		# this is the fold for which the new data will be added
		curr_fold = which(apply(curr_contr_datasets==min_cell,1,any))[1]
		# update the relevant data structures
		cvVec[curr_inds] = curr_fold
		fold_samples[curr_fold,] = fold_samples[curr_fold,]+apply(curr_sub_y,2,sum)
		fold_datasets[curr_fold,curr_contributions] = fold_datasets[curr_fold,curr_contributions]+1
	}
	#print (fold_datasets)
	return (cvVec)
}

#getMostSpecificDOs<-function(termslist,anceslists){
#	if (length(termslist)==0){return (c())}
#	for (currdo in termslist){termslist = removeAncestorsFromSet(currdo,anceslists,termslist)}
#	return (termslist)
#}
#removeAncestorsFromSet<-function(go,ancslist,s){
#	return (setdiff(s,ancslist[[go]]))
#}

#############################################################
# This is an important function that runs cross-validation analysis.
# x,y: feature and label matrices. Rows are samples.
# batch_info: partitioning of samples to folds. In our analyses this vector
# should be obtained by considering the mapping of samples to their studies.
# We want to avoid splitting a study between folds.
# classification_function - the multilabel classification function to use
# numCores - in unix we can use many CPUs (e.g., one for each fold).
# otherRunRes: the models from the result of another LDO run using the same batch_info vector
# setting this to an object (i.e., not NULL) means that the base classifiers and the internal CV split are taken from this object
getMultitaskLeaveDatasetOutClassificationPerformance_multicore<-function(
	x,y,batch_info,classification_function=multitaskChainClassifier,prediction_args=list(),numCores = 1,otherRunRes=NULL,...){

	# remove classes that are represented in at most a single fold only
	to_r = which(apply(y,2,sum)==0)
	for (j in 1:dim(y)[2]){
		curr_table = table(y[,j],batch_info)
		num_batches_for_each_class = apply(curr_table!=0,1,sum)
		# a class we cannot test: the class or its complement is not represented in at least two batches
		if (length(which(num_batches_for_each_class<=1))>0){to_r=c(to_r,j)}
	}
	if (length(to_r)>0){y = y[,-to_r]}
	print (paste("number of removed classes:",length(to_r)))
	batches = unique(batch_info)
	print ("calling the CV process")
	preds_obj_res = mclapply(batches,applyMultitaskClassificationTest,mc.cores=numCores,x=x,y=y,batch_info=batch_info,
						classification_function=classification_function,prediction_args=prediction_args,otherRunRes = otherRunRes,...)
	names(preds_obj_res) = batches
	preds_lists = list()
	models = list(preds_obj_res)
	for (currn in batches){
		if (class(preds_obj_res[[currn]])=="try-error"){
			print (currn)
			print (preds_obj_res[[currn]]);stop()
		}
		preds_lists[[currn]] = preds_obj_res[[currn]]$predObj
	}
	gc()
	print ("after CV")
	# for each fold we now have the classification predictions, we merge them to one big matrix
	predictions = preds_lists[[batches[1]]]
	#print(class(predictions))
	if (is.null(dim(predictions))){
		currnames = names(predictions)
		predictions = matrix(predictions ,nrow=1)
		colnames(predictions) = currnames
	}
	for (j in 2:length(preds_lists)){
		predictions = rbind(predictions,preds_lists[[batches[j]]])
	}
	print(dim(predictions))
	print(dim(y))
	length(intersect(rownames(y),rownames(predictions)))
	predictions = predictions[rownames(y),]
	performance = getMultitaskPerformanceMeasures(predictions,y[rownames(predictions),])
	return(c(performance,list(predictions=predictions,y=y[rownames(predictions),],cv_batch_info=batch_info)))
}

applyMultitaskClassificationTest<-function(b,x,y,batch_info,classification_function=multitaskChainClassifier,prediction_args,otherRunRes=NULL,...){
	try({
		cvTestObj = getNextMultitaskTestData(b,x,y,batch_info)
		train_data=cvTestObj$train_data;tr_y=cvTestObj$tr_y
		test_data=cvTestObj$test_data
		if (is.null(otherRunRes)){
			model = classification_function(train_data,tr_y,...)
		}
		else{
			model = classification_function(train_data,tr_y,cvVec = otherRunRes[[b]]$cvVec,baseClassifiers = otherRunRes[[b]]$baseClassifiers,...)
		}
		args = c(list(model,test_data),prediction_args)
		predObj = do.call(predict,args=args)
		# make sure that the name of the sample is kept in case there is only one sample in the test data
		if (dim(test_data)[1]==1){rownames(predObj)=rownames(test_data)}
		print (paste("done with batch",b,sep = " "))
		print("#########Are the test data and prediction result are in the same order?")
		print(all(rownames(predObj)==rownames(test_data)))
		print("#####################")
		return (list(predObj=predObj))
	})
}

getNextMultitaskTestData<-function(b,x,y,batch_info){
	curr_inds = which(batch_info == b)
	tr_data = x[,-curr_inds];tr_y = y[-curr_inds,]
	te_data = x[,curr_inds]; te_y = y[curr_inds,]
	if (length(curr_inds)==1){
		te_data = matrix(te_data,ncol=1)
		colnames(te_data) = colnames(x)[curr_inds]
		rownames(te_data) = rownames(x)
		#print (te_data)
	}
	train_data= t(tr_data);test_data = t(te_data)
	obj = list(b=b,curr_inds=curr_inds,train_data=train_data,tr_y=tr_y,test_data=test_data,te_y=te_y,batch_info=batch_info)
	class(obj) = "MultitaskCVTest"
	return (obj)
}
#############################################################

# perform a validation within the dataset only if the class
# is > thr of the samples, defualt is 0.2
getMultitaskPerformancePositivesVsNegatives<-function(preds,y,currsamp2dataset){
	auprs = c()
	rocs = c()
	bgc_auprs = c()
	bgc_rocs = c()
	for (j in 1:ncol(y)){
		currname = colnames(y)[j]
		rocs[currname] = getPosVsNegsSeparation(preds[,j],y[,j],currsamp2dataset,roc=T,useBGCs=F)
		auprs[currname] = getPosVsNegsSeparation(preds[,j],y[,j],currsamp2dataset,roc=F,useBGCs=F)
		bgc_rocs[currname] = getPosVsNegsSeparation(preds[,j],y[,j],currsamp2dataset,roc=T,useBGCs=T)
		bgc_auprs[currname] = getPosVsNegsSeparation(preds[,j],y[,j],currsamp2dataset,roc=F,useBGCs=T)

	}
	return (list(rocs=rocs,auprs=auprs,bgc_rocs=bgc_rocs,bgc_auprs=bgc_auprs))
}

getPosVsNegsSeparation<-function(preds,y,d,roc=T,useBGCs=F){
	curr_datasets = unique(d[which(y==1)])
	in_datasets_inds = sapply(d,is.element,set=curr_datasets)
	inds = sample(which((in_datasets_inds & y==0) | y==1))
	if (!useBGCs){return (calcAupr(preds[inds],y[inds],roc=roc))}
	inds = sample(which((!in_datasets_inds & y==0) | y==1))
	return (return (calcAupr(preds[inds],y[inds],roc=roc)))
}

getPrecitionAndRecall <- function(thr=0.5,preds,y){
	TP = length(which(preds>=thr & y=="1"))
	FP = length(which(preds>=thr & y!="1"))
	FN = length(which(preds<thr & y=="1"))
	precision=1
	if((TP+FP)>0){precision = TP/(TP+FP)}
	recall=0
	if((TP+FN)>0){recall = TP/(TP+FN)}
	F_measure=0
	if ((precision+recall)>0){F_measure = 2*((precision*recall)/(precision+recall))}
	results = c(precision=precision, recall=recall, F=F_measure,TP=TP,FP=FP)
	return (results)
}

getRecallAtPrecision <- function(preds,y,thr=0.1){
	breaks = hist(preds,breaks=100,plot=F)[1]
	thresholds = unlist(breaks)
	res = sapply(thresholds,getPrecitionAndRecall,preds=preds,y=y)
	par(mfrow=c(1,3))
	plot(res['precision',])
	plot(res['recall',])
	plot(res['precision',],res['recall',])
	prerec = t(res[c("recall", "precision"),])[length(thresholds):1,]
	aupr=-1
	f <- approxfun(prerec)
	aupr = integrate(f, 0, 1,stop.on.error=F)$value
	try({aupr = integrate(f, 0, 1)$value})
	thr_ind = which(res[1,]<thr)
	maxRecallAthThr = 0
	if (length(thr_ind)>0){maxRecallAthThr = max(res['recall',thr_ind])}
	return (c(maxRecallAthThr = maxRecallAthThr ,thr=thr, aupr=aupr))
}

calcAupr <- function(pred, gs,roc=F,useAbs=T) {
	ord.idx = NULL
	if (useAbs){
		ord.idx <- order(abs(pred), decreasing = T)
	}
	else{
		ord.idx <- order(pred, decreasing = T)
	}
	prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
	rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
	fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate
	prec <- c(prec[1], prec)
	rec <- c(0, rec)
	fpr <- c(0, fpr)
	if (!roc){
		aupr <- areaUnderCurve(rec, prec)
		return (aupr)
	}
	#plot(fpr,rec,type='l')
	auroc <- areaUnderCurve(fpr, rec)
	return(auroc)
}

areaUnderCurve <- function(x, y) {
	dx <- diff(x)
	my <- y[1:(length(y) - 1)] + diff(y) / 2
	return(sum(dx * my))
}

getPrecisionAndRecall_multitask <- function(preds,y,thr=0.1,func = getRecallAtPrecision){
	results = c()
	for (i in 1:dim(preds)[2]){
		results = rbind(results,func(preds=preds[,i],y=y[,i],thr=thr))
	}
	rownames(results) = colnames(preds)
	return (results)
} 

getOR<-function(x,y,q=0.6666){
	thr = quantile(x,q)
	above_thr = which(x>thr)
	p_above = sum(y[above_thr])/length(above_thr)
	p_below = sum(y[-above_thr])/length(y[-above_thr])
	or = p_above/p_below
	return (or)
}

getROCPval<-function(x,y){
	if (length(unique(y))<2){return (1)}
	x0 = x[which(y==0)]
	x1 = x[which(y==1)]
	# test if x1>x0
	wtest = wilcox.test(y=x0,x=x1,alternative='greater')
	return (wtest$p.value)
}

getORTwoThrs<-function(x,y,p=0.333,q=0.6666){
	thr1 = quantile(x,q)
	above_thr = which(x>=thr1)
	thr2 = quantile(x,p)
	below_thr = which(x<=thr2)
	p_above = sum(y[above_thr])/length(above_thr)
	p_below = sum(y[-above_thr])/length(y[-above_thr])
	or = p_above/p_below
	return (or)
}

# calculate cases in which there is a contradiction
# in prediction between parents and descendants
calcNumContradictions<-function(preds,y,doancs,thr=0.5){
	n = ncol(y);m=nrow(y)
	dnames = colnames(y)
	numConts1 = matrix(ncol=n,nrow=m)
	numConts1[,] = F
	numConts2 = matrix(ncol=n,nrow=m)
	numConts2[,] = F
	total_preds = 0
	for (i in 1:n){
		d1 = dnames[i]
		for (j in 1:n){
			d2 = dnames[j]
			if (i==j || !isAncestor(d1,d2,doancs)){next}
			#print (paste(d1,d2))
			numConts1[,i] = numConts1[,i] | analyzeContradictionsOfDiseasePair(preds[,i],preds[,j],thr)
			numConts2[,i] = numConts2[,i] | analyzeContradictionsVsRealLabels(preds[,i],preds[,j],y[,i],thr)
			total_preds = total_preds + nrow(y)
		}
	}
	return (c(theoreticalErrors = sum(numConts1),realErrors = sum(numConts2),totalChecked =	n*m ))
}

# here we know that 1 is an ancestor of 2
analyzeContradictionsOfDiseasePair<-function(preds1,preds2,thr){
	return (preds2>=thr & preds1<thr)
}

analyzeContradictionsVsRealLabels<-function(preds1,preds2,y1,thr){
	return (y1==1 & preds2>=thr & preds1<thr)
}

analyzeContradictionsOfASample<-function(preds,dnames,doancs,thr){
	numConts = 0
	n = length(preds)
	for (i in 1:n){
		for (j in 1:n){
			if (i==j){next}
			numConts = numConts + isContradiction(preds[i],preds[j],thr,dnames[i],dnames[j],doancs)
		}
	}
	return (numConts)
}

# is 1 ancestor of 2 and pred2=1 and pred1=0?
isContradiction<-function(pred1,pred2,thr,d1,d2,doancs){
	if (is.null(doancs)){doancs = as.list(DOANCESTOR)}
	if (!isAncestor(d1,d2,doancs)){return (FALSE)}
	if (pred2>=thr && pred1<thr){return (TRUE)}
	return (FALSE)
	
}

isAncestor<-function(d1,d2,doancs){
	return (is.element(d1,set=doancs[[d2]]))
}


#### New methods for our selection methods #####
Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(p.val)
}

getStouferPval<-function(zs){
	Z  <- sum(zs)/sqrt(length(zs))
	p.val <- 1-pnorm(abs(Z))
	return (p.val)
}
# Return the joint q-values for a set of datasets
get_datasets_qvals<-function(currd1,dataset_disease_zscores,disease,standardize = F){
	currd1_mat = c()
	for (dataset in currd1){
		if (is.element(disease,set = names(dataset_disease_zscores[[dataset]]))){
			currd1_mat = cbind(currd1_mat,dataset_disease_zscores[[dataset]][[disease]])
		}
	}
	#if (length(currd1_mat) == 0){return (NULL)}
	print (dim(currd1_mat))
	if (standardize){currd1_mat = apply(currd1_mat,2,function(x){(x-mean(x))/sd(x)})}
	curr_pvals1 = apply(currd1_mat,1,getStouferPval)
	return (p.adjust(curr_pvals1,method='fdr'))
}
# get the roc scores of a set of datasets for a specific disease
get_roc_scores<-function(currd1,curr_others1,x,y,disease,d){
	curr_samples = names(d)[is.element(d,set=union(currd1,curr_others1))]
	curr_samples = intersect(curr_samples,colnames(x))
	curr_x = x[,curr_samples]
	curr_y = y[curr_samples,disease]
	curr_d = d[curr_samples]
	curr_datasets = unique(d[which(curr_y==1)])
	curr_in_datasets_inds = is.element(curr_d,set=curr_datasets)
	curr_negs = which(curr_in_datasets_inds & curr_y==0)
	curr_positives  = which(curr_y==1)
	curr_bgcs = which(!curr_in_datasets_inds & curr_y==0)
	feature_pos_neg_rocscores = apply(curr_x[,c(curr_positives,curr_negs)],1,calcAupr,gs=curr_y[c(curr_positives,curr_negs)],roc=T)
	feature_pos_neg_rocscores = pmax(feature_pos_neg_rocscores,1-feature_pos_neg_rocscores)
	feature_pos_bgc_rocscores = apply(curr_x[,c(curr_positives,curr_bgcs)],1,calcAupr,gs=curr_y[c(curr_positives,curr_bgcs)],roc=T)
	feature_pos_bgc_rocscores = pmax(feature_pos_bgc_rocscores,1-feature_pos_bgc_rocscores)
	return (list(feature_pos_neg_rocscores,feature_pos_bgc_rocscores))
}

selectedFeaturesClassifier<-function(x,y,classification_function=randomForest,fset,...){
	classifier = classification_function(x[,fset],y,...)
	obj = list(inds=fset,classifier=classifier)
	class(obj) = "featureSelectionClassifier"
	return (obj)
}

predict.selectedFeaturesClassifier<-function(obj,newx,...){
	args = c(list(obj$classifier,newx[,obj$inds]),...)
	return (do.call(predict,args=args))
}

binary_relevance_with_ourFS<-function(x,y,d,dataset_disease_zscores,numFeaturesPerType = 100,
		baseClassifier=liblinear_wrapper,baseClassifierArgs=list(),baseClassifierPredArgs=list()){
	disease2features = list()
	disease2qvals = list()
	disease2rocs = list()
	disease2model = list()
	d = d[rownames(y)]
	for (disease in colnames(y)){
		curr_y = y[,disease]
		curr_datasets = unique(d[which(curr_y==1)])
		curr_other_datasets = setdiff(unique(d),curr_datasets)
		in_datasets_inds = is.element(d,set=curr_datasets)
		negs = which(in_datasets_inds & curr_y==0)
		positives  = which(curr_y==1)
		bgcs = which(!in_datasets_inds & curr_y==0)
		qvals = get_datasets_qvals(curr_datasets,dataset_disease_zscores,disease,standardize = F)
		qvals = qvals[colnames(x)]
		#print (quantile(qvals))
		rocscores = get_roc_scores(curr_datasets,curr_other_datasets,t(x),y,disease,d)
		disease2qvals[[disease]] = qvals
		disease2rocs[[disease]] = rocscores
		# select features
		set1 = names(sort(rocscores[[1]],decreasing=T))[1:numFeaturesPerType]
		if (length (curr_other_datasets) > 0){
			set2 = names(sort(rocscores[[2]],decreasing=T))[1:numFeaturesPerType]
		}
		else{ set2 = c()}
		set3 = names(sort(qvals,decreasing=F))[1:numFeaturesPerType]
		disease2features[[disease]] = union(set1,set2)
		disease2features[[disease]] = union(disease2features[[disease]],set3)
		#print (disease2features[[disease]])
		print (paste(term2name(disease),length(disease2features[[disease]]),sep="     "))
		args = c(list(x=x[,disease2features[[disease]]],y=as.factor(as.character(curr_y))),baseClassifierArgs)
		disease2model[[disease]] = do.call(baseClassifier,args=args)
	}
	obj = list(disease2model = disease2model ,disease2features = disease2features,baseClassifierPredArgs=baseClassifierPredArgs )
	class(obj) <- "binary_relevance_with_ourFS"
	return (obj)
}

predict.binary_relevance_with_ourFS<-function(obj,newx,...){
	if (is.null(dim(newx))){newx = matrix(newx,nrow=1)}
	numlabels = length(obj$disease2model)
	numsamples = dim(newx)[1]
	numclasses = length(obj$disease2model)
	preds = matrix(rep(0,times=numclasses*numsamples),nrow=numsamples,ncol=numclasses)
	rownames(preds) = rownames(newx)
	colnames(preds) = names(obj$disease2model)
	for (disease in names(obj$disease2model)){
		args = c(list(obj$disease2model[[disease]],newx[,obj$disease2features[[disease]]]),obj$baseClassifierPredArgs,...)
		predsMat = do.call(predict,args=args)
		predsMat = getPredProbabilities(predsMat)
		currClassInd = which(colnames(predsMat)=="1")
		preds[,disease] = predsMat[,currClassInd]
	}
	return (preds)
}

remove_lables_for_classification<-function(y,min_samps=10,rm_doancs=T){
	labels2remove = names(which(apply(y,2,sd) < 1e-5 | colSums(y)<min_samps))
	if (rm_doancs){
		for(i in 1:ncol(y)){
			do1 = colnames(y)[i]
			for (j in 1:ncol(y)){
				do2 = colnames(y)[j]
				if (i==j){next}
				if (!is.element(do1,set=do_ancs[[do2]])){next}
				if (all(y[,do1]==y[,do2])){
					labels2remove = union(labels2remove,do2)
				}
			}
		}
	}
	labels2keep = setdiff(colnames(y),labels2remove)
	y = y[,labels2keep]
	return (y)
}

remove_lables_for_classification_use_datasets<-function(y,d,min_d=3){
	labels2remove = c()
	for(i in 1:ncol(y)){
		samps = rownames(y)[y[,i]==1]
		numd = length(unique(d[samps]))
		if (numd < min_d){labels2remove = c(labels2remove,i)}
	}
	y = y[,-labels2remove]
	return (y)
}

remove_features_for_gene_analysis<-function(x,minSum=10){
	# FILTER 1
	fsums = colSums(x)
	to_remove = which(fsums<10)
	x = x[,-to_remove]
	# FILTER 2: keep the genes
	x_genes = colnames(x)[!grepl(colnames(x),pattern=';')]
	x = x[,x_genes]
	return (x)
}

# Transform standard label matrix y to contain negatives
getLabelMatrixWithNegatives<-function(y,d,neg_val = -1){
	if (!is.null(rownames(y))){d = d[rownames(y)]}
	newy  = apply(y,2,getVectorWithLabels,d=d,neg_val=neg_val)
	rownames(newy) = rownames(y)
	colnames(newy) = colnames(y)
	return (newy)
}

# here index i: y[i] is the label of sample i and
# d[i] is its dataset
getVectorWithLabels<-function(y,d,neg_val = -1){
	pos_inds = which(y==1)
	#print (pos_inds )
	pos_inds_datasets = unique(d[pos_inds])
	#print (pos_inds_datasets)
	is_in_dataset = which(is.element(d,set = pos_inds_datasets))
	neg_inds = setdiff(is_in_dataset,pos_inds)
	newy = rep(0,length(y))
	newy[pos_inds] = 1
	newy[neg_inds] = neg_val
	return(newy)
}

# LP (label powerset) multilabel classifier
# This is similar to multilabel via multiclass
# ... are arguments for multitaskViaMulticlassClassifier
LPClassifierWithNegatives<-function(x,y,...){
	classifier = multitaskViaMulticlassClassifier(x,y,...)
	class(classifier) = "LPClassifierWithNegatives"
	return (classifier)
}

predict.LPClassifierWithNegatives<-function(obj,newx,...){
	if (is.null(dim(newx))){newx = matrix(newx,nrow=1)}
	args = c(list(obj$model,newx[,obj$features]),obj$baseClassifierPredArgs,...)
	# predObj is a matrix: rows are samples 
      predObj = do.call(predict,args=args)
	preds = t(apply(predObj,1,LPClassifierWithNegatives_transformPrediction,xnames=colnames(predObj)))
	print (dim(preds))
	# fix in case we predict a single sample
	if (is.null(dim(preds))){preds = matrix(preds,nrow=1)}
	rownames(preds) = rownames(newx)
	colnames(preds) = obj$original_tasks
	return (preds)
}

LPClassifierWithNegatives_transformPrediction<-function(x,xnames){
	names_matrix = sapply(xnames,function(x)strsplit(x,split=',')[[1]])
	mode(names_matrix) = 'numeric'
	names_matrix[names_matrix==-1] = 0
	return(names_matrix%*%x)
}

# Wrappers for our Java alagorithms 
# trx, try: the training set
# tex: the test set on which we want to predict
# fs = do feature selection? if fs<=0 then use all features
PCTsRFClassifier<-function(x,y,labelweights = NULL, jarpath = "GetPCTRFModel.jar", xmx = "-Xmx8g",ntree=500,mtry=NULL,nodesize=1,bagsize=100,alpha=0.5,fs=1000){

	features = colnames(x)
	if (fs >0){
		feature_scores = apply(x,2,sum)
		features = names(sort(feature_scores,decreasing=T))[1:fs]
		x = x[,features]
		print (dim(x))
		gc()
	}

	print ("features selected")

	num_models_in_wd = sum(grepl(pattern="(rfmodel.*mulan)|(rfmodel.*csv)",dir(getwd()),perl=T))
	rand_suffix = runif(1,min=1,max=1000)
	outfile = paste("rfmodel_",num_models_in_wd+1,"_",rand_suffix,".mulan",sep="")
	csvoutfile = paste("rfmodel_training_",num_models_in_wd+1,"_",rand_suffix,".csv",sep="")

	print (paste("printing: " ,csvoutfile))

	# concatanate x and y and print to a csv file
	feature_names = colnames(x)
	label_names = colnames(y)
	xnames = paste("x_",1:length(feature_names),sep="")
	ynames = paste("y_",1:length(label_names),sep="")
	if (!is.null(labelweights)){
		ynames = paste (ynames,labelweights,sep=";")
	}
	trdata = cbind(x,y)
	colnames(trdata) = c(xnames,ynames)
	write.csv(trdata,file = csvoutfile, row.names=F,quote=F)

	if (is.null(mtry)){
		mtry = floor(sqrt(ncol(x)))
	}
	
	system(paste("export _JAVA_OPTIONS=",xmx,sep=''))
	
	cmd = paste("java -jar -XX:-UseGCOverheadLimit ",xmx," ",jarpath," ",csvoutfile," ",outfile," mtry=",mtry," alpha=",alpha," ntree=",ntree," nodesize=",nodesize," bagpercent=",bagsize,sep="")
	print (cmd)
	system(cmd)
	print (outfile)
	file.remove(csvoutfile)

	obj = list(feature_names=feature_names,label_names=label_names,xnames=xnames,ynames=ynames,
		labelweights=labelweights,wd=getwd(),modelfile=outfile,mtry=mtry,ntree=ntree,nodesize=nodesize,alpha=alpha)
	class(obj) = "PCTsRFClassifier"
	return (obj)
}
predict.PCTsRFClassifier<-function(obj,newdata,jarpath="TestPCTRFModel.jar",xmx = "-Xmx5000m",deleteModelFile=T){

	if (is.null(dim(newdata))){newdata=matrix(newdata,nrow=1)}
	newdata = newdata[,obj$feature_names]

	num_tests_in_wd = sum(grepl(pattern="_testset.csv$",dir(getwd()),perl=T))
	rand_suffix = runif(1,min=1,max=1000)
	outfile = paste("rfmodel_run_",num_tests_in_wd+1,"_",rand_suffix,"_testset_preds.csv",sep="")
	csvoutfile = paste("rfmodel_run_",num_tests_in_wd+1,"_",rand_suffix,"_testset.csv",sep="")

	# print to csv
	feature_names = colnames(newdata)
	xnames = paste("x_",1:length(feature_names),sep="")
	colnames(newdata) = xnames
	write.csv(newdata,file = csvoutfile, row.names=F,quote=F)
	
	cmd = paste("java -jar ",xmx,jarpath,csvoutfile,obj$modelfile,outfile)
	system(cmd)
	file.remove(csvoutfile)

	preds = as.matrix(read.csv(outfile,header=F))
	mode(preds) = 'numeric'
	colnames(preds) = obj$label_names
	file.remove(outfile)
	rownames(preds) = rownames(newdata)

	if (deleteModelFile){
		file.remove(obj$modelfile)
	}
	return (preds)
}


# Wrappers for our Java runnables based on Mulan 
# fs = do feature selection? if fs<=0 then use all features
# additionals = additional string that specifies the command line 
# parameters
MulanClassifier<-function(x,y,labelweights = NULL, jarpath = "GetMLkNNModel.jar", xmx = "-Xmx8g",fs=1000,file_rm = T,additionals = "k=20"){
	features = colnames(x)
	if (fs >0){
		feature_scores = apply(x,2,sum)
		features = names(sort(feature_scores,decreasing=T))[1:fs]
		x = x[,features]
		print (dim(x))
		gc()
	}
	print ("features selected")

	num_models_in_wd = sum(grepl(pattern="(rfmodel.*mulan)|(rfmodel.*csv)",dir(getwd()),perl=T))
	rand_suffix = runif(1,min=1,max=1000)
	outfile = paste("rfmodel_",num_models_in_wd+1,"_",rand_suffix,".mulan",sep="")
	csvoutfile = paste("rfmodel_training_",num_models_in_wd+1,"_",rand_suffix,".csv",sep="")
	print (paste("printing: " ,csvoutfile))

	# concatanate x and y and print to a csv file
	feature_names = colnames(x)
	label_names = colnames(y)
	xnames = paste("x_",1:length(feature_names),sep="")
	ynames = paste("y_",1:length(label_names),sep="")
	if (!is.null(labelweights)){
		ynames = paste (ynames,labelweights,sep=";")
	}
	trdata = cbind(x,y)
	colnames(trdata) = c(xnames,ynames)
	write.csv(trdata,file = csvoutfile, row.names=F,quote=F)

	system(paste("export _JAVA_OPTIONS=",xmx,sep=''))
	
	cmd = paste("java -jar -XX:-UseGCOverheadLimit ",xmx," ",jarpath," ",csvoutfile," ",outfile," ",additionals,sep="")
	print (cmd);system(cmd)
	if (file_rm){file.remove(csvoutfile)}

	obj = list(feature_names=feature_names,label_names=label_names,xnames=xnames,ynames=ynames,
		labelweights=labelweights,wd=getwd(),modelfile=outfile,additionals=additionals)
	class(obj) = "MulanClassifier"
	return (obj)
}
predict.MulanClassifier<-function(obj,newdata,jarpath="GetModelPredictionsOnTestSet.jar",xmx = "-Xmx5000m",deleteModelFile=T,deleteCSV=F){

	if (is.null(dim(newdata))){newdata=matrix(newdata,nrow=1)}
	newdata = newdata[,obj$feature_names]

	ynames = paste("y_",1:length(obj$label_names),sep="")
	if (!is.null(obj$labelweights)){
		ynames = paste (ynames,obj$labelweights,sep=";")
	}
	dummy_y_mat = matrix("?",nrow=nrow(newdata),ncol=length(ynames))
	colnames(dummy_y_mat) = ynames
	newdata = cbind(newdata,dummy_y_mat)

	num_tests_in_wd = sum(grepl(pattern="_testset.csv$",dir(getwd()),perl=T))
	rand_suffix = runif(1,min=1,max=1000)
	outfile = paste("rfmodel_run_",num_tests_in_wd+1,"_",rand_suffix,"_testset_preds.csv",sep="")
	csvoutfile = paste("rfmodel_run_",num_tests_in_wd+1,"_",rand_suffix,"_testset.csv",sep="")

	# print to csv
	feature_names = colnames(newdata)
	xnames = paste("x_",1:length(feature_names),sep="")
	colnames(newdata) = xnames
	write.csv(newdata,file = csvoutfile, row.names=F,quote=F)
	
	cmd = paste("java -jar -XX:-UseGCOverheadLimit ",xmx,jarpath,csvoutfile,obj$modelfile,outfile)
	system(cmd)
	if (deleteCSV){file.remove(csvoutfile)}

	preds = as.matrix(read.csv(outfile,header=F))
	mode(preds) = 'numeric'
	colnames(preds) = obj$label_names
	file.remove(outfile)
	rownames(preds) = rownames(newdata)

	if (deleteModelFile){
		file.remove(obj$modelfile)
	}
	return (preds)
}


chisq_test_numeric_vs_nominal <-function(y,x,bins=20){
	disc_x = cut(x,breaks=bins,labels=F)
	tab = table(disc_x,y)
	suppressWarnings({p = (chisq.test(tab)$p.value)})
	return (p)
}

chisq_analyze_one_feature_vs_matrix<-function(x,y){
	pvals = apply(y,2,chisq_test_numeric_vs_nominal,x=x)
	return (min(pvals))	
}

chisq_get_feature_scores<-function(x,y){
	scores = apply(x,2,chisq_analyze_one_feature_vs_matrix,y=y)
	return (scores)
}

# Tests
#y = matrix(ncol=3,nrow=120)
#y[,] = 0
#d = rep("A",nrow(y))
#for (j in 1:ncol(y)){
#	currn = nrow(y)/ncol(y)
#	curr_samps = ((j-1)*currn+1):(j*currn)
#	y[curr_samps,j] = 1
#	d[curr_samps] = paste("A")
#}
#y1 = getLabelMatrixWithNegatives(y,d)
#all(y1==y)

# RF using the ranger package
try({library('ranger',lib.loc='rlibs/')})
try({library(ranger)})
# printErr - binary, print the error if ranger fails
ranger_wrapper<-function(x,y,printErr=T,...){
	o = list()
	succ = F
	tryobj = try({
		o[[1]] = ranger(y~.,data=data.frame(x,y,check.names=F),write.forest=T,classification=T,
			probability=T,save.memory=T,min.node.size=3,importance='impurity',...)
		succ=T
		# correct the forest variable names
		#newnames = o[[1]]$forest$independent.variable.names
		#newnames = substring(newnames,3)
		#o[[1]]$forest$independent.variable.names = newnames
	})
	if (!succ){
		if(printErr){print(tryobj)}
		o[[1]] = randomForest(x,y)
	}
	print(class(o[[1]]))
	class(o) = 'ranger_wrapper'
	return (o)
}
predict.ranger_wrapper<-function(o,d){
	if (class(o[[1]]) == "randomForest"){return(predict(o[[1]],d,type='prob'))}
	return(predict(o[[1]],d,verbose=F)$predictions)
}

# p - positives, n = negatives, then:
simple_sampling_based_learner<-function(x,y,d=patient2dataset,func = featureSelectionClassifier,reps=10,max_num_samples = 1000,...){
	positives = rownames(x)[y=="1"]
	datasets = unique(d[positives])
	negatives = rownames(x)[y=="0" & is.element(d[rownames(x)],set=datasets)]
	bgcs = setdiff(rownames(x),positives);bgcs = setdiff(bgcs,negatives)
	classifiers = list()
	samp_size = min(c(length(positives),max_num_samples))
	for (i in 1:reps){
		newbgcs = NULL
		if(length(bgcs)>0){
			newbgcs = sample(bgcs)[1:min(samp_size,length(bgcs))]
		}
		pos_sample = sample(positives)[1:samp_size]
		neg_sample = NULL
		if(length(negatives)>0){
			neg_sample = sample(negatives)[1:min(samp_size,length(negatives))]
		}
		print(paste("number of positives, negatives, and bgcs:"))
		print(c(length(pos_sample),length(neg_sample),length(newbgcs)))
		curr_samples = c(pos_sample,neg_sample,newbgcs)
		curr_inds = is.element(rownames(x),curr_samples)
		newy = y[curr_inds]
		newx = x[curr_inds,]
		classifiers[[i]] = func(newx,as.factor(newy),...)
	}
	obj = classifiers
	class(obj) = "simple_sampling_based_learner"
	return(obj)
}
predict.simple_sampling_based_learner<-function(obj,x,...){
	preds = c()
	counts = 0
	for (i in 1:length(obj)){
		if (length(preds)==0){
			try({
				preds = predict(obj[[i]],x);
				counts = counts + 1
			})
		}
		else{
			try({
				preds = preds + predict(obj[[i]],x);
				counts = counts + 1
			})
		}
	}
	preds = preds / counts
	return (preds)
}

# mocd = 1:nrow(xTrain);names(mocd)=rownames(xTrain)
# o = simple_sampling_based_learner(xTrain,yTrain,d=mocd,func=liblinear_binary_svm_wrapper)
# predict(o,xTrain[1:10,])

SMOTE_sampling_based_learner<-function(x,y,d=patient2dataset,func = featureSelectionClassifier,
		reps=1,numFeatures=100,ntree=500,perBGCRemoval=0.5,binarize=T,...){
	positives = rownames(x)[y=="1"]
	datasets = unique(d[positives])
	negatives = rownames(x)[y=="0" & is.element(d[rownames(x)],set=datasets)]
	n = length(positives)+length(negatives)
	if (n >= nrow(x)/2){
		rfs = list()
		rfs[[1]] = func(x,as.factor(y),classification_function = ranger_wrapper,
			numFeatures=numFeatures,num.tree=ntree*reps,...)
		obj = list(rf = rf,importance=importance)
		class(obj) = "simple_sampling_based_learner"
		return(obj)
	}
	bgcs = setdiff(rownames(x),positives)
	bgcs = setdiff(bgcs,negatives)
	rfs = list()
	importance = rep(0,ncol(x));names(importance) = colnames(x)
	for (i in 1:reps){
		curry = as.factor(y)
		#bgcs to exclude, was added to reduce SMOTE's running time
		bgcs_to_rem = sample(bgcs)[1:(perBGCRemoval*length(bgcs))]
		to_rem = is.element(rownames(x),set=bgcs_to_rem)
		newData = SMOTE(curry~.,data=data.frame(x=x[!to_rem,],curry=curry[!to_rem]))
		newData = as.matrix(newData)
		newy = newData[,ncol(newData)]
		print("newy after SMOTE:");print(table(newy))
		newx = newData[,-ncol(newData)];colnames(newx) = colnames(x)
		mode(newx) = 'numeric'
		# SMOTE does not return a binary dataset: we binarize it before moving on
		if(binarize){
			newx = newx>0.5
			mode(newx) = 'numeric'
		}
		rfs[[i]] = func(newx,as.factor(newy),classification_function = ranger_wrapper,
			numFeatures=numFeatures,num.tree=ntree,...)
		imps = rfs[[i]]$classifier[[1]]$variable.importance
		importance[names(imps)] = importance[names(imps)]+imps
	}
	obj = list(rfs = rfs,importance=importance)
	class(obj) = "simple_sampling_based_learner"
	return(obj)
}
# Test on the cosmic data with the second label: DOID:3119 (>1000 samples):
# > system.time(SMOTE(y~.,data=data.frame(x=x[1:1000,],y=y[1:1000])))
#   user  system elapsed
# 113.024   0.992 114.006
# > system.time(SMOTE(y~.,data=data.frame(x=x,y=y)))
#     user    system   elapsed
# 10146.596   687.664 10833.485

# fast sampling-based algorithm using randomForest
# too slow...
# ranger is better anyway
simple_sampling_based_learner_randomForest<-function(x,y,d=patient2dataset,reps=10,classifier = featureSelectionClassifier,numFeatures=1000,ntree=100,...){
	positives = rownames(x)[y=="1"]
	datasets = unique(d[positives])
	negatives = rownames(x)[y=="0" & is.element(d[rownames(x)],set=datasets)]
	n = length(positives)+length(negatives)
	print (paste("num_pos:",length(positives),", num positives + num negatives: ",n))
	if (n >= nrow(x)/2){
		rf = classifier(x,as.factor(y),ntree=reps*ntree,numFeatures=numFeatures,...)$classifier
		importance = rf$importance[,1]
		obj = list(rf = rf,importance=importance)
		class(obj) = "simple_sampling_based_learner"
		return(obj)
	}
	bgcs = setdiff(rownames(x),positives)
	bgcs = setdiff(bgcs,negatives)
	rf = NULL
	importance = rep(0,ncol(x));names(importance) = colnames(x)
	for (i in 1:reps){
		newbgcs = sample(bgcs)[1:n]
		newy = y[c(positives,negatives,newbgcs)]
		newx = x[c(positives,negatives,newbgcs),]
		curr_rf = classifier(newx,as.factor(newy),ntree=ntree,numFeatures=numFeatures,...)$classifier
		if(is.null(rf)){
			rf = curr_rf
		}
		else{
			curr_imp = curr_rf$importance
			importance[rownames(curr_imp)] = importance[rownames(curr_imp)] + curr_imp[,1]
			curr_rf$importance = rf$importance
			rf = randomForest::combine(rf,curr_rf)
		}
	}
	obj = list(rf = rf,importance=importance)
	class(obj) = "simple_sampling_based_learner"
	return(obj)
}
predict.simple_sampling_based_learner_randomForest<-function(obj,x,...){predict(obj$rf,x,type='prob')}

#newrf = simple_sampling_based_learner(x,y,d,reps=3)
#predict(newrf,x)
#simple_binary_relevance(x2,y2[,1:5],baseClassifier = simple_sampling_based_learner, baseClassifierArgs=list(d=patient2dataset,reps=2,numFeatures=10))
#rfs = simple_sampling_based_learner(x2,as.factor(y2[,1]),d=patient2dataset,reps=5,ntree=5,numFeatures=100)
#preds = predict(rfs,x2)
#plot(preds[,1],y2[,1])

getMultitaskPerformanceMeasures<-function(preds,y,useROCR=T,asMat=F,...){
	class2ROC = c()
	class2AUPR = c()
	class2Fmeasure = c()
	class2Accuracy = c()
	class2Size = c()
	m = dim(preds)[2]
	for (i in 1:m){
		curr_task = colnames(y)[i]
		curr_labels = as.factor(y[,i])
		curr_preds = preds[,i]
		to_r = which(is.na(curr_preds))
		if (length(to_r)>0){
			curr_labels = curr_labels[-to_r]
			curr_preds = curr_preds[-to_r]
		}
		if (length(unique(curr_labels))<2){next}
		class2ROC[curr_task] = getROCScore(curr_preds,curr_labels,useROCR=useROCR)
		class2Size[curr_task] = sum(y[,i],na.rm=T)
		class2Accuracy[curr_task] = length(which((curr_preds>=0.5 & curr_labels=="1")|(curr_preds<0.5 & curr_labels=="0")))/length(curr_labels)
		curr_precision_and_recall = getPrecitionAndRecall(preds=curr_preds,y=curr_labels)
		class2Fmeasure[curr_task] = curr_precision_and_recall['F'] 
		class2AUPR[curr_task] = calcAupr(curr_preds,as.numeric(as.character(curr_labels)))
	}
	if (asMat){
		summmat = rbind(class2AUPR,class2ROC,class2Accuracy,class2Size)
		return (t(summmat))
	}
	return (list(class2ROC =class2ROC, class2Accuracy=class2Accuracy,class2Size=class2Size,class2Fmeasure=class2Fmeasure,class2AUPR=class2AUPR))	
}

