#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, glmnet, randomForest, Hmisc"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-x", "--file_x"), type="character", default="-",
     help="Matrix for feature (row for cell, col for feature) [%default] \n\t\tUse '-' if passing 'x' with a UNIX pipe.", metavar="character"),
    make_option(c("-y", "--file_y"), type="character", default="-",
     help="Matrix for drug response (row for drug, col for cell) [%default] \n\t\tUse '-' if passing 'y' with a UNIX pipe.", metavar="character"),
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    make_option(c("-d", "--drug"), type="character", default="1",
     help="A string deciding which drug will be used [%default] \n\t\tBoth row number and row name are acceptable.", metavar="character"),
    make_option(c("-m", "--model"), type="character", default="EN",
     help="A string deciding which statistical model will be used [%default] \n\t\tOptions include 'EN' (Elastic Net) and 'RF' (Random Forests).", metavar="character"),
    make_option(c("--nsteps"), type="integer", default="1000",
     help="Number of models for building [%default]", metavar="integer"),
    make_option(c("--nfolds"), type="integer", default="10",
     help="Number of folds for each model [%default]", metavar="integer"),
    make_option(c("--balancing"), type="logical", default="FALSE", action = "store_true",
     help="Only if this value is TRUE the trainSet will be balanced? [%default]", metavar="logical"),

    # fread
    make_option(c("--col_sep"), type="character", default="\t",
     help="Separator between columns [\\t] \n\t\t'auto' represent for '[,\\t |;:]'. See also 'data.table::fread'.", metavar="character"),
    make_option(c("--row_num"), type="integer", default="-1",
     help="Number of rows to read [%default] \n\t\t'-1' means all. '0' is a special case that just returns the column names.", metavar="integer"),
     make_option(c("--col_name"), type="character", default="auto",
     help="The first data line contain column names, 'auto'|TRUE|FALSE [%default] \n\t\tDefaults according to whether every non-empty field on the first data\n\t\tline is type character.", metavar="character"),
    make_option(c("--na_str"), type="character", default="NA",
     help="String(s) which are to be interpreted as NA values [%default]", metavar="character"),
    make_option(c("--row_skip"), type="character", default="0",
     help="Row can be taken as the first data row [%default] \n\t\tskip='string' searches for 'string' in the file and starts on that row.", metavar="character"),
    make_option(c("--row_name"), type="character", default="1",
     help="Column number (or name) can be taken as the row names [%default]", metavar="character")
)

opt_parser = OptionParser(usage="usage: %prog [options] -x <file> -y <file> -z <file> -o <dir> -p 2 --cv --iter_y \n\tBy default <file> is -' (stdin)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

x = opt$options$file_x
y = opt$options$file_y
o = opt$options$dir_out
d = opt$options$drug
stat_model = opt$options$model
nSteps = opt$options$nsteps
nFolds = opt$options$nfolds
balancing = opt$options$balancing

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

if(grepl("^[[:digit:]]+$", d)) d=as.numeric(d)

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl('^[[:digit:]]+$', ra)) ra=as.numeric(ra)

if(x=="-") x = "file:///dev/stdin"
if(y=="-") y = "file:///dev/stdin"

# read
features = fread(x, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(features) = features[,ra]
features = as.matrix(features[,-ra])
features = features[apply(features, 1, function(x) sum(is.na(x))!=length(x)),]
features = features[,apply(features, 2, function(x) !any(is.na(x)))]

drug_res = fread(y, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(drug_res) = drug_res[,ra]
drug_res = as.matrix(drug_res[,-ra])

# main
library(randomForest)
library(Hmisc)
library(glmnet)

if(!dir.exists(o)) dir.create(o)

drug_id <- ifelse(is.numeric(d), rownames(drug_res)[d], d)

#####################
### LIB_BALANCING ###
#####################

# ------------------------------------------------------------------------------------------
# DESC: Balance a vector of response data, by splitting the data in n bins of equal size. 
#       This function resamples the data so that each bin has afterwards the size of the 
#       largest bin. Empty bins remain empty. In case the modulo of largest bin and bin to 
#       balance is larger than zero, randomly samples will be taken from the bin to balance
#       for equalling size with the largest bin. 
# IN:   res   =>  Vector containing the drug response data. This vector needs to have the 
#                 the COSMIC ids as name.
#       nBins =>  Number of splits (default = 10).
# OUT:  Returns a balanced response vector. 
# ------------------------------------------------------------------------------------------
balance <- function(res, nBins=10) {
  
  if (length(res) < 10) 
    stop("The response vector neesds at least 10 elements!")
  
  if (!is.numeric(res)) 
    stop("Response need to be a numeric vector!")
  
  if (is.null(names(res)))
    stop("Response vector needs names, which are the ids of each entry!")
  
  if (length(unique(names(res))) != length(names(res)))
    stop("The names of the response vector being not unique IDs!")
  
  membership <- as.numeric(cut(res, nBins))
  counts<-table(membership)
  
  # how often each bin will be resampled
  nResample_new <- floor(max(counts) / counts)
  
  # how many random samples taken from each bin
  nDraw_new <- max(counts) %% counts
  
  bRes <- c()
  for (member in names(counts)) {
    binRes <- res[membership == member]
    bRes <- c(bRes, rep(binRes, nResample_new[member]))
    # only if a random number of sample has to be drawn to fill up the bin
    if (nDraw_new[member]>0) {
      bRes <- c(bRes, sample(binRes, nDraw_new[member]))
    }
  }
  
  # hist is ignoring "breaks" parameter, therefore alternative visualization. 
  # plot(cut(res, nBins), las=2)
  # plot(cut(bRes, nBins), las=2)
  
  return(sample(bRes))
}

###################
### LIB_METRICS ###
###################

# obs and pred need to be numeric for all functions!

# function to calculate the concordance index
cIDX <- function(pred, obs) {
  # calculate the concordance index 
  rcorr.cens(pred, obs)[1]
}

# function to calculate the root mean square error
# observed (independent variable, x)
# predicted (dependent variable, y)
RMSE <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE)) 
}

# coefficient of determination / R squared
R2 <- function(obs, pred) {
  # total variation in the observed data
  SSyy = sum((obs-mean(obs))^2)
  # sqared error of prediction
  SSE = sum((obs-pred)^2)
  1 - SSE / SSyy
}

# The proportion of variance explained (PVE)
PVE <- function(obs, pred) {
  # covariance sum of squares
  SSxy <- sum((obs-mean(obs)) * (pred-mean(pred)))
  # sum of squares for variable X
  SSx  <- sum((obs-mean(obs))^2)
  # sum of squares for variable Y
  SSy  <- sum((pred-mean(pred))^2)
  
  numerator <- SSxy^2
  denominator <- SSx*SSy
  if (denominator > 0)
    return(numerator/ denominator)
  else
    return(0)
}

# #######################################################################################################################################
# RandomForest Function
# #######################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------
#DESC:  The function "trainRF" trains a "randomForest" model with the data from "train_feat". 
#       This function will return a list containing the best model and a sublist with the c-indxes for the train-, xTrain-set.
#
#IN:    train_feat      ==> train values obtained from crossvalidating the observated feature data.
#                           the model will be trained with this values
#       train_obs       ==> observed drug_response for the values of train_feat
#       xTrain_feat     ==> xTrain values obtained from crossvalidating the observed feature data.
#                           necessary to figure out when to stop the training of the model with train_feat
#       xTrain_obs      ==> observed drug_response for the values of xTrain_feat
#       nTreeOffset     ==> the number of trees the first generated model will have
#       nTreeStep       ==> number of trees to grow(or to be added to an existing ensemble of trees)
#                           default value = 1
#       maxTrees        ==> maximum value for the value of nTrees
#                           if nTrees bigger than maxTrees the training will be stoppped
#                           default value = 1000
#
#OUT:   returns a list that contains the best model created and a sublist with the c-indexes for the train-,xTrain-set and the trees
#       grown by the best performing model
#-----------------------------------------------------------------------------------------------------------------------------------------

trainRF <-function(train_feat, train_obs, xTrain_feat, xTrain_obs, nTreeOffset=10, nTreeStep=10, maxTrees=1000) {
  
  # check if xTrain_feat has more than two elements and is of type "numeric"
  if(nrow(xTrain_feat)<2 || class(xTrain_feat)=="numeric"){
    stop("xTrain_feat contains less elements than 2. So it is not possible to calculate a Cindex value, that's why the
         code will be exited")
  }
  #######################################################################
  # check that there are no NA values in any committed parameter
  #######################################################################
  if(any(is.na(xTrain_feat))){
    stop("xTrain_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(xTrain_obs))){
    stop("xTrain_obs contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_feat))){
    stop("train_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_obs))){
    stop("train_obs contains NA values which are causing problems! This should never happen...")
  }
  
  # check if train_obs is unique!!! if it is unique throw Exception
  if(length(unique(train_obs))==1){
    stop("train_obs is unique!!! This should never happen...")
  }
  
  # check if there's variance in train_feat, if not throw Exception
  if(length(unique((as.vector(train_feat)))) == 1){
    stop("No variance in train_feat, all values are unique!")
  }
  
  # HINT:
  # concordance index cannot be calculated in case the xTrain observations are all the same. In this 
  # particular case, no rank can be build, and therefore, no cross-training is possible. As a solution,
  # the number of trees will be set to the max tree number (default=1000). 
  if(length(unique(xTrain_obs))==1) {
    return(list(model=randomForest(train_feat, y=train_obs, ntree=nTreeOffset), 
                performance=NA,
                nTree_best_performer=nTreeOffset,
                nTrees=NA))
  } 
  else {
    # set the parameter in charge of stopping the training to false, so that the training can start
    stopTrain <- FALSE
    # create a model, train it with the values of train_feat and calculate the c-index
    rf_old <- randomForest(train_feat, y=train_obs, ntree=nTreeOffset)
    cIdx_train <- cIDX(predict(rf_old, train_feat), train_obs)[1]
    
    # feed the model with the values from xtrain and calculate the c-index
    cIdx_xTrain <- cIDX(predict(rf_old, xTrain_feat), xTrain_obs)[1]
    # set cIdx_xTrain_old  to the value of cidx_xTrain 
    cIdx_xTrain_old <- cIdx_xTrain
    
    nTree <- nTreeOffset
    
    # train as long as stopTrain is false and the numberr of trees in the model is smaller than maxTrees
    while (!stopTrain  & rf_old$ntree < maxTrees) {
      # add nTreeStep trees to the existing model
      rf_new <- grow(rf_old, nTreeStep)
      # feed the model with the values from xtrain and calculate the c-index, define result as cIdx_xTrain_new
      cIdx_xTrain_new <-cIDX(predict(rf_new, xTrain_feat), xTrain_obs)[1]
      # only if there was an improved of the c-index between the old one and the new one
      if (cIdx_xTrain_old > cIdx_xTrain_new) {
        # stop training
        stopTrain <- TRUE
      } else {
        # set actual model to the previous one(rf_old) and the actual c-index of it to the previous one as well
        rf_old <- rf_new
        cIdx_xTrain_old <- cIdx_xTrain_new
      }
      # save all calculated c-indexes for the xTrain-set and the train-set in one list each
      cIdx_xTrain <- c(cIdx_xTrain, cIdx_xTrain_new)
      cIdx_train <- c(cIdx_train, cIDX(predict(rf_new, train_feat), train_obs)[1])
      
      nTree <- c(nTree, rf_new$ntree)
    }
  }
  
  return(list(model=rf_old, 
              performance=list(train_cIdx=cIdx_train, xTrain_cIdx=cIdx_xTrain),
              nTree_best_performer=rf_old$ntree,
              nTrees=nTree))
}

# #######################################################################################################################################
# ElasticNet Function
# #######################################################################################################################################

#-----------------------------------------------------------------------------------------------------------------------------------------
#DESC:  The function "trainEN" trains a "elasticNet" model with the data from "train_feat". 
#       This function will return a list containing the best model, the lambda calculated for the best model, 
#       the cIndex for the train and xTrain set and the weights of the features.
#
#IN:    train_feat      ==> train values obtained from crossvalidating the observated feature data.
#                           the model will be trained with this values
#       train_obs       ==> observed drug_response for the values of train_feat
#       xTrain_feat     ==> xTrain values obtained from crossvalidating the observed feature data.
#                           neccessary to find the best value for lambda
#       xTrain_obs      ==> observed drug_response for the values of xTrain_feat
#
# OUT:  model           ==> Elastic net model
#       lambda          ==> Chosen lambda for further predictions (parameter which has been optimize)
#       xTrain_cIDX     ==> Concordance indices of the cross-train data
#       train_cIDX      ==> Concordance indices of the train data
#       featWeights     ==> weigths for each feature of the chosen model
#-----------------------------------------------------------------------------------------------------------------------------------------

trainEN <- function(train_feat, train_obs, xTrain_feat, xTrain_obs) {
  
  # check if xTrain_feat has more than two elements and is of type "numeric"
  if(nrow(xTrain_feat)<2 || class(xTrain_feat)=="numeric"){
    stop("xTrain_feat contains less elements than 2. So it is not possible to calculate a Cindex value, that's why the
         code will be exited")
  }
  #######################################################################
  # check that there are no NA values in any committed parameter
  #######################################################################
  if(any(is.na(xTrain_feat))){
    stop("xTrain_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(xTrain_obs))){
    stop("xTrain_obs contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_feat))){
    stop("train_feat contains NA values which are causing problems! This should never happen...")
  }
  
  if(any(is.na(train_obs))){
    stop("train_obs contains NA values which are causing problems! This should never happen...")
  }
  
  # check if train_obs is unique!!! if it is unique throw Exception
  if(length(unique(train_obs))==1){
    stop("No variance in train_obs!!! This should never happen...")
  }
  
  # check if there's variance in train_feat, if not throw Exception
  if(length(unique((as.vector(train_feat)))) == 1){
    stop("No variance in train_feat, all values are unique!")
  }
  
  # HINT:
  # concordance index cannot be calculated in case the xTrain observations are all the same. In this 
  # particular case, no rank can be build, and therefore, no cross-training is possible. As a solution,
  # lambda is set to the first quantile of fit$lambda. 
  # alpha is used with default settings
  if(length(unique(xTrain_obs))==1) {
    
    fit <- glmnet(train_feat, train_obs)
    lambda <- fit$lambda[length(fit$lambda) - length(fit$lambda) / 4]                
    
    return(list(model=fit,
                lambda=lambda,
                xTrain_cIDX=NA,
                train_cIDX=NA,
                featWeights=fit$beta[,lambda==fit$lambda]))
  } else{
    # finding optimal EN mixing parameter (alpha) 
    # and optimal number of input features (lambda)
    xTrain_cIDX <- NA
    for (a in seq(1, 0, -0.1)) {
      # train an elasticNet model with the train_feat and train_obs
      current_fit <- glmnet(train_feat, train_obs, alpha=a)
      
      # predict with the xTrain_feat
      current_xTrain_pred <- predict(current_fit, type="response",newx=xTrain_feat)
      
      # find best lambda by calculating the xTrain_cidx
      current_xTrain_cIDX <- c()
      for (i in 1:ncol(current_fit$beta)) {
        current_xTrain_cIDX <- c(current_xTrain_cIDX, cIDX(current_xTrain_pred[,i], xTrain_obs))
      }
      
      # choose current model with this alpha if best fit
      if (length(xTrain_cIDX) == 1 || (max(xTrain_cIDX) < max(current_xTrain_cIDX))) {
        xTrain_cIDX <- current_xTrain_cIDX
        fit <- current_fit
        lambda <- current_fit$lambda[which(current_xTrain_cIDX == max(current_xTrain_cIDX))[1]]
      }
    }
    
    # train set performance across differnt lambda's
    train_pred <- predict(fit, type="response",newx=train_feat)
    train_cIDX <- c()
    for (i in 1:ncol(fit$beta)) {
      train_cIDX <- c(train_cIDX, cIDX(train_pred[,i], train_obs))
    }
  }
  
  return(list(model=fit,
              lambda=lambda,
              xTrain_cIDX=xTrain_cIDX,
              train_cIDX=train_cIDX,
              featWeights=fit$beta[,lambda==fit$lambda]))
}

######################
### main_bootstrap ###
######################

drug_res <- drug_res[drug_id,]
cells <- as.character(names(drug_res[!is.na(drug_res)]))

# using only cells with drug response and complete feature
cell_withFeat <- rownames(features)[!is.na(rowSums(features))]
cells <- intersect(cells, cell_withFeat)

# DON'T DO TRAINING WITH res<2*nFold and stop the script for this drug
if (length(cells) < nFolds){
  stop(paste("less than nFolds cell lines have a response for drug_id:", drug_id, sep=""))
} else{

  # define size of one bin for boot strapping (test set size)
  binSize <- ceil(length(cells) / nFolds)
  
  # for evaluating performance, save output
  predMat <- matrix(NA, ncol=binSize, nrow=nSteps)
  obsMat <- matrix(NA, ncol=binSize, nrow=nSteps)
  cellIDmat <- matrix(NA, ncol=binSize, nrow=nSteps)
  
  # create a matrix "impMat" with the featurenames of the features matrix
  impMat <- matrix(NA, ncol=ncol(features), nrow=nSteps)
  colnames(impMat) <- colnames(features)
  
  # train 1000 models if possible
  for (i in 1:(nSteps+1)){

    # train as long as i is smaller than 1000
    if(i<=nSteps){
      
      # sample COSMIC ID'S
      ids <- sample(cells)
         
      # devide ids into three samples
      test <- ids[1:binSize]
      xTrain <- ids[(binSize+1):(2*binSize)]
      train <- ids[(2*binSize+1):length(ids)]
      
      #####################################################################################################################
      # TEST section for the bootstrapping
      #####################################################################################################################
      
      #TEST if all COSMIC ID's have been used
      if(sum(length(test),length(train),length(xTrain))!=length(cells)) 
          stop(paste("not all COSMIC ID's are either used in train, xTrain or test in step", i, sep=""),call. =FALSE)
      
      #TEST that there are no equal COSMIC ID's in train and xtrain
      if(length(intersect(train,xTrain))>0) 
        stop(paste("train and xTrain contain the same values in step:", i, sep=""))
      
      #TEST that there are no equal COSMIC ID's in train and test
      if(length(intersect(train,test))>0) 
        stop(paste("train and test contain the same values in step:", i, sep=""))
      
      #TEST that there are no equal COSMIC ID's in xTrain and test
      if(length(intersect(xTrain,test))>0) 
        stop(paste("xTrain and test contain the same values in step:", i, sep=""))
      
      #TEST if there's variance in the train set
      if(length(unique(drug_res[as.character(train)]))==1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      #TEST if there's variance in the xTrain set
      if(length(unique(drug_res[as.character(xTrain)]))==1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      ##############################################################################################################
      # Train statistical model
      ##############################################################################################################
      
      # prepare trainset
      train_obs <- drug_res[as.character(train)]
      # check if trainset should be balanced or not
      if(balancing){
        train_obs <- balance(train_obs)
      }
      train_feat <- features[names(train_obs), ]
      rownames(train_feat) <- names(train_obs)
      
      # check if there's variance in train_feat 
      if(length(unique((as.vector(train_feat)))) == 1){
        # go to the next model for this drug with new sampled cell lines
        next
      }
      
      # prepare cross-trainset
      xTrain_obs <- drug_res[as.character(xTrain)]
      xTrain_feat <- features[names(xTrain_obs), ]
      rownames(xTrain_feat) <- names(xTrain_obs)
      
      # prepare testset
      test_obs <- drug_res[as.character(test)]
      test_feat <- features[names(test_obs), ]
      rownames(test_feat) <- names(test_obs)
      
      # #######################################################################
      # Train a randomForest model
      # #######################################################################
      if(stat_model=="RF"){
        # train a randomForest model
        trainOut <- trainRF(train_feat, train_obs, xTrain_feat, xTrain_obs)
        
        # predict with test set
        test_pred <- predict(trainOut$model, test_feat)
        
        # define best model importance
        best_model_importance <- trainOut$model$importance
      }
      
      # #######################################################################
      # Train an elasticNet model
      # #######################################################################
      if(stat_model=="EN"){
        # train an elasticNet model
        trainOut  <- trainEN(train_feat, train_obs, xTrain_feat, xTrain_obs)
        
        # predict with test set
        test_pred <- predict(trainOut$model,type="response",newx=test_feat)[,trainOut$model$lambda == trainOut$lambda]
        
        # define best model importance
        best_model_importance <- as.matrix(trainOut$featWeights) * apply(features[cells,], 2, sd)
      }
      
      # #######################################################################
      # store output matrix
      # #######################################################################
      
      predMat[i, ] <- test_pred
      obsMat[i, ] <- test_obs
      cellIDmat[i, ] <- names(test_obs)
    
      # #######################################################################
      # Fill feat importance matrix
      # #######################################################################
 
      impMat[i, rownames(best_model_importance)] <- best_model_importance
          
    } else{  
      # #######################################################################
      # Create folder structure for each drug
      # #######################################################################
      
      # calculate the mean and the standard deviation of the impMat
      vMean <- apply(impMat, 2, mean)
      vStd <- apply(impMat, 2, sd)
      
      # save the standard deviation and the mean in a matrix
      effectSize <- cbind(vMean, vStd)
      
      path_to_predMat <- file.path(o, 'predMat.tsv')
      write.table(predMat, path_to_predMat, sep='\t', row.names=F, col.names=F, quote=F)
    
      path_to_obsMat <- file.path(o, 'obsMat.tsv')
      write.table(obsMat, path_to_obsMat, sep='\t', row.names=F, col.names=F, quote=F)
    
      path_to_cellIDmat <- file.path(o, 'cellIDmat.tsv')
      write.table(cellIDmat, path_to_cellIDmat, sep='\t', row.names=F, col.names=F, quote=F)
    
      path_to_effectSize <- file.path(o, 'effectSize.tsv')
      write.table(t(c('vFeature', colnames(effectSize))), path_to_effectSize, sep='\t', row.names=F, col.names=F, quote=F)
      write.table(effectSize, path_to_effectSize, sep='\t', row.names=T, col.names=F, quote=F, append=T)
      
      if (!all(is.na(predMat))) {
          performance <- cor(as.vector(predMat), as.vector(obsMat), method='pearson', use="complete.obs")
          pVal <- cor.test(as.vector(predMat), as.vector(obsMat), method='pearson', use="complete.obs")$p.value
          VAR <- var(as.vector(predMat))
          STE <- sd(predMat, na.rm=T) / sqrt(prod(dim(cellIDmat)))
          # pooled variance
          frame <- data.frame(cell_id=as.vector(cellIDmat), pred=as.vector(predMat))
          nCell <- aggregate(pred~cell_id, frame,  length)
          variance <- aggregate(pred~cell_id, frame,  var)
          POOLED_VAR <- sum((nCell$pred-1) * variance$pred) / (sum(nCell$pred) - nrow(nCell))
          POOLED_STE <- sqrt(POOLED_VAR) / sqrt(prod(dim(cellIDmat)))

          perfTab <- data.frame(Variate = c('performance', 'pVal', 'VAR', 'STE', 'POOLED_STE', 'POOLED_VAR'), 
                            Value = c(performance, pVal, VAR, STE, POOLED_STE, POOLED_VAR))
          colnames(perfTab)[2] = drug_id
          path_to_perfTab <- file.path(o, 'perfTab.tsv')
          write.table(perfTab, path_to_perfTab, sep='\t', row.names=F, col.names=T, quote=F)
      } else {
          print("All elements of 'predMat' are NA")
      }
    }
  }
}

# save the options
optTab <- data.frame(Options = names(unlist(opt)), Values = unlist(opt))
colnames(optTab)[2] = drug_id
path_to_optTab <- file.path(o, 'optTab.tsv')
write.table(optTab, path_to_optTab, sep='\t', row.names=F, col.names=T, quote=F)

print('finished!')

### THE END ###
