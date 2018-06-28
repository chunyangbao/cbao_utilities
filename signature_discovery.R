#!/usr/bin/env Rscript

library("optparse")
library("data.table")

V="Version: 1.0"
D="Depends: R (>= 3.1.0), optparse, data.table, gplots, ggplot2, reshape2"

a = commandArgs(trailingOnly=TRUE)

option_list = list(
    make_option(c("-o", "--dir_out"), type="character", default=".",
     help="Output directory [%default] \n\t\tBy default <o> is '.' (current directory).", metavar="character"),
    # somatic_discovery
    make_option(c("--n_iter"), type="integer", default="50",
     help="Number of independent simulations [%default] ", metavar="integer"),
    make_option(c("--kcol"), type="integer", default="96",
     help="Number of initial signatures [%default]", metavar="integer"),
    make_option(c("--tol"), type="character", default="1.e-07",
     help="Tolerance for convergence [%default]", metavar="character"),
    make_option(c("--a0"), type="integer", default="10",
     help="hyper-parameter [%default]", metavar="integer"),
    make_option(c("--tumor_type"), type="character", default="tumor.type",
     help="Please specify your cohort name here [%default]", metavar="character"),
    make_option(c("--prior"), type="character", default="L1KL",
     help="Choose pirors for W and H [%default] \n\t\tL1KL (expoential priors); L2KL (half-normal priors)", metavar="character"),
    make_option(c("--hyper"), type="logical", default="FALSE", action = "store_true",
     help="Reduce the effect of hyper-mutant samples in the signature discovery? [%default]", metavar="logical"),
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

opt_parser = OptionParser(usage="usage: %prog [options] <lego96.tsv>\n\tBy default <lego96.tsv> is '-' (stdin). This matrix should contain \n\tmutation counts along 96 tri-nucleotide mutation contexts (rows) \n\tacross samples (columns). Rownames of the lego matrix should be \n\t4-letters. ex: CGTA (C to G mutation at 5'-TCA-3'contexts)", option_list=option_list, description = paste(V, D, sep="\n"))
opt = parse_args(opt_parser, args=a, positional_arguments=TRUE)

o = opt$options$dir_out
n.iter = opt$options$n_iter
Kcol = opt$options$kcol
tol = as.numeric(opt$options$tol)
a0 = opt$options$a0
tumor.type = opt$options$tumor_type
prior = opt$options$prior
hyper = opt$options$hyper

cs = opt$options$col_sep
rn = opt$options$row_num
cn = as.logical(opt$options$col_name)
na = opt$options$na_str
rs = opt$options$row_skip
ra = opt$options$row_name

m = opt$args[1]

if (is.na(m)|m=='-') m = 'file:///dev/stdin'

if(grepl("^[[:digit:]]+$", rs)) rs=as.numeric(rs)
if(grepl("^[[:digit:]]+$", ra)) ra=as.numeric(ra)
# read
fm = fread(m, sep=cs, nrows=rn, header=cn, na.string=na, skip=rs, data.table=FALSE, check.names=FALSE)
rownames(fm) = fm[,ra]
fm = as.matrix(fm[,-ra])
lego96 = fm

library(gplots)
library(ggplot2)
library(reshape2)

OUTPUT <- paste0(o, "/")
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)

######################################################################################################
####### Mutation Signature Profiling using Bayesian NMF algorithms 
######################################################################################################
####### For details on the implementation 
####### see J Kim, Mouw K, P Polak et al, Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors 
####### Nat. Genet. DOI: 10.1038/ng.3557
####### For details on the original algorithms 
####### see Tan, V.Y. & Févotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
####### IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).
######################################################################################################

##########################################################
###################### R functions  ######################
##########################################################

## Bayesian NMF algorithm with exponential priors for W and H
BayesNMF.L1KL <- function(V0,n.iter,a0,tol,K,K0,phi) {
        eps <- 1.e-50
        del <- 1.0
        active_nodes <- colSums(V0) != 0
        V0 <- V0[,active_nodes]
        V <- V0-min(V0) + eps
        Vmin <- min(V)
        Vmax <- max(V)
        N <- dim(V)[1]
        M <- dim(V)[2]
        W <- matrix(runif(N * K)*sqrt(Vmax),ncol=K)
        H <- matrix(runif(M * K)*sqrt(Vmax),ncol=M)
        V.ap <- W %*% H + eps
        I <- array(1,dim=c(N,M))

        C <- N + M + a0 + 1
        b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
        lambda.bound <- b0/C
        lambda <- (colSums(W)+rowSums(H)+b0)/C
        lambda.cut <- 1.5*lambda.bound

        n.like <- list()
        n.evid <- list()
        n.error <- list()
        n.lambda <- list()
        n.lambda[[1]] <- lambda

        iter <- 2
        count <- 1
        while ((del >= tol) & (iter < n.iter)) {
                H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+phi/lambda,M),ncol=M) + eps)
                V.ap <- W %*% H + eps
                W <- W * ((V/V.ap) %*% t(H))/t(matrix(rep(rowSums(H)+phi/lambda,N),ncol=N) + eps)
                V.ap <- W %*% H + eps
                lambda <- (colSums(W) + rowSums(H) + b0) / C
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                like <- sum(V*log(V/V.ap)+V.ap-V)
                n.like[[iter]] <- like
                n.evid[[iter]] <- like+phi*sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
                n.lambda[[iter]] <- lambda
                n.error[[iter]] <- sum((V-V.ap)^2)
                if (iter %% 100 == 0) {
                        cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
                }
                iter <- iter+1
        }
        return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

## Bayesian NMF algorithm with hlaf-normal priors for W and H
BayesNMF.L2KL <- function(V0,n.iter,a0,tol,K,K0,phi) {
        eps <- 1.e-50
        del <- 1.0
        active_nodes <- colSums(V0) != 0
        V0 <- V0[,active_nodes]
        V <- V0-min(V0)
        Vmax <- mean(V)
        N <- dim(V)[1]
        M <- dim(V)[2]

        W <- matrix(runif(N * K)*Vmax,ncol=K)
        H <- matrix(runif(M * K)*Vmax,ncol=M)
        V.ap <- W%*%H+eps

        C <- (N+M)/2+a0+1
        b0 <- 3.14*(a0-1)*mean(V)/(2*K0)
        lambda <- (0.5*colSums(W*W)+0.5*rowSums(H*H)+b0)/C
        lambda.bound <- b0/C
        lambda.cut <- b0/C*1.25

        I <- array(1,dim=c(N,M))
        n.like <- list()
        n.evid <- list()
        n.error <- list()
        n.lambda <- list()
        n.lambda[[1]] <- lambda
        iter <- 2
        count <- 1
        while ((del >= tol) & (iter < n.iter)) {
                H <- H*(t(W)%*%(V/V.ap)/(t(W)%*%I+phi*H/matrix(rep(lambda,M),ncol=M)+eps))^0.5
                V.ap <- W%*%H+eps
                W <- W*((V/V.ap)%*%t(H)/(I%*%t(H)+phi*W/t(matrix(rep(lambda,N),ncol=N))+eps))^0.5
                lambda <- (0.5*colSums(W*W)+0.5*rowSums(H*H)+b0)/C
                V.ap <- W%*%H+eps
                del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
                like <- sum(V * log((V+eps)/(V.ap+eps)) + V.ap - V)
                n.like[[iter]] <- like
                n.evid[[iter]] <- like + phi*sum((0.5*colSums(W^2)+0.5*rowSums(H^2)+b0)/lambda+C*log(lambda))
                n.lambda[[iter]] <- lambda
                n.error[[iter]] <- sum((V-V.ap)^2)
                if (iter %% 100 == 0) {
                        cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
                }
                iter <- iter+1
        }
        return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

######### Handling hypermutant samples; See J Kim et al Nat. Genet. DOI: 10.1038/ng.3557 for details.
get.lego96.hyper <- function(lego96) {
	x <- lego96
	for (i in 1:100) {
	        SNV <- colSums(x)
	        q1 <- quantile(SNV,prob=1/4)
	        q3 <- quantile(SNV,prob=3/4)
	        sample.hyper <- colnames(x)[SNV > (median(SNV)+1.5*(q3-q1))]
	        if (length(sample.hyper)==0) break
	       	lego96.hyper <- as.matrix(x[,(colnames(x) %in% sample.hyper)])
	       	colnames(lego96.hyper) <- sample.hyper
	       	lego96.nonhyper <- x[,!(colnames(x) %in% sample.hyper)]
	        lego96.hyper1 <- apply(lego96.hyper,2,function(x) x/2)
	        lego96.hyper2 <- lego96.hyper1
	        colnames(lego96.hyper1) <- paste(colnames(lego96.hyper1),1,sep="__")
	        colnames(lego96.hyper2) <- paste(colnames(lego96.hyper2),2,sep="__")
	        x <- cbind(lego96.nonhyper,lego96.hyper1,lego96.hyper2)
	}
	return(x)
}

######### Visualizing the activity of signatures
plot.activity.barplot <- function(H.mid,H.norm,scale,tumor.type) {
        .theme_ss <- theme_bw(base_size=14) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
                axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
                axis.text = element_text(size = 12*scale, family = "mono"))
        ordering <- order(colSums(H.mid),decreasing=T)
        H.mid <- H.mid[,ordering]
        rownames(H.mid) <- paste("W",seq(1:nrow(H.mid)),sep="")
        H.norm <- H.norm[,ordering]
        rownames(H.norm) <- paste("W",seq(1:nrow(H.norm)),sep="")
	sample.ordering <- colnames(H.mid)
        x1 <- melt(H.mid)
        x2 <- melt(H.norm)
        colnames(x1) <- c("Signature","Sample","Activity")
        colnames(x2) <- c("Signature","Sample","Activity")
        x1[,"class0"] <- c("Counts")
        x2[,"class0"] <- c("Fractions")
        df2 <- rbind(x1,x2)
        df2$class0 <- factor(df2$class0,c("Counts","Fractions"))
	df2$Sample <- factor(df2$Sample,sample.ordering)
        scale <- 1
        p = ggplot(df2,aes(x=factor(Sample),y=Activity,fill=factor(Signature)))
        p = p+geom_bar(stat="identity",position='stack',color='black',alpha=0.9)
        p = p + scale_fill_manual(values=c("red","cyan","yellow","blue","magenta","gray50","orange","darkgreen","brown","black",rainbow(10)[4:10]))
        p = p + facet_grid(class0 ~ ., scale = "free_y")
        p = p + ggtitle(paste("Siganture Activities in",tumor.type,sep=" "))
        p = p + theme(plot.title=element_text(lineheight=1.0,face="bold",size=14*scale))
        p = p + xlab("Samples") + ylab("Signature Activities")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.text.x = element_text(angle=90,vjust=0.5,size=8*scale,face="bold",colour="black"))
        p = p + theme(axis.text.y = element_text(size=10*scale,face="bold",colour="black"))
        p = p + theme(legend.title=element_blank())
        p = p + .theme_ss
        p = p + theme(legend.position="top")
        return(p)
}

######### Collecting data from several independent runs
get.stats.simulation <- function(tumor,n.iter,OUTPUT) {
        W.ALL <- list()
        H.ALL <- list()
        K.ALL <- list()
        lambda.ALL <- list()
        evid.ALL <- list()
        for (i in 1:n.iter) {
                x <- load(paste(OUTPUT,paste(method,a0,i,"RData",sep="."),sep=""))
                W <- res[[1]]
                H <- res[[2]]
                evid <- res[[4]]
                evid <- evid[[length(evid)]]
                lambda <- res[[5]]
                lambda <- lambda[[length(lambda)]]
                lambda.min <- min(lambda)
                lambda <- lambda-lambda.min
                lambda.norm <- lambda/sum(lambda)
                index <- lambda.norm > 0.01
                K <- sum(index)
                W.ALL[[i]] <- W
                H.ALL[[i]] <- H
                lambda.ALL[[i]] <- lambda
                K.ALL[[i]] <- K
                evid.ALL[[i]] <- evid
        }
        return(list(W.ALL,H.ALL,lambda.ALL,K.ALL,evid.ALL))
}

######### Visualizing the profile of signatures
plot.lego.observed.barplot <- function(mat,title) {
        scale <- 1.5
        .theme_ss <- theme_bw(base_size=20) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=10*scale, family="mono"),
                axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
                axis.text = element_text(size = 16*scale, family = "mono"))
        p = ggplot(mat)
        p = p + geom_bar(aes_string(x="context3",y="signature",fill="base.sub"),stat="identity",position="identity",colour="gray50")
        p = p + facet_grid(class ~ base.sub, scale = "free_y")
        p = p + .theme_ss
        p = p + scale_fill_manual(values=c("cyan","red","yellow","purple","green","blue","black","gray")) #p = p + scale_fill_brewer(palette = "Set1")
        p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
        p = p + ggtitle(title)
        p = p + xlab("Motifs") + ylab("Contributions")
        p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
        p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
        return(p)
}

get.df.solution <- function(W,H,lambda,tumor) {
        lambda <- lambda[[length(lambda)]]
        lambda.min <- min(lambda)
        lambda <- (lambda-lambda.min)/max(lambda)
        norm.lambda <- lambda/sum(lambda)
        index <- lambda > 0.01
        x <- colSums(W)/max(colSums(W))
        index <- x > 0.01
        K <- sum(index)
        if (K != 1) {
                W <- W[,index]
                H <- H[index,]
                norm.W <- t(apply(W,1,function(x) x/sum(x)))
                norm.H <- apply(H,2,function(x) x/sum(x))
        }
        colnames(W) <- paste("W",seq(1:ncol(W)),sep="")
        rownames(H) <- paste("H",seq(1:ncol(W)),sep="")

        context96 <- rownames(W)[1:96]
        context96 <- sub("->","",context96)
        context96 <- gsub("[.]","",context96)
        index96 <- c(order(substring(context96,1,2),decreasing=F))
        W.SNP <- as.vector(W[index96,])
        seq <- c(context96[index96])
        context4 <- rep(seq,K)
        context3 <- paste(substring(context4,3,3),"-",substring(context4,4,4),sep="")
        base.sub <- paste(substring(context4,1,1),"->",substring(context4,2,2),sep="")
        signature.class <- as.vector(t(matrix(t(rep(1:K,96)),nrow=K)))
        signature.class <- paste("W",signature.class,sep="")
        df <- data.frame(W.SNP,context4,context3,base.sub,signature.class,rep(tumor,nrow(W)*ncol(W)))
        colnames(df) <- c("signature","context4","context3","base.sub","class","tumor")
        return(df)
}
##########################################################
##########################################################
##########################################################

##################################
############### Choose pirors for W and H
############### Default = L1KL (expoential priors); L2KL (half-normal priors)
##################################
if (prior=="L1KL") {
	method <- paste("L1KL.lego96",tumor.type,sep=".")
} else {
	method <- paste("L2KL.lego96",tumor.type,sep=".")
}
##################################

##################################
############## Default = FALSE ; TRUE - to reduce the effect of hyper-mutant samples in the signature discovery
##################################
if (hyper) {
	lego96 <- get.lego96.hyper(lego96)
	method <- paste(method,"hyper",sep=".")
}
##################################

##########################################################
###################### Running the algorithms ############
##########################################################
for (i in 1:n.iter) {
	if (prior=="L1KL") {
		res <- BayesNMF.L1KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
	} else {
		res <- BayesNMF.L2KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
	}
	save(res,file=paste(OUTPUT,paste(method,a0,i,"RData",sep="."),sep=""))
}

##########################################################
###################### Analysis ##########################
##########################################################
res.WES <- get.stats.simulation(tumor.type,n.iter,OUTPUT)

############## frequency figure 
pdf(file=paste(OUTPUT,paste(method,a0,"signature.freq.pdf",sep="."),sep=""),width=4,height=5)
s1 <- 1.5
s2 <- 2.0
par(mfrow=c(1,1))
par(mar=c(5,5,2,1))
        barplot(table(unlist(res.WES[[4]])),cex=s1,cex.axis=s1,cex.main=s1,cex.names=s1,cex.lab=s1,xlab="# of signatures",ylab="Freq.",main=paste(tumor.type,sep="."))
dev.off()

signature.freq = table(unlist(res.WES[[4]]))
tsf = as.data.frame(signature.freq,stringsAsFactors=F)
colnames(tsf) = c("n_sigs","count")
msf = tsf[which(tsf[,2]==max(tsf[,2])),1]
write.table(tsf,file=paste(OUTPUT,paste(method,a0,"signature.freq",msf,"tsv",sep="."),sep=""),sep='\t',row.names=F,col.names=T,quote=F)

##########################################################
############## select the best solution (maximum posteria solution) for given K
##########################################################
tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
        tmpX <- res.WES[[5]]
        tmpX[tmpK != unique.K[i]] <- 0
        MAP[[i]] <- which.min(tmpX)
}

tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
        tmpX <- res.WES[[5]]
        tmpX[tmpK != unique.K[i]] <- 0
        MAP[[i]] <- which.min(tmpX)
}
MAP <- unlist(MAP)
names(MAP) <- unique.K
MAP.nontrivial <- MAP[names(MAP)!=1]
##########################################################
##########################################################

n.K <- length(MAP.nontrivial)
if (n.K > 0) { 
for (j in 1:n.K) {
        load(paste(OUTPUT,paste(method,a0,MAP.nontrivial[[j]],"RData",sep="."),sep=""))
        W <- res[[1]]
        H <- res[[2]]
        W1 <- W
        H1 <- H
        W.norm <- apply(W,2,function(x) x/sum(x))
        for (i in 1:ncol(W)) {
                H1[i,] <- H1[i,]*colSums(W)[i]
                W1[,i] <- W1[,i]*rowSums(H)[i]
        }
        lambda <- res[[5]]
        df <- get.df.solution(W1,H,lambda,tumor.type)
        K <- length(table(df$class))

        ############# Signature plot
        width <- 16
        height <- ifelse(K==1,3,K*2)
        pdf(file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"signature.pdf",sep="."),sep=""),width=width,height=height)
                p <- plot.lego.observed.barplot(df,paste("Mutation Signatures in ",tumor.type,sep=""))
                plot(p)
        dev.off()

        index.nonzero <- colSums(W) != 0
        lambda <- unlist(res[[5]][length(res[[5]])])
        lambda <- lambda/min(lambda)
        lambda <- lambda[index.nonzero]
        names(lambda) <- paste("W",seq(1:length(lambda)),sep="")
        K <- sum(index.nonzero)
        if (K != 1) {
                W0.tumor <- W[,index.nonzero]
                H0.tumor <- H[index.nonzero,]
                x <- W0.tumor
                y <- H0.tumor
                x.norm <- apply(x,2,function(x) x/sum(x))
                W1.tumor <- x.norm
                for (j in 1:K) y[j,] <- y[j,]*colSums(x)[j]
                H1.tumor <- y
                y.norm <- apply(y,2,function(x) x/sum(x))
                H2.tumor <- y.norm
        } else {
                stop("No non-trivial solutations; All simulations converged to a trivial solution with one signature")
        }

        ############# Reconstructing the activity of signatures
        W.mid <- W1.tumor
        H.mid <- H1.tumor
        H.norm <- H2.tumor
        if (length(grep("__",colnames(H.mid))) != 0) {
                hyper <- colnames(H.mid)[grep("__",colnames(H.mid))]
                H.hyper <- H.mid[,colnames(H.mid) %in% hyper]
                H.nonhyper <- H.mid[,!(colnames(H.mid) %in% hyper)]
                sample.hyper <- colnames(H.hyper)
                sample.hyper <- sapply(sample.hyper,function(x) strsplit(x,"__")[[1]][[1]])
                unique.hyper <- unique(sample.hyper)
                n.hyper <- length(unique.hyper)
                x.hyper <- array(0,dim=c(nrow(H.hyper),n.hyper))
                for (i in 1:n.hyper) {
                        x.hyper[,i] <- rowSums(H.hyper[,sample.hyper %in% unique.hyper[i]])
                }
                colnames(x.hyper) <- unique.hyper
                rownames(x.hyper) <- rownames(H.mid)
                H.mid <- (cbind(H.nonhyper,x.hyper))
                H.norm <- apply(H.mid,2,function(x) x/sum(x))
        }
	W.norm <- apply(W.mid,2,function(x) x/sum(x))

	##########################################################
        ############# W.norm = extracted signatures normalized to one
	############# H.mid = activity of signatures across samples (expected mutations associated with signatures)
	############# H.norm = normalized signature activity 
	##########################################################
        WH <- list(W.norm,H.mid,H.norm)
        save(WH,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"WH.RData",sep="."),sep=""))

        ############# Activity plot
	p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
	pdf(file = paste(OUTPUT,paste(method,a0,"activity.barplot1",K,"pdf",sep="."),sep=""),width=15,height=12)
	        plot(p1)
	dev.off()
    
    ### write (by Chunyang Bao) ###
    W_norm = W.norm
    H_mid = H.mid
    H_norm = H.norm
    
    colnames(W_norm) = paste0("W",seq(1:ncol(W_norm)))
    rownames(H_mid) = paste0("W",seq(1:nrow(H_mid)))
    rownames(H_norm) = paste0("W",seq(1:nrow(H_norm)))
    
    write.table(t(c('motifs',colnames(W_norm))),file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"W_norm.tsv",sep="."),sep=""),sep='\t',row.names=F,col.names=F,quote=F)
    write.table(W_norm,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"W_norm.tsv",sep="."),sep=""),sep='\t',row.names=T,col.names=F,quote=F,append=T)
    write.table(t(c('signatures',colnames(H_mid))),file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"H_mid.tsv",sep="."),sep=""),sep='\t',row.names=F,col.names=F,quote=F)
    write.table(H_mid,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"H_mid.tsv",sep="."),sep=""),sep='\t',row.names=T,col.names=F,quote=F,append=T)
    write.table(t(c('signatures',colnames(H_norm))),file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"H_norm.tsv",sep="."),sep=""),sep='\t',row.names=F,col.names=F,quote=F)
    write.table(H_norm,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"H_norm.tsv",sep="."),sep=""),sep='\t',row.names=T,col.names=F,quote=F,append=T)

}
}



### THE END ###
