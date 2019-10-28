# a comment

##' getSI function
##'
##' An internal function to compute susceptible/infected given a matrix of infection times: assumes there are no ties in these times.
##'
##' @param mat matrix of infection times: assumes there are no ties in these times.
##' @return a list of infected/susceptibles at each time point

getSI <- function(mat){
    m <- mat[mat>0]
    o <- order(m,decreasing=TRUE)
    m <- m[o]
    n <- length(m)

    S <- list()
    I <- list()
    IDX <- rep(NA,n)
    for(i in 1:n){
        S[[i]] <- which(mat<m[i])
        I[[i]] <- which(mat>=m[i])
        if(i<n){
            IDX[i] <- which(S[[i]]==which(mat==m[i+1])) # which of the susceptibles became infected
        }
    }
    return(list(S=S,I=I,idx=IDX))
}


##' gammaFun function
##'
##' An internal two-parameter function to evaluate the measure of dependence between cells a given distance apart. Function must work on vector inputs.
##'
##' On the to-do list is extend this to more general functions with an arbitrary number of parameters
##'
##' @param dst distance between cells
##' @param theta paramter theta
##' @param phi parameter phi
##' @return the function gamma evaluated at those distances

gammaFun <- function(dst,theta,phi){
    return(theta + exp(-dst/phi))
    #return(theta + 1/(phi*((dst^2+1)/10)^2))
    #return(theta + exp(-dst^2/phi))
}


##' distFun function
##'
##' An internal function to create a function to evaluate the transmission probabilities from cell i to cell j - includes model for susceptibility and also infectivity, as well as distance decay
##'
##' @param mat matrix of hcv levels vectorised using function as.vector, see prepData function
##' @return a function to evaluate the transmission probabilities from cell i to cell j
distFun <- function(mat){
    sqR <- as.integer(sqrt(length(mat)))
    dm <- as.matrix(dist(expand.grid(1:sqR,1:sqR)))

    f <- function(i,j,theta,phi,beta,inf,sus,N){

        gamma <- gammaFun(dm[i,j,drop=FALSE],theta,phi)

        if(sum(N)==0){
            return(gamma)
        }
        else{

            alpha_inf <- rep(1,length(i))
            alpha_sus <- rep(1,length(j))
            if(N[1]>0){
                alpha_inf <- as.vector(exp(inf[i,,drop=FALSE]%*%beta[1:N[1]]))
            }
            if(N[2]>0){
                alpha_sus <- as.vector(exp(sus[j,,drop=FALSE]%*%beta[(N[1]+1):sum(N)]))
            }

            out <- outer(alpha_inf,alpha_sus,FUN="*")
            return(out * gamma)

        }

    }

    return(f)
}


##' partrans function
##'
##' An internal function to backtransform the parameters
##'
##' @param param a vector of transformed parameters (theta,phi,beta_1,...,beta_d)
##' @return a vector of back-transformed parameters

partrans <- function(param){
    param[1:2] <- exp(param[1:2])
    return(param)
}


##' logPL function
##'
##' An internal function to calculate the partial likelihood
##'
##' @param param vector of parameters on the transformed scale (note they are back-transformed as part of this routine)
##' @param si list of susceptible/infected, see ?getSI
##' @param distfun function created using the function distFun
##' @param sus covariate data relating to susceptibility
##' @param inf covariate data relating to infectivity
##' @param N vector of length 2: the number of covariates associated with respectively infectivity and susceptibility
##' @return the log-partial likelihood evaluated at the input parameters

logPL <- function(param,si,distfun,inf,sus,N){
    param <- partrans(param)
    theta <- param[1]
    phi <- param[2]
    beta <- param[-c(1:2)]
    nt <- length(si$S)
    rho_mats <- lapply(1:nt,function(i){return(distfun(si$I[[i]],si$S[[i]],theta,phi,beta,inf,sus,N))})
    den <- sapply(rho_mats,sum)
    num <- sapply(1:(nt-1),function(i){return(sum(rho_mats[[i]][,si$idx[i]]))})
    return(sum(log(num/den[1:(nt-1)])))
}


##' asHCVdata function
##'
##' A function to convert a dataset into one that can be used with HCV's internal/visible functions
##'
##' @param dat  a dataframe with the HCV data, sorted in 'expand.grid' order within patient and sample number
##' @param hcv a character string, the name of the variable containing the HCV levels
##' @param subjectid a character string, the name of the variable containing the subject identifier
##' @param gridnum a character string, the name of the variable containing the sample number for each individual
##' @param x a character string, the name of the variable containing the x-index of the grid
##' @param y a character string, the name of the variable containing the y-index of the grid
##' @return ...
##' @export

asHCVdata <- function(dat,hcv,subjectid,gridnum,x,y){
    names(dat)[which(names(dat)==hcv)] <- "hcv"
    names(dat)[which(names(dat)==subjectid)] <- "subjectid"
    names(dat)[which(names(dat)==gridnum)] <- "gridnum"
    names(dat)[which(names(dat)==x)] <- "x"
    names(dat)[which(names(dat)==y)] <- "y"
    class(dat) <- c("data.frame","HCVdataset")
    return(dat)
}


##' prepData function
##'
##' An internal function to prepare data into a form ready for analysis. Many computations with this model can be performed once at the start and stored in memory in order to increase computational speed; this function manages that process.
##'
##' @param dat a HCVdataset object, see ?asHCVdata
##' @param iForm formula for infectivity
##' @param sForm formula for susceptibility
##' @param inits initial values for optimiser, the default is NULL, where initial values are chosen automatically (in an ad-hoc fashion)
##' @return the data and model inputs ready to run the optimiser.

prepData <- function(dat,iForm=NULL,sForm=NULL,inits=NULL){

    if(!inherits(dat,"HCVdataset")){
        stop("dat must inherit class 'HCVdataset', see ?asHCVdata")
    }

    sid <- sort(unique(dat$subjectid))
    nsub <- length(sid)
    gid <- list()
    for(i in 1:length(sid)){
        gid[[i]] <- sort(unique(dat$gridnum[dat$subject==sid[i]]))
    }

    if(is.null(inits)){
        inits <- c(0,0)
    }

    matList <- list()
    covList_inf <- list()
    covList_sus <- list()

    count <- 1

    for(i in 1:nsub){
        for(j in 1:length(gid[[i]])){

            ssid <- dat$subjectid==sid[i] & dat$gridnum==gid[[i]][j]
            matList[[count]] <- dat$hcv[ssid]

            covList_inf[[count]] <- 1
            covList_sus[[count]] <- 1

            N <- c(0,0)

            if(!is.null(iForm)){
                covList_inf[[count]] <- model.matrix(iForm,dat[ssid,])
                covList_inf[[count]] <- covList_inf[[count]][,-1,drop=FALSE] # remove intercept
                if(count==1){
                    inits <- c(inits,rep(0,ncol(covList_inf[[count]])))
                    N[1] <- ncol(covList_inf[[count]])
                }
            }

            if(!is.null(sForm)){
                covList_sus[[count]] <- model.matrix(sForm,dat[ssid,])
                covList_sus[[count]] <- covList_sus[[count]][,-1,drop=FALSE] # remove intercept
                if(count==1){
                    inits <- c(inits,rep(0,ncol(covList_sus[[count]])))
                    N[2] <- ncol(covList_sus[[count]])
                }
            }

            count <- count + 1

        }
    }

    return(list(N=N,inits=inits,matList=matList,covList_inf=covList_inf,covList_sus=covList_sus))

}

##' optimPL function
##'
##' A function to maximise the partial likelihood.
##'
##' @param dat a HCVdataset object, see ?asHCVdata
##' @param iForm formula for infectivity
##' @param sForm formula for susceptibility
##' @param inits initial values to pass to optim
##' @param pd optional prepData output, used in function bootCI
##' @param ... addition arguments passed to optim function
##' @return output from optim, run on the (log) partial likelihood function
##' @export
optimPL <- function(dat,iForm=NULL,sForm=NULL,inits=NULL,pd=NULL,...){

    if(is.null(pd)){
        pd <- prepData(dat=dat,iForm=iForm,sForm=sForm,inits=inits)
    }

    si <- lapply(pd$matList,getSI)
    df <- lapply(pd$matList,distFun)
    n <- length(pd$matList)

    logPLmulti <- function(param,si,distfun,inf,sus){
        param <- rep(list(param),n)
        ans <- mapply(logPL,param,si,distfun,inf,sus,MoreArgs=list(N=pd$N))
        return(sum(ans))
    }

    return(optim(pd$inits,logPLmulti,si=si,distfun=df,inf=pd$covList_inf,sus=pd$covList_sus,control=list(fnscale=-1,...),hessian=TRUE))
}


##' simulateProcess function
##'
##' A function to simulate a realisation of a process, only works with one individual grid of data ATM
##'
##' @param param vector of back transformed parameters
##' @param startcell index of cell with first infection
##' @param nmax number of cells to be infected
##' @param nrow number of rows in grid
##' @param ncol number of columns in grid
##' @param iForm formula for infectivity
##' @param sForm formula for susceptibility
##' @param dat a HCVdataset object, see ?asHCVdata, this will contain the covariate data for simulation
##' @return an output grid of HCV values: larger values represent older infections
##' @export

simulateProcess <- function(param,startcell,nmax,nrow=10,ncol=10,iForm=NULL,sForm=NULL,dat=NULL){

    # if(!inherits(dat,"HCVdataset")){
    #     stop("dat must inherit class 'HCVdataset', see ?asHCVdata")
    # }

    theta <- param[1]
    phi <- param[2]
    beta <- param[-c(1,2)]

    sus <- 1
    inf <- 1
    N <- c(0,0)
    if(!is.null(iForm)){
        inf <- model.matrix(iForm,dat)
        inf <- inf[,-1,drop=FALSE] # remove intercept
        N[1] <- ncol(inf)
    }
    if(!is.null(sForm)){
        sus <- model.matrix(sForm,dat)
        sus <- sus[,-1,drop=FALSE] # remove intercept
        N[2] <- ncol(sus)
    }

    gr <- matrix(0,nrow,ncol)
    df <- distFun(gr)
    inftime <- nmax
    gr[startcell] <- inftime
    inftime <- inftime - 1

    S <- (1:100)[-startcell]
    I <- startcell

    count <- 2

    while(inftime>=1){
        rho_mats <- df(I,S,theta,phi,beta,inf,sus,N)
        den <- sum(rho_mats)
        num <- colSums(rho_mats)

        if(length(S)>1){
            infidx <- sample(S,1,prob=num/den)
        }
        else{
            infidx <- S
        }

        S <- S[-which(S==infidx)]
        I <- c(I,infidx)
        gr[infidx] <- inftime

        count <- count + 1
        inftime <- inftime - 1
    }
    return(gr)
}



##' bootCI function
##'
##' A function to compute boostrapped confidence intervals for the partial likelihood model
##'
##' @param param backtransformed parameters
##' @param nboot number of bootstrap samples
##' @param dat a HCVdataset object, see ?asHCVdata
##' @param iForm formula for infectivity
##' @param sForm formula for susceptibility
##' @param nr integer number of rows in grid
##' @param nc integer number of columns in grid
##' @return matrix containing bootstrapped parameter estimates
##' @export

bootCI <- function(param,nboot,dat,iForm=NULL,sForm=NULL,nr=NULL,nc=NULL){

    pd <- prepData(dat=dat,iForm=iForm,sForm=sForm)

    #browser()

    # sid <- sort(unique(dat$subjectid))
    # nsub <- length(sid)
    # gid <- list()
    # for(i in 1:length(sid)){
    #     gid[[i]] <- sort(unique(dat$gridnum[dat$subject==sid[i]]))
    # }
    # subidx <- rep(si)

    nmats <- length(pd$matList)

    pb <- txtProgressBar(1,nboot,1)
    start <- Sys.time()
    pars <- matrix(NA,nboot,length(param))
    for(i in 1:nboot){
        idx <- sample(1:nmats,1)
        #print(idx)
        startcell <- which(pd$matList[[idx]]==max(pd$matList[[idx]]))
        nmax <- sum(pd$matList[[idx]]>0)
        if(is.null(nr)){
            nr <- as.integer(sqrt(length(pd$matList[[idx]])))
        }
        if(is.null(nc)){
            nc <- as.integer(sqrt(length(pd$matList[[idx]])))
        }

        if(is.null(iForm)){
            pd$covList_inf[[idx]] <- t(t(rep(1,nr*nc)))
            names(pd$covList_inf[[idx]]) <- "dummy_inf"
        }
        if(is.null(sForm)){
            pd$covList_sus[[idx]] <- t(t(rep(1,nr*nc)))
            names(pd$covList_sus[[idx]]) <- "dummy_sus"
        }
        #browser()
        datss <- as.data.frame(cbind(pd$covList_inf[[idx]],pd$covList_sus[[idx]]))

        sp <- simulateProcess(param,startcell,nmax,nrow=nr,ncol=nc,iForm=iForm,sForm=sForm,dat=datss)
        datss$hcv <- as.numeric(sp)
        datss$subjectid <- 1
        datss$gridnum <- 1
        grd <- expand.grid(1:nr,1:nc)
        datss$x <- grd[,1]
        datss$y <- grd[,2]
        datss <- asHCVdata(dat=datss,hcv="hcv",subjectid="subjectid",gridnum="gridnum",x="x",y="y")

        p <- optimPL(dat=datss,iForm=iForm,sForm=sForm)$par
        pars[i,] <- p

        setTxtProgressBar(pb,i)
    }
    end <- Sys.time()
    print(end-start)
    close(pb)
    return(pars)
}


# not needed:
# ij2idx <- function(i,j,M){
#     return(i + (j-1)*M)
# }
