# Ejecutar este
########################################
#qsub -l h_vmem=4G R_job17.com
#qsub -l h_vmem=4G R_job.com
#Results are in HCV/synthetic/res
########################################
ruta<-"/home/users/moragasp/HCV"
#ruta <- "C:/Users/moragasp/HCV"
setwd(ruta)

########################################
# CHOOSE
# real data
typedata <- ""
# synthetic data
typedata <- "synthetic"

# without covariates
modelnum <- 1
# with covariates to model alpha and beta
modelnum <- 2
# I have changed 2executepartiallikelihood.R and 1partiallikelihood.R:
# name model is modelnum, optim() accepts 3 initial values, partial likelihood formula is changed with covariates
# need to change 4plotTransmissitivityfunction.R, 2-1LikelihoodRatioTest.R
########################################

source(system.file("inst","0librariesandpath.R",package="HCV"))

#source("code/0librariesandpath.R")


#taskid goes from 1 to 11. 12,13,14,15 will do all samples each individual. 16 will do all samples all individuals,
#taksid 17-21 (each samples synthetic individual), 22 (all samples synthetic individual)
taskid<-as.numeric(Sys.getenv("SGE_TASK_ID"))

#taskid<-17

###############################################
# each sample of individuals
if(taskid<=11 | (17<= taskid & taskid<=21)){
    print(taskid)

    #####################
    # INI MAXIMIZE LIKELIHOOD
    #####################
    # Analysis with individual i and sample j

    # real data
    if(taskid<=11){
        #individual
        i<-ceiling(taskid/3)
        #sample
        j<-taskid-(3*(i-1))
    }

    # synthetic data
    if(17<= taskid & taskid<=21){
        #individual
        i<-1
        #sample 1 to 5
        j<-taskid-16
    }

    d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]
    d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
    #From 1 until 100 (they are in ordered!)
    d$numcell<-1:100
    #Numbers are 1:10, first column from bottom to top
    #Numbers are 11:20, second column from bottom to top

    #optim performs minimization
    #PL return multiplied by (-1) to maximize
    #valPL<-PL(c(1,3), d)



    nameguardar<-paste0("model",modelnum, "ind",i,"sample",j,"replicate0")
    initialvalues <- c(1,1)
    if(modelnum == 2){
        initialvalues <- c(1,1,1)
    }
    res<-optim(initialvalues, PL, d=d, hessian=TRUE)
    #control=list(fnscale=-1)
    save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
    write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))

    summary(res)
    res$hessian
    res

    #####################
    # END MAXIMIZE LIKELIHOOD
    #####################



    #####################
    # INI Parametric Bootstrap. Simulate point pattern using maximum likelihood estimates
    #####################


    # Generate sample with infection times (Generate R samples)
    R<-100
    for(r in 1:R){
        print(paste("Replicate", r))
        d<-fnGenerateSampleDataset(i,j,"ind")

        # INI BORRAR
        ## Choose individual i and sample j (they are selected above)
        #d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]
        #d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
        #d$numcell<-1:100
        ## num cells infected in the sample
        #(n<-length(which(d$ranktimeinfected!=0)))

        ## maximum likelihood estimates
        #nameguardar<-paste0("model",1,"ind",i,"sample",j,"replicate0")
        #p<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
        #mlephi<-exp(p)[1]
        #mlerho<-exp(p)[2]
        #(params<-c(mlephi, mlerho))
        #d<-fnGenerateSample(d, n, params)
        # END BORRAR

        # Estimate parameters in the simulated pattern
        nameguardar<-paste0("model",modelnum,"ind",i,"sample",j,"replicate",r)
        initialvalues <- c(1,1)
        if(modelnum == 2){
            initialvalues <- c(1,1,1)
        }
        res<-optim(initialvalues, PL, d=d, hessian=TRUE)
        save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
        write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
    }



    #####################
    # END Parametric Bootstrap. Simulate point pattern using maximum likelihood estimates
    #####################

}

###############################################


if( (taskid>=12 & taskid<=15) | taskid == 22){
    #############################
    # INI Pool data with all repeated samples of individual i
    #############################

    # real data
    numreplic<-3
    if(taskid==12){i<-1}
    if(taskid==13){i<-2}
    if(taskid==14){i<-3}
    if(taskid==15){
        i<-4
        numreplic<-2
    }

    # synthetic data
    if(taskid==22){
        i<-1
        numreplic<-5
    }

    dlist<-list()
    for(index in 1:numreplic){
        j<-index
        d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]
        d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
        d$numcell<-1:100
        dlist[[j]]<-d
    }

    print("start maximize PLAllSamples replicate0")
    nameguardar<-paste0("model",modelnum,"ind",i,"allsamples","replicate0")
    initialvalues <- c(1,1)
    if(modelnum == 2){
        initialvalues <- c(1,1,1)
    }
    res<-optim(initialvalues, PLAllSamples, d=dlist, hessian=TRUE)
    save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
    write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
    #############################
    # END Pool data with all repeated samples of individual i
    #############################

    #####################
    # INI Parametric Bootstrap
    #####################

    # Generate sample with infection times (Generate R samples)
    R<-100
    for(r in 1:R){
        print(paste("Replicate", r))
        #d<-fnGenerateSample(d, n, params)
        dlist<-list()
        for(index in 1:numreplic){
            j<-index
            d<-fnGenerateSampleDataset(i,j,"all")
            dlist[[j]]<-d
        }
        ################
        # Estimate parameters in the simulated pattern
        nameguardar<-paste0("model",modelnum,"ind",i,"allsamples","replicate",r)
        initialvalues <- c(1,1)
        if(modelnum == 2){
            initialvalues <- c(1,1,1)
        }
        res<-optim(initialvalues, PLAllSamples, d=dlist, hessian=TRUE)
        save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
        write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
    }




    #####################
    # END Parametric Bootstrap
    #####################


}

# This is using all samples of all individuals together. Do not do it

if(taskid==16){

    dlist<-list()
    #indicelist es para la posicion de i y j
    indicelist<-0
    numreplic<-3
    for(i in 1:4){
        if(i==4){numreplic<-2}
        for(index in 1:numreplic){
            indicelist<-indicelist+1
            j<-index
            d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]
            d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
            d$numcell<-1:100
            dlist[[indicelist]]<-d
        }
    }



    print("start maximize PLAllSamples (all samples all ind) replicate0")
    nameguardar<-paste0("model",modelnum,"allind","allsamples","replicate0")
    initialvalues <- c(.0001,.1)
    if(modelnum == 2){
        initialvalues <- c(.0001,.1,1)
    }
    res<-optim(initialvalues, PLAllSamples, d=dlist, hessian=TRUE)
    save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
    write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))


    # For the moment I do not do Parametric Bootstrap


}
