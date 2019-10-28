library(HCV)

library(reshape2)
library(ggplot2)
library(viridis)
library(grid)
library(gridExtra)

# source("code/1fncalculatetimeinfected.R")
# source("code/1fnpartiallikelihood.R")
# source("code/3simulatepatternsandexecutepartiallikelihood.R")




# code is read from ruta (does not have typedata).
# Data is read using typedata. Plots and tables are saved using typedata
ruta <- paste0(ruta, "/", typedata)


# Read real data
if(typedata == ""){
    dwhole <- read.csv(paste0(ruta,"/data/HCV-mono-average_cells.csv"))
    dwhole$logifitlevel <- log(dwhole$ifit.level)
    dwhole$roundlogifitlevel <- round(log(dwhole$ifit.level), 1)
}


# Read synthetic data
if(typedata == "synthetic"){

    # Load synthetic data. There are 5 10x10 biopsies of the same "patient" and one 100x100 "biopsy" generated in silico.
    #"I have asked Ashish Goyal, the postdoc at LANL, to generate some synthetic data for Paula to analyze along the lines of what she had done
    #before with the actual biopsy data.
    #This was something we had thought could be a good next step. These data were generated with a stochastic simulation for a specific set of parameters.
    #I am sending 5 data files, corresponding to 5 10x10 biopsies of the same "patient" (similar to the actual data) and then for comparison,
    #one 100x100 "biopsy" generated in silico. The information in the files is the same as usual, (x,y) coordinates and HCV RNA level."

    # Put 5 synthetic data in the same format as real data. With the same columns


    dsynwhole <- NULL
    for(numsyn in 1:5){
        dsyn <- read.csv(paste0(ruta,"/data/Dataset_", numsyn, ".csv"))

        # I set subid equal to 20000. In this way, when I do sort(unique(dwhole$subid)) 20000 is the last one
        dsynOneDataset <- data.frame(subid = 20000,
                           gridN = numsyn,
                           xpos = dsyn$X,
                           ypos = dsyn$Y,
                           hcv.level = dsyn$HCV.RNA..copies.mL.,
                           il28b.gt = NA,
                           ifit.level = NA)

        dsynwhole <- rbind(dsynwhole, dsynOneDataset)
    }
    dwhole <- dsynwhole
}


print("finish 0librariesandpath.R")
head(dwhole)

# Generate sample dataset of each sample j of individual i with the parameters obtained by mle
# It calls fnGenerateSample
# If indorall=="all" the mle are the ones obtained pooling the samples of individual i
# If indorall=="ind" the mle are the ones obtained for individual i and sample j
##' fnGenerateSampleDataset function
##'
##' A function to
##'
##' @param i X
##' @param j X
##' @param indorall X
##' @return ...
##' @export
fnGenerateSampleDataset<-function(i,j,indorall){
    # 1. read real dataset
    # Choose individual i and sample j (they are selected above)
    d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]
    d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
    d$numcell<-1:100
    # num cells infected in the sample
    (n<-length(which(d$ranktimeinfected!=0)))
    # 2. read maximum likelihood estimates of the real data (replicate0)
    if(indorall=="ind"){
        nameguardar<-paste0("model",1,"ind",i,"sample",j,"replicate0")
    }
    if(indorall=="all"){
        nameguardar<-paste0("model",1,"ind",i,"allsamples","replicate0")
    }
    params<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
    #params in log scale
    #borrar mlephi<-exp(p)[1]
    #borrar mlerho<-exp(p)[2]
    #borrar (params<-c(mlephi, mlerho))
    # 3. Generate sample with the parameters (conditional on time n)
    d<-fnGenerateSample(d, n, params)
    return(d)
}
