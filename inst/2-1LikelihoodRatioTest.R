#optim performs minimization
#PL returns likelihood multiplied by (-1) to maximize (optim does minimization)
#valPL<-PL(c(1,3), d)
#When I get partial likelihood from function PL I multiply by (-1) to get the PL


############################################################################
# Get mles when I have only one sample (not pooled data)
#H_1: phi!=0. max likelihood if phi!=0
nameguardarphinot0<-paste0(nameguardar,"phinot0.dat")
res<-optim(c(1,1), PL, d=d, hessian=TRUE)
save(res,      file=paste0(ruta,"/res/object",    nameguardarphinot0,".dat"))
write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardarphinot0,".txt"))


#H_0: phi=0. max likelihood if phi=0
nameguardarphi0<-paste0(nameguardar,"phi0.dat")
#one-dimensional optimization (use method Brent)
res<-optim(c(1), PL, d=d, hessian=TRUE, "Brent", method="Brent", lower=log(1e-100), upper=log(1e+10))
save(res,      file=paste0(ruta,"/res/object",    nameguardarphi0,".dat"))
write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardarphi0,".txt"))
############################################################################

##############################
##############################
##############################
##############################

#real data
vectaskid <- 1:15

if(typedata == "synthetic"){
    vectaskid <- 17:22
}

matrizresultado<-NULL
for(taskid in vectaskid){
    print(taskid)


    #Convergence equal to 0 indicates successful completion of the function, and 1 indicates that the iteration limit maxit had been reached.

    #Convergence
    #N. parametric bootstrap simulations did not converge.

    # Returns number of parametric bootstrap for which optim function to calculate mles did not converge
    fnCalculateHowManySimulationParametricBootstrapDoNotConverge<-function(nameguardar){
        # For the moment I have not done parametric bootstrap pooling all data of all ind and all samples
        if(nameguardar=="model1allindallsamplesreplicate0"){
            return("parametric bootstrap not done")
        }
        else{
            vecresconvergence<-NULL
            for(repl in 1:100){
                sub("replicate0",paste0("replicate",repl),nameguardar)
                load(file=paste0(ruta,"/res/object", nameguardar,".dat"))
                vecresconvergence<-c(vecresconvergence,res$convergence)
            }
        }
        return(length(which(vecresconvergence!=0)))
    }





    # Calculate LRT statistic and pvalue
    fnCalculateLRS<-function(loglikH1,loglikH0,resconvergence){
        D <- 2*(loglikH1-loglikH0)
        pvalor <- 1-pchisq(D,1)
        signif<-"Not signif"
        if(pvalor<0.05){signif<-"Signif"}

        #res$convergence is an integer code.
        #0 indicates successful completion. 1 indicates that the iteration limit maxit had been reached.
        if(resconvergence==0){resconvergence<-"Converged"}
        if(resconvergence==1){resconvergence<-"Iteration limit reached"}

        #resultado<-paste("individual", i, "sample", j,"loglikH1:", round(loglikH1,2), "loglikH0:", round(loglikH0,2),
        #                 "D:", round(D,2), "pvalor", round(pvalor,2), signif,resconvergence)

        resultado<-paste("&", i, "&", j,"&", round(loglikH1,2), "&", round(loglikH0,2),
                        "&", round(D,2), "&", round(pvalor,2), "&", signif, "&", resconvergence)

        return(resultado)
    }



    if(taskid<=11 | (17<= taskid & taskid<=21)){

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
            #sample
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



        nameguardar<-paste0("model",1,"ind",i,"sample",j,"replicate0")

        #############################################################
        #############################################################
        #Test H_0: phi=0 vs. phi!=0
        #PL returns likelihood multiplied by (-1) to maximize (optim does minimization)
        # hat is why a put negative sign in the loglikH1
        #H_0. max likelihood if phi!=0
        print(nameguardar)
        load(file=paste0(ruta,"/res/object", nameguardar,".dat"))
        loglikH1<- - res$value
        #H_0. max likelihood if phi==0
        #The value of the PL for any value of rho is the same
        loglikH0<- - PL(c(1), d)
        resultado<-fnCalculateLRS(loglikH1,loglikH0,res$convergence)
        numberbootstrapmlesdidnotconverge<-fnCalculateHowManySimulationParametricBootstrapDoNotConverge(nameguardar)
        #resultado<-paste(resultado, "numberbootstrapdidnotconverge",numberbootstrapmlesdidnotconverge)
        resultado<-paste(resultado, "&" ,numberbootstrapmlesdidnotconverge, "\\\\")
        print(resultado)
        #############################################################
        #############################################################


    }

    ###############################################



    if( (taskid>=12 & taskid<=15) | taskid == 22){
        #############################
        # INI Pool data with all repeated samples of individual i
        #############################

        #real data
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


        nameguardar<-paste0("model",1,"ind",i,"allsamples","replicate0")
        j<-"all"

        #############################################################
        #############################################################
        #Test H_0: phi=0 vs. phi!=0
        #PL returns likelihood multiplied by (-1) to maximize (optim does minimization)
        # hat is why a put negative sign in the loglikH1
        #H_0. max likelihood if phi!=0
        print(nameguardar)
        load(file=paste0(ruta,"/res/object", nameguardar,".dat"))
        loglikH1<- - res$value
        #H_0. max likelihood if phi==0
        #The value of the PL for any value of rho is the same
        loglikH0<- - PLAllSamples(c(1), dlist)
        resultado<-fnCalculateLRS(loglikH1,loglikH0,res$convergence)
        numberbootstrapmlesdidnotconverge<-fnCalculateHowManySimulationParametricBootstrapDoNotConverge(nameguardar)
        #resultado<-paste(resultado, "Number bootstrap simulations did not converge",numberbootstrapmlesdidnotconverge)
        resultado<-paste(resultado, "&" ,numberbootstrapmlesdidnotconverge, "\\\\")
        print(resultado)
        #############################################################
        #############################################################


    }



    if(taskid==16){
        #############################
        # INI Pool data with all repeated samples of individual i
        #############################

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




        nameguardar<-paste0("model",1,"allind","allsamples","replicate0")
        i<-"all"
        j<-"all"

        #############################################################
        #############################################################
        #Test H_0: phi=0 vs. phi!=0
        #PL returns likelihood multiplied by (-1) to maximize (optim does minimization)
        # hat is why a put negative sign in the loglikH1
        #H_0. max likelihood if phi!=0
        print(nameguardar)
        load(file=paste0(ruta,"/res/object", nameguardar,".dat"))

        # No he hecho el 16. Se hace en 2executepartiallikelihood.R Eso hay q borrarlo cuando lo haga
        ### I get mles here different from the other data
        #res<-optim(c(.000000001,.1), PLAllSamples, d=dlist, hessian=TRUE)
        #save(res,      file=paste0(ruta,"/res/object",    nameguardar,".dat"))
        #write(res$par, file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
        ###

        loglikH1<- - res$value
        #H_0. max likelihood if phi==0
        #The value of the PL for any value of rho is the same
        loglikH0<- - PLAllSamples(c(1), dlist)
        resultado<-fnCalculateLRS(loglikH1,loglikH0,res$convergence)
        # For the moment I do not parametric bootstrap pooling all samples and all individuals
        #numberbootstrapmlesdidnotconverge<-fnCalculateHowManySimulationParametricBootstrapDoNotConverge(nameguardar)
        numberbootstrapmlesdidnotconverge<-"not done"
        #resultado<-paste(resultado, "Number bootstrap simulations did not converge",numberbootstrapmlesdidnotconverge)
        resultado<-paste(resultado, "&" ,numberbootstrapmlesdidnotconverge, "\\\\")
        print(resultado)
        #############################################################
        #############################################################


    }

    matrizresultado<-rbind(matrizresultado,paste(taskid,resultado))
    matrizresultado
}

matrizresultado


tablalatexresultado<-paste(
#"\\begin{tabular}{c ccccccccc}
#taskid&individual&sample&loglikH1&loglikH0&D&p-value&Significance&Convergence&N. bootstrap samples did not converge \\\\ \\hline",
paste(matrizresultado,collapse="")
#,"\\end{tabular}"
)


write(tablalatexresultado, file=paste0(ruta,"/res/tablewithLRTresults", "model1",".txt"))
