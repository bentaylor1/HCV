#####################
# INI CONSTRUCT 95% CONFIDENCE INTERVAL AND DO PLOT TRANSMISSITIVITY FUNCTION
#####################
#initial phi (decay) = 1
#initial rho (spontaneous) = 3
#p<-c(-2.248119, -14.037492)
#hessian<-matrix(c(315.65513, -33.505716, -33.50572, 3.563991), nrow=2)

#vecnumsamplesi contains the number of samples of each individual: 3,3,3and2.
#The last element is the number of individuals in the dataset.
#We need to know that when we pool all samples of all individuals. for will be i in 1:5
#For the moment we only maximize pooled PL for the real data. We do not do parametric bootstrap.
#We plot transmissivity function using weighted averages. We do for(i in 1:5){ down.


numindividuals<-4
vecnumsamplesi<-c(3,3,3,2,4)
if(typedata == "synthetic"){
    vecnumsamplesi<-c(5,1)
    numindividuals<-1
}
for(i in 1:numindividuals){
    numsamplesi<-vecnumsamplesi[i]

    for(j in 0:numsamplesi){
        # Save table with parameters
        vectablelatexparamphirho<-NULL

        print(paste("individual",i))
        print(paste("sample",j))
        #when j=0 I maximize partial likelihood pooling all samples of individual i


        #i<-4
        #j<-1
        r<-0

        ###############################

        nameguardarplot<-paste0("model",1,"ind",i,"sample",j)
        # read mle estimates for the real data (replicate r=0)
        nameguardar<-    paste0("model",1,"ind",i,"sample",j,"replicate",r)
        if(j==0){
            nameguardarplot<-paste0("model",1,"ind",i,"allsamples")
            nameguardar<-    paste0("model",1,"ind",i,"allsamples","replicate",r)
        }

        #I need load object if I calculate 95% CI with the hessian
        #load(file=paste0(ruta,"/res/object", nameguardar,".dat"))
        ##p<-res$par
        p<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))


        compute95CIusinghessian<-FALSE

        #################################################
        #plot gamaij for the real data (replicate 0)
        #dat contains 2 columns: distance and gammaij for each distance
        #data.frame estimates has the parameters mle and 95%CI
        #95% CI calculated using hessian
        if(compute95CIusinghessian==TRUE){

            hessian<-res$hessian
            # 1. The Hessian is the matrix of second derivatives of the likelihood with respect to the parameters
            # 2. The Information matrix is the negative of the expected value of the Hessian matrix
            #    Si fn devuelve -likelihood no tengo q multiplicar por (-1) el hessian.
            information_matrix <- hessian #-hessian
            # 3. The variance of an ML estimator is the inverse of the Information matrix:
            fisher_info<-solve(information_matrix)

            # 4. The standard errors of the estimator are the square roots of the diagonal terms in
            #    the variance-covariance matrix.
            prop_sigma<-sqrt(diag(fisher_info))

            upper<-p+1.96*prop_sigma
            lower<-p-1.96*prop_sigma
            interval<-data.frame(value=p, lower=lower, upper=upper)
            estimates<-exp(interval)
            phi<-estimates[1,1]
            rho<-estimates[2,1]

            dist<-seq(0,2,length.out=100)
            gamma<- rho+exp(-dist/phi)
            dat<-data.frame(dist=dist, gamma=gamma)

        }

        #################################################
        #plot gamaij for the real data (replicate 0) and the other replicates r=1,2,...
        #dat contains 2+r columns: distance and gammaij for each distance for real data and r=1,2,...
        #data.frame estimates has the parameters mle and 95%CI
        #95% CI calculated as quantiles of the mle of the parameters for the replicates
        if(compute95CIusinghessian==FALSE){
            dist<-seq(0,10,length.out=100)
            #gamma<- rho+exp(-dist/phi)
            dat<-data.frame(dist=dist)
            matestimates<-NULL



            # Real data
            if(i==1 & j==0){R<-93}
            if(i==1 & j==1){R<-93}
            if(i==1 & j==2){R<-93}
            if(i==1 & j==3){R<-93}

            if(i==2 & j==0){R<-93}
            if(i==2 & j==1){R<-93}
            if(i==2 & j==2){R<-93}
            if(i==2 & j==3){R<-93}

            if(i==3 & j==0){R<-93}
            if(i==3 & j==1){R<-93}
            if(i==3 & j==2){R<-93}
            if(i==3 & j==3){R<-93}

            if(i==4 & j==0){R<-93}
            if(i==4 & j==1){R<-93}
            if(i==4 & j==2){R<-93}

            # Synthetic data
            if(typedata == "synthetic"){
                if(i==1 & j==0){R<-100}
                if(i==1 & j==1){R<-100}
                if(i==1 & j==2){R<-100}
                if(i==1 & j==3){R<-100}
                if(i==1 & j==4){R<-100}
                if(i==1 & j==5){R<-100}
            }

            for(r in 0:R){
                print(r)
                nameguardar<-paste0("model",1,"ind",i,"sample",j,"replicate",r)
                if(j==0){
                    nameguardar<-paste0("model",1,"ind",i,"allsamples","replicate",r)
                }
                p<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))

                interval<-data.frame(value=p)
                estimates<-exp(interval)
                phi<-estimates[1,1]
                rho<-estimates[2,1]

                gamma<- rho+exp(-dist/phi)
                dat<-cbind(dat,gamma)

                matestimates<-rbind(matestimates,t(estimates))

            }
            #no quito real data (first row)
            #borrar estimates<-cbind(matestimates[1,],t(apply(matestimates[-1,],2,quantile,c(0.025,0.975))))
            estimates<-cbind(matestimates[1,],t(apply(matestimates,2,quantile,c(0.025,0.975))))
        }



        ##############################
        # INI to calculate pooled mle, save variance of the parameter mles and gammas obtained in each bootstrap sample
        ##############################
        # calculate estimates of the pooled data (all samples of an individual)
        # as a weighted average of the parameters of each sample with weights proportional to the variance
        # that I obtianed from the bootstrap samples

        #save variance parameters (phi, rho) for the ind i sample j obtained from bootstrap samples
        #first row of matestimates are the parameters mle obtained in real data
        #I keep the row of the real data
        meanparambootstrap<-apply(matestimates,2,mean) #instead of this use the mle in the real data?
        varparambootstrap<-apply(matestimates,2,var)

        #save variance gamma(dist) for each distance for the ind i sample j obtained from bootstrap samples
        #column 1 is distance, column 2 is gamma(dist) for the real data,
        #the rest of columns are gamma(dist) for the simulations 1,...,R
        #I get variance for each of the distances in vector dist
        meangammabootstrap<-apply(dat[,-c(1)],1,mean)
        vargammabootstrap<-apply(dat[,-c(1)],1,var)


        write(meanparambootstrap, file=paste0(ruta,"/res/meanparambootstrap", nameguardarplot,".txt"))
        write(varparambootstrap,  file=paste0(ruta,"/res/varparambootstrap", nameguardarplot,".txt"))
        write(meangammabootstrap, file=paste0(ruta,"/res/meangammabootstrap", nameguardarplot,".txt"))
        write(vargammabootstrap,  file=paste0(ruta,"/res/vargammabootstrap", nameguardarplot,".txt"))



        ##############################
        # END to calculate pooled mle, save variance of the parameter mles and gammas obtained in each bootstrap sample
        ##############################



        #################################################
        # INI PLOT

        #apply melt to be able to plot all lines in the same plot
        names(dat)=c("dist", 1:(ncol(dat)-1))
        datmelted<-melt(dat, id.vars="dist")

        # Plot transmissivity function gamma_{ij}
        titulo<-bquote(paste(ind," ",.(i) ," ",sample," ",.(j)," ",
            phi==.(signif(estimates[1,1],3)),
            " (", .(signif(estimates[1,2],3)),",", .(signif(estimates[1,3],3)),"), ",
            rho==.(signif(estimates[2,1],3)),
            " (", .(signif(estimates[2,2],3)),",", .(signif(estimates[2,3],3)),") "))
        gp<-ggplot(data=datmelted, aes(x=dist, y=value, group=variable))  +
            geom_line(col="gray")+
            geom_line(data=datmelted[which(datmelted$variable==1),],aes(x=dist, y=value, group=variable),size=1.2)+
            theme_bw() +  ggtitle(titulo) +# theme(text=element_text(size = 8)) +
            scale_x_continuous(breaks= 0:10, labels= 0:10)+#, limits)+
            xlab("distance") + ylab(expression(paste("transmissivity   ",rho+exp(-distance/phi))))
        #+  coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept=rho, col="red") +
        #annotate("text", x = 0.1, y = .05, label = paste("rho == ", signif(rho,3)),parse=TRUE)
        png(paste0(ruta,"/plots/transmissitivity",nameguardarplot,".png"))
        plot(gp)
        dev.off()

        # END PLOT
        #################################################

        #####################
        # END CONSTRUCT 95% CONFIDENCE INTERVAL AND DO PLOT TRANSMISSITIVITY FUNCTION
        #####################
        #http://www.cookbook-r.com/Graphs/Bar_and_line_graphs_(ggplot2)/

        # No escribo j. Solo voy a utilizar cuando j=0 (all samples)
        vectablelatexparamphirho<-c(vectablelatexparamphirho, i, "&",#, j, "&",
        paste(signif(estimates[1,1],3),
          " (",signif(estimates[1,2],3),",", signif(estimates[1,3],3),") &",
          signif(estimates[2,1],3),
          " (", signif(estimates[2,2],3),",", signif(estimates[2,3],3),")\\\\"))
        write(vectablelatexparamphirho,paste0(ruta,"/res/tablelatexparamphirho",nameguardarplot,".txt"))


        ####
        #end loops i and j
    }
}


######################


##############################################
# DE AQUI PARA ARRIBA HAY CAMBIOS HECHOS PARA SYNTHETIC DATA

# DE AQUI PARA ABAJO TODAVIA NO LO HE CAMBIADO
##############################################



############################################################################
# INI Calculate estimates of the pooled datasets using weighted average
# (all samples in one individual, or all samples for all individuals)
# by doing the weighted average of the estimates at each sample
# with weights proportional to the variance obtained in the bootstrap samples

#If I want to pool all samples do loop only for individual i and use this.
#if(j==0){
#  nameguardar<-paste0("model",1,"ind",i,"allsamples","replicate",r)
#}

# i = 1-4 pooling all samples of each individual, i = 5 all samples all individuals
for(i in 1:5){
    numsamplesi<-vecnumsamplesi[i]


    vecestparam1<-NULL
    vecestparam2<-NULL
    matestgamma<-NULL
    vecweightsparam1<-NULL
    vecweightsparam2<-NULL
    matweightsgamma<-NULL
    nameguardarweighted<-paste0("model",1,"ind",i,"allsamples","weighted")
    nameguardarallsim<-paste0("model",1,"ind",i,"allsamples","allsim")
    if(i==5){
        nameguardarweighted<-paste0("model",1,"allind","allsamples","weighted")
        nameguardarallsim<-paste0("model",1,"allind","allsamples","allsim")
    }

    #for(j in 0:numsamplesi){  # j=0 is for all samples
    #I start with 0 because I put real data also
    for(j in 0:numsamplesi){
        print(paste("individual",i))
        print(paste("sample",j))

        # load parameter mles for the real data for individual i and sample j
        r<-0 # real data
        nameguardar<-paste0("model",1,"ind",i,"sample",j,"replicate",r) #r=0 for getting the estimates
        nameguardarplot<-paste0("model",1,"ind",i,"sample",j) #for getting the variance
        if(i==5){
            #here change i by j
            nameguardar<-paste0("model",1,"ind",j,"allsamples","replicate",r) #r=0 for getting the estimates
            nameguardarplot<-paste0("model",1,"ind",j,"allsamples") #for getting the variance
        }

        p<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
        interval<-data.frame(value=p)
        estimates<-exp(interval)
        phi<-estimates[1,1]
        rho<-estimates[2,1]
        gamma<- rho+exp(-dist/phi)


        #estimates
        estparam<-c(phi,rho)
        estgamma<-gamma
        vecestparam1<-c(vecestparam1,estparam[1])
        vecestparam2<-c(vecestparam2,estparam[2])
        matestgamma<-rbind(matestgamma,estgamma)

        # weights
        varparam<-scan(file=paste0(ruta,"/res/varparambootstrap", nameguardarplot,".txt"))
        vargamma<-scan(file=paste0(ruta,"/res/vargammabootstrap", nameguardarplot,".txt"))
        vecweightsparam1<-c(vecweightsparam1,varparam[1])
        vecweightsparam2<-c(vecweightsparam2,varparam[2])
        matweightsgamma<-rbind(matweightsgamma,vargamma)


    }

    param1weighted<-sum(vecestparam1*1/vecweightsparam1)/sum(1/vecweightsparam1)
    param2weighted<-sum(vecestparam2*1/vecweightsparam2)/sum(1/vecweightsparam2)
    # multiply column of one matrix by column of other matrix, divide by sum columns
    gammaweighted<-colSums(matestgamma*1/matweightsgamma)/colSums(1/matweightsgamma)

    # mean and CI putting together all parameters obtained in all simulations
    param1meanallsim<-mean(vecestparam1)
    param1quant1allsim<-quantile(vecestparam1,0.025)
    param1quant2allsim<-quantile(vecestparam1,0.975)

    param2meanallsim<-mean(vecestparam2)
    param2quant1allsim<-quantile(vecestparam2,0.025)
    param2quant2allsim<-quantile(vecestparam2,0.975)



    # save
    write(param1weighted,  file=paste0(ruta,"/res/param1weighted", nameguardarweighted,".txt"))
    write(param2weighted,  file=paste0(ruta,"/res/param2weighted", nameguardarweighted,".txt"))
    write(gammaweighted,  file=paste0(ruta,"/res/gammaweighted", nameguardarweighted,".txt"))
    #save mean and quantiles all simulations
    write(param1meanallsim,    file=paste0(ruta,"/res/param1meanallsim", nameguardarallsim,".txt"))
    write(param1quant1allsim,  file=paste0(ruta,"/res/param1quant1allsim", nameguardarallsim,".txt"))
    write(param1quant2allsim,  file=paste0(ruta,"/res/param1quant2allsim", nameguardarallsim,".txt"))
    write(param2meanallsim,    file=paste0(ruta,"/res/param2meanallsim", nameguardarallsim,".txt"))
    write(param2quant1allsim,  file=paste0(ruta,"/res/param2quant1allsim", nameguardarallsim,".txt"))
    write(param2quant2allsim,  file=paste0(ruta,"/res/param2quant2allsim", nameguardarallsim,".txt"))
}




# END Calculate estimates of the pooled datasets using weighted average
############################################################################







############################################################################
# INI PLOT Compare estimates obtained with bootstrap sampling and estimates obtained by weighted average

# i = 1-4 pooling all samples of each individual, i = 5 all samples all individuals
for(i in 1:5){
    print(i)
    j<-0
    nameguardar<-paste0("model",1,"ind",i,"allsamples","replicate0")
    nameguardarweighted<-paste0("model",1,"ind",i,"allsamples","weighted")
    nameguardarallsim<-paste0("model",1,"ind",i,"allsamples","allsim")
    if(i==5){
        nameguardar<-paste0("model",1,"allind","allsamples","replicate0")
        nameguardarweighted<-paste0("model",1,"allind","allsamples","weighted")
        nameguardarallsim<-paste0("model",1,"allind","allsamples","allsim")
    }
    nameguardarplot<-nameguardarweighted

    # 1. transimissivity obtained using the parameter estimates obtained by maximizing the PL of the pooled data
    # For the moment I have not fitted all samples all ind. I put this gamma=0
    if(i!=5){
        p<-scan(file=paste0(ruta,"/res/paramsmle", nameguardar,".txt"))
        interval<-data.frame(value=p)
        estimates<-exp(interval)
        phi<-estimates[1,1]
        rho<-estimates[2,1]
        gamma<- rho+exp(-dist/phi)
    }
    else{
        estimates<-""
        phi<-""
        rho<-""
        gamma<-dist*0
    }


    # 3. weighted at each distance
    param1weighted<-scan(file=paste0(ruta,"/res/param1weighted", nameguardarweighted,".txt"))
    param2weighted<-scan(file=paste0(ruta,"/res/param2weighted", nameguardarweighted,".txt"))
    gammaweightedateachdistance<-scan(file=paste0(ruta,"/res/gammaweighted", nameguardarweighted,".txt"))


    # 2. transimissivity obtained using the parameter estimates obtained by averaging the parameter mles in each of the samples
    gammapluginaveragemles<- param2weighted+exp(-dist/param1weighted)
    print(c(param2weighted,param1weighted))

    # 4. pooling all simulations
    param1meanallsim<-scan(file=paste0(ruta,"/res/param1meanallsim", nameguardarallsim,".txt"))
    param2meanallsim<-scan(file=paste0(ruta,"/res/param2meanallsim", nameguardarallsim,".txt"))
    gammaallsim<- param2meanallsim+exp(-dist/param1meanallsim)
    print(c(param1meanallsim,param2meanallsim))

    param1quant1allsim<-scan(file=paste0(ruta,"/res/param1quant1allsim", nameguardarallsim,".txt"))
    param1quant2allsim<-scan(file=paste0(ruta,"/res/param1quant2allsim", nameguardarallsim,".txt"))

    param2quant1allsim<-scan(file=paste0(ruta,"/res/param2quant1allsim", nameguardarallsim,".txt"))
    param2quant2allsim<-scan(file=paste0(ruta,"/res/param2quant2allsim", nameguardarallsim,".txt"))


    # Data

    dat<-data.frame("dist"=dist,"1"=gamma,"2"=gammapluginaveragemles, "3"=gammaweightedateachdistance,
    "4"=gammaallsim)

    #apply melt to be able to plot all lines in the same plot
    names(dat)=c("dist", 1:(ncol(dat)-1))
    datmelted<-melt(dat, id.vars="dist")





    if(i!=5){
        # Plot transmissivity function gamma_{ij}
        titulo<-bquote(atop(paste(ind," ",.(i) ," ",sample," ",.(j),"\n ",
            phi==.(signif(estimates[1,1],3))," ",
            rho==.(signif(estimates[2,1],3))," "),
        paste(paste(phi,"w")==.(signif(param1weighted,3))," ",
            paste(rho,"w")==.(signif(param2weighted,3))," ",
            paste(phi,"all")==.(signif(param1meanallsim,3))," ",
            "(",.(signif(param1quant1allsim,3)),",",.(signif(param1quant2allsim,3)),") ",
            paste(rho,"all")==.(signif(param2meanallsim,3)),
            "(",.(signif(param2quant1allsim,3)),",",.(signif(param2quant2allsim,3)),") ")
        ))
    }
    if(i==5){
        titulo<-bquote(paste("all individuals, all samples \n",
        paste(phi,"w")==.(signif(param1weighted,3))," ",
        paste(rho,"w")==.(signif(param2weighted,3))
        ))#,




        titulo<-bquote(atop(paste(ind," ",.(i) ," ",sample," ",.(j),"\n "),
        paste(paste(phi,"w")==.(signif(param1weighted,3))," ",
            paste(rho,"w")==.(signif(param2weighted,3))," ",
            paste(phi,"all")==.(signif(param1meanallsim,3))," ",
            "(",.(signif(param1quant1allsim,3)),",",.(signif(param1quant2allsim,3)),") ",
            paste(rho,"all")==.(signif(param2meanallsim,3)),
            "(",.(signif(param2quant1allsim,3)),",",.(signif(param2quant2allsim,3)),") ")
        ))

    }
    gp<-ggplot(data=datmelted)

    if(i!=5){
        etiquetas<-c("mles","weighted avg mles","weighted avg gamma_ij at each distance","allsim")
        gp<-gp+
            geom_line(data=datmelted[which(datmelted$variable==1),],aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            geom_line(data=datmelted[which(datmelted$variable==2),],aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            geom_line(data=datmelted[which(datmelted$variable==3),],
            aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            geom_line(data=datmelted[which(datmelted$variable==4),],
            aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            scale_linetype_manual(name="Calculated with",labels=etiquetas,values = 1:4)+
            scale_size_manual(name="Calculated with",labels=etiquetas,values = c(1.2,1.1,1,1))+
            scale_color_manual(name="Calculated with",labels=etiquetas,values = c("black","red","green","orange"))
    }

    if(i==5){
        #gp<-gp+
        # geom_line(data=datmelted[which(datmelted$variable==2),],aes(x=dist, y=value, group=variable),size=1.2)


        etiquetas<-c("weighted avg mles","weighted avg gamma_ij at each distance","allsim")
        gp<-gp+
            geom_line(data=datmelted[which(datmelted$variable==2),],aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            geom_line(data=datmelted[which(datmelted$variable==3),],
            aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            geom_line(data=datmelted[which(datmelted$variable==4),],
            aes(x=dist, y=value, group=variable,linetype=variable,size=variable,color=variable))+
            scale_linetype_manual(name="Calculated with",labels=etiquetas,values = 2:4)+
            scale_size_manual(name="Calculated with",labels=etiquetas,values = c(1.1,1,1))+
            scale_color_manual(name="Calculated with",labels=etiquetas,values = c("red","green","orange"))

    }
    gp<-gp+theme_bw() +  ggtitle(titulo) +# theme(text=element_text(size = 8)) +
        theme(legend.position = "bottom", legend.direction = "horizontal")+
        scale_x_continuous(breaks= 0:10, labels= 0:10)+#, limits)+
        xlab("distance") + ylab(expression(paste("transmissivity   ",rho+exp(-distance/phi))))
        #+  coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept=rho, col="red") +
        #annotate("text", x = 0.1, y = .05, label = paste("rho == ", signif(rho,3)),parse=TRUE)
    png(paste0(ruta,"/plots/transmissitivity",nameguardarplot,".png"))
    plot(gp)
    dev.off()
}

# END PLOT Compare estimates obtained with bootstrap sampling and estimates obtained by weighted average
############################################################################
