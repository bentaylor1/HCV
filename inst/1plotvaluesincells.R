
# PLOT hcv.level, timeinfected and ranktimeinfected
#######################################

# Samples are 10 $\times$ 10 squared cells.
# Plot 1 (infected), 0 (no infected)
#m = matrix(rbinom(n=100,size=1,prob=.5),10,10)
#image(m, col=c("black","gray"), axes=FALSE)

# CHOOSE
vble<-"hcv.level"
vble<-"ranktimeinfected"
vble<-"timeinfected"
vble <- "ifit.level"
vble <- "logifitlevel"
vble <- "roundlogifitlevel"

####################################

# same colours for all samples
# if it is 0, put white

# choose colours for real data and synthetic data
colsleyenda <- c("grey",rev(viridis(256, 1)))
if(vble=="hcv.level"){
    maxvbleallindandsamples<-51
    maxvbleallindandsampleshcv<-51

    if(typedata == "synthetic"){
        maxvbleallindandsamples<-62
        maxvbleallindandsampleshcv<-62
    }
}


if(vble=="ifit.level"){
    maxvbleallindandsamples<-55364
}

if(vble=="logifitlevel"){
    maxvbleallindandsamples<-11
}




if(vble=="roundlogifitlevel"){
    maxvbleallindandsamples<-11
}


if(vble=="ranktimeinfected"){
    maxvbleallindandsamples<-61
    if(typedata == "synthetic"){
        maxvbleallindandsamples<-61
    }
}

if(vble=="timeinfected"){
    maxvbleallindandsamples<-52
    if(typedata == "synthetic"){
        maxvbleallindandsamples<-63
    }
}

###########################################


pl<-list()
l<-1
for(i in 1:length(sort(unique(dwhole$subid)))){

    if(vble=="hcv.level"){
        # Calculate timeinfected and ranktimeinfected taking into account data in dwhole
        dwhole<-fnAddColumnsTimeinfectedAndRankTimeinfected(dwhole)
        maxvble<-max(dwhole[,vble])
        print(maxvble)
    }

    numsamples<-3
    if(i==4){
        numsamples<-2
    }

    # If synthetic (subid = 20000), numsamples is 5
    if(sort(unique(dwhole$subid))[i] == "20000"){
        numsamples<-5
    }

    for(j in 1:numsamples){
        # Select one individual and sample
        d<-dwhole[which(dwhole$subid==sort(unique(dwhole$subid))[i] & dwhole$gridN==j),]

        if(vble!="hcv.level"){
            # Calculate timeinfected and ranktimeinfected taking into account data in d
            d<-fnAddColumnsTimeinfectedAndRankTimeinfected(d)
            maxvble<-max(d[,vble])
            print(maxvble)
        }


        # need to make sure xpos and ypos ordered
        m<-d[,vble]

        mm<-structure(m, .Dim=c(10L,10L), .Dimnames = list(paste("y", 1:10), paste("x", 1:10)))
        melted <- melt(mm)

        # legend colours defined at the beginning. Put same colours for all samples. If it is 0, put white.
        pl[[l]] <- ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile(color="black") +
          theme_void()+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(),
                              axis.title.x=element_blank(), axis.title.y=element_blank()) +
          geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
          #scale_fill_gradient(low= "white", high="grey")
          #scale_fill_viridis(low= "white", limits=c(0, 51))
          scale_fill_gradientn(colours = colsleyenda,limits=c(0, maxvbleallindandsamples))

        l<-l+1
    }

}


# I do this to get legend
pl[[l]] <- ggplot(melted, aes(x = Var2, y = Var1, fill = value)) + geom_tile(color="black") +
    theme_void()+
    theme(legend.position="left",axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    scale_fill_gradientn(colours = colsleyenda, limits=c(0, maxvbleallindandsamples)) +
    coord_equal()

gg <- ggplotGrob(pl[[l]]  + theme(legend.position = "left"))$grobs
leyenda <- gg[[which(sapply(gg, function(x) x$name) == "guide-box")]]


# Use coord_equal() to get geom_tile() to draw square rather than rectangular cells

# real data
if(unique(d$subid) != "20000"){
    png(paste0(ruta,"/plots/","cells",sub("\\.", "", vble),".png"), width=800, height=800)
    grid.arrange(pl[[1]]+coord_equal(), pl[[2]]+coord_equal(), pl[[3]]+coord_equal(),
                 pl[[4]]+coord_equal(), pl[[5]]+coord_equal(), pl[[6]]+coord_equal(),
                 pl[[7]]+coord_equal(), pl[[8]]+coord_equal(), pl[[9]]+coord_equal(),
                 pl[[10]]+coord_equal(), pl[[11]]+coord_equal(),
                 leyenda, ncol = 3, name="dd")
    dev.off()
}

# synthetic data
if(unique(d$subid) == "20000"){
    png(paste0(ruta,"/plots/","cells",sub("\\.", "", vble),".png"), width=800, height=800)
    grid.arrange(pl[[1]]+coord_equal(), pl[[2]]+coord_equal(), pl[[3]]+coord_equal(),
               pl[[4]]+coord_equal(), pl[[5]]+coord_equal(),
               leyenda, ncol = 3, name="dd")
    dev.off()

}

#http://rstudio-pubs-static.s3.amazonaws.com/2852_379274d7c5734f979e106dcf019ec46c.html


# Histogram of HCV levels (removing levels that are equal to 0)
d2<-dwhole[which(dwhole$hcv.level!=0),]
png(paste0(ruta,"/plots/","hist","hcvlevelnot0.png"), width=800, height=800)
ggplot(data=d2, aes(hcv.level)) +
    theme_bw() + theme(legend.position = "none") +
    geom_histogram(breaks=seq(0, maxvbleallindandsampleshcv, by=1), fill="grey", col="black", aes(fill=..count..)) +
    labs(x="HCV level", y="") +
    facet_grid(subid ~ gridN)
    #scale_fill_gradient("Count", low = "green", high = "red")
dev.off()
