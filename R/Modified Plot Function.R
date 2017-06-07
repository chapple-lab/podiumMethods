plot.OSC.results<-function(obj,plot="scores",groups=NULL,lgroups=NULL,nOSCcomp=NULL,alpha=1,glines=F){
  require(ggplot2)
  #plot = one of: c("RMSEP","scores","loadings","multiLoadings","delta.weights,"summary")
  #groups is a factor to show group visuyalization in scores plot
  switch(plot,
         RMSEP   		=  .local<-function(obj){
           #bind info and RMSEP
           comps<-obj$total.LVs
           ocomps<-obj$OSC.LVs
           plot.obj<-obj$RMSEP
           bound<-do.call("rbind",lapply(1:length(comps),function(i)
           {
             out<-as.data.frame(cbind(plot.obj[[i]][,1],c(0:comps[i]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
             colnames(out)<-c("RMSEP","component","model")
             out
           }))
           bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))

           #custom theme
           .theme<- theme(
             axis.line = element_line(colour = 'gray', size = .75),
             panel.background = element_blank(),
             plot.background = element_blank()
           )
           #plot
           p<-ggplot(data=bound, aes(x=component, y=RMSEP, group=model,color=model)) + geom_line(size=1,alpha=.5) + geom_point(size=2,alpha=alpha)+.theme
           print(p)
         },
#Multi Scores
         scores 			=	.local<-function(obj){
           comps<-obj$total.LVs
           ocomps<-obj$OSC.LVs
           plot.obj<-obj$scores
           nullGroups=F

           if(is.null(groups))
           {
             nullGroups=T
             groups<-rep("No Groups",nrow(plot.obj[[1]][,]))
           }

           bound<-do.call("rbind",lapply(1:length(comps),function(i)
           {
             out<-as.data.frame(cbind(plot.obj[[i]][,1:2],unlist(groups),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
             colnames(out)<-c("Comp1","Comp2","groups","model")
             out
           }))
           bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))

           if(!nullGroups)
           {
             #calculate convex hull for polygons for each group
             data.obj <- split(bound, bound$model)
             tmp.obj <- lapply(1:length(data.obj), function(i){
               obj<-data.obj[[i]]
               s2<-split(obj,obj[,3])
               do.call(rbind,lapply(1:length(s2),function(j){
                 tmp<-s2[[j]]
                 tmp[chull(tmp[,1:2]),]
               }))
             })
             chull.boundaries <- do.call("rbind", tmp.obj)
           }

           #custom theme
           if(glines)
           {
             .theme<- theme(
               axis.line = element_line(colour = 'gray', size = .75),
               panel.background = element_blank(),
               panel.border = element_rect(colour="gray",fill=NA),
               panel.grid.minor = element_line(colour = "gray80", linetype = "dotted"),
               panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
               plot.background = element_blank()
             )}
           else
           {
             .theme<- theme(
               axis.line = element_line(colour = 'gray', size = .75),
               panel.background = element_blank(),
               panel.border = element_rect(colour="gray",fill=NA),
               plot.background = element_blank()
             )}

           if(!is.null(nOSCcomp))
           {
             models = split(bound,bound$model)
             bound = models[[(nOSCcomp+1)]]
             bound$model = factor(bound$model)
             bounds = split(chull.boundaries,chull.boundaries$model)
             chull.boundaries = bounds[[(nOSCcomp+1)]]
#              chull.boundaries$model = factor(chull.boundaries$model)
           }

           #make plot
           p<-ggplot(data=bound, aes(x=Comp1, y=Comp2, group=groups,color=groups)) + #geom_density2d(aes(group=groups))+
             geom_hline(aes(yintercept=0),color="gray60",linetype="dashed")+geom_vline(aes(xintercept=0),color=I("gray60"),linetype=2)+facet_grid(. ~ model)

           if(nullGroups) { p<-p+geom_point(size=2)+.theme }
           else { p<-p+geom_polygon(data=chull.boundaries,aes(x=Comp1,y=Comp2,fill=groups),alpha=.5) +geom_point(size=2,alpha=alpha)+.theme }

           print(p)
         },
         loadings 		= 	.local<-function(obj){ # will only plot first component for each model
           comps<-obj$total.LVs
           ocomps<-obj$OSC.LVs
           plot.obj<-obj$loadings
           bound<-do.call("rbind",lapply(1:length(comps),function(i)
           {
             out<-as.data.frame(cbind(plot.obj[[i]][,1:2],rownames(plot.obj[[i]]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
             colnames(out)<-c("Comp1","Comp2","variable","model")
             out
           }))
           bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))

           #custom theme
           .theme<- theme(
             axis.line = element_line(colour = 'gray', size = .75),
             panel.background = element_blank(),
             legend.position = "none",
             plot.background = element_blank()
           )

           #make plot
           p<-ggplot(data=bound, aes(x=variable,y=Comp1, fill=variable)) + geom_bar(stat = "identity") + coord_flip() + #geom_density2d(aes(group=groups))+
             facet_grid(. ~ model) +.theme
           print(p)
         },
#Multi Loadings
        multiLoadings   		=	.local<-function(obj){
          comps<-obj$total.LVs
          ocomps<-obj$OSC.LVs
          plot.obj<-obj$loadings
          if(is.null(groups)){groups<-rep("No Groups",nrow(plot.obj[[1]][,]))}
          bound<-do.call("rbind",lapply(1:length(comps),function(i)
          {
            out<-as.data.frame(cbind(plot.obj[[i]][,1:2],unlist(groups),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
            colnames(out)<-c("Comp1","Comp2","groups","model")
            out
          }))
          bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))


          #custom theme
          if(glines)
          {
          .theme<- theme(
            axis.line = element_line(colour = 'gray', size = .75),
            panel.background = element_blank(),
            panel.border = element_rect(colour="gray",fill=NA),
            panel.grid.minor = element_line(colour = "gray80", linetype = "dotted"),
            panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
            plot.background = element_blank()
          )}
          else
          {
            .theme<- theme(
              axis.line = element_line(colour = 'gray', size = .75),
              panel.background = element_blank(),
              panel.border = element_rect(colour="gray",fill=NA),
              plot.background = element_blank()
          )}

          if(!is.null(nOSCcomp))
          {
            models = split(bound,bound$model)
            bound = models[[(nOSCcomp+1)]]
            bound$model = factor(bound$model)
          }

          #make plot
          p<-ggplot(data=bound, aes(x=Comp1, y=Comp2, group=groups,color=groups))+
            geom_hline(aes(yintercept=0),color="gray60",linetype="dashed")+geom_vline(aes(xintercept=0),color=I("gray60"),linetype=2)+facet_grid(. ~ model)

          p<-p+geom_point(size=2,alpha=alpha)+.theme
          print(p)
        },
         delta.weights 	= 	.local<-function(obj){ # will only plot first component for each model
           comps<-obj$total.LVs
           ocomps<-obj$OSC.LVs
           plot.obj<-obj$loading.weights
           bound<-do.call("rbind",lapply(2:(length(ocomps)),function(i)
           {
             out<-as.data.frame(cbind(plot.obj[[1]][,1]-plot.obj[[i]][,1],names(plot.obj[[i]][,1]),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep="")))
             colnames(out)<-c("delta_weight","variable","model")
             out
           }))
           bound[,1]<-as.numeric(as.matrix(bound[,1]))

           #theme
           .theme<- theme(
             axis.line = element_line(colour = 'gray', size = .75),
             panel.background = element_blank(),
             legend.position = "none",
             plot.background = element_blank()
           )
           #make plot
           p<-ggplot(data=bound, aes(x=variable,y=delta_weight, fill=variable)) + geom_bar(stat = "identity") + coord_flip() + #geom_density2d(aes(group=groups))+
             facet_grid(. ~ model) +.theme
           print(p)
         },
        summary   		=	.local<-function(obj){
          comps<-obj$total.LVs
          ocomps<-obj$OSC.LVs
          plot.obj<-obj$scores
          plot.obj2 = obj$loadings
          nullGroups=F

          if(is.null(lgroups)) {lgroups = rep("Loadings",nrow(plot.obj2[[1]][,]))}

          if(is.null(groups))
          {
            nullGroups=T
            groups<-rep("No Groups",nrow(plot.obj[[1]][,]))
          }

          bound<-do.call("rbind",lapply(1:length(comps),function(i)
          {
            out<-as.data.frame(cbind(plot.obj[[i]][,1:2],unlist(groups),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs",sep=""),i))
            out2<-as.data.frame(cbind(plot.obj2[[i]][,1:2],unlist(lgroups),paste(comps[i]," LVs and ",ocomps[i]," OSC LVs Loadings",sep=""),i))
            out = rbind(out,out2)
            colnames(out)<-c("Comp1","Comp2","groups","model","plotCode")
            out
          }))
#           boundOut <<- bound
          bound[,1:2]<-as.numeric(as.matrix(bound[,1:2]))

          if(!nullGroups)
          {
            #calculate convex hull for polygons for each group
            filtered = bound[!grepl("Loadings",bound$model),] #remove Loadings from data
            data.obj <- split(filtered, filtered$model)
            tmp.obj <- lapply(1:length(data.obj), function(i){
              obj<-data.obj[[i]]
              s2<-split(obj,obj[,3])
              do.call(rbind,lapply(1:length(s2),function(j){
                tmp<-s2[[j]]
                tmp[chull(tmp[,1:2]),]
              }))
            })
            chull.boundaries <- do.call("rbind", tmp.obj)
          }

          #custom theme
          if(glines)
          {
            .theme<- theme(
              axis.line = element_line(colour = 'gray', size = .75),
              panel.background = element_blank(),
              panel.border = element_rect(colour="gray",fill=NA),
              panel.grid.minor = element_line(colour = "gray80", linetype = "dotted"),
              panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
              plot.background = element_blank()
            )}
          else
          {
            .theme<- theme(
              axis.line = element_line(colour = 'gray', size = .75),
              panel.background = element_blank(),
              panel.border = element_rect(colour="gray",fill=NA),
              plot.background = element_blank()
            )}

          if(!is.null(nOSCcomp))
          {
            models = split(bound,bound$model)
            bound = models[[(nOSCcomp+1)]]
            bound$model = factor(bound$model)
            bounds = split(chull.boundaries,chull.boundaries$model)
            chull.boundaries = bounds[[(nOSCcomp+1)]]
            #              chull.boundaries$model = factor(chull.boundaries$model)
          }

          #prep data
          loadingData = bound[grepl("Loadings",bound$model),]
          scoreData = bound[!grepl("Loadings",bound$model),]


          #make plot
          p<-ggplot(data=bound, aes(x=Comp1, y=Comp2, group=groups,color=groups, fill=groups)) + #geom_density2d(aes(group=groups))+
            geom_hline(aes(yintercept=0),color="gray60",linetype="dashed")+geom_vline(aes(xintercept=0),color=I("gray60"),linetype=2)+
            facet_wrap(plotCode ~ model,scales="free") + geom_point(data=loadingData,size=2,alpha=alpha)


          if(nullGroups) { p<-p+geom_point(data=scoreData,size=2)+.theme }
          else { p<-p+geom_polygon(data=chull.boundaries,alpha=.5)+geom_point(data=scoreData,size=2,alpha=1)+.theme }

          print(p)
        }
  )
  .local(obj)
}
