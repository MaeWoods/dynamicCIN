library(gtools)
library(ggplot2)
source("multiplot.R")

#####
## Parameters
#####

n.resim <- 100
path <- "./simulated-data-sets/"
path.resim <- "./resimulated-data-sets/"


wass.dist <- function(d1,d2){
   qpoints <- seq(0.05,0.95,0.01)
   q1 <- quantile( d1, prob=qpoints )
   q2 <- quantile( d2, prob=qpoints )
   del <- (q1-q2)
   return( sum(del^2) )
}

if(0){
   #####
   ## Access and plot the data
   #####
   data <- read.table("data-divergence-gastric.txt",header=T)$pdiv

   pdf("plot-data-density.pdf") 
   dp <- data.frame( divergence=data )
   p <- ggplot(dp, aes(x = divergence)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + xlim(-2, 5) + xlab("fractional chromosomal divergence") + ylab("")
   print(p)
   dev.off()


   #####
   ## Access the simulations and do the ABC
   #####
   files <- list.files(path=path, pattern="results-S*")
   files <- mixedsort(sort(files))
   nfiles <- length(files)
   cat("n particles total:", nfiles, "\n", sep=" ")

   dists <- rep(0,nfiles)
   for(i in 1:nfiles){
      f <- paste(path,files[i],sep="")
      sim <- read.csv(f)$RelLen
      #dists[i] <- ks.test(sim,data)$statistic
      dists[i] <- wass.dist(sim,data)
   }

   pdf("plot-fit-distances.pdf")
   dp <- data.frame( distance=dists )
   p <- ggplot(dp, aes(x = distance)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + ylab("") 
   print(p)
   dev.off()

   #####
   ## Plot the best 10 fitting simulated data sets
   #####

   top.sims <- order(dists)[1:10]

   plots <- vector("list",10)
   for(i in 1:10){
      f <- paste(path,files[top.sims[i]],sep="")
      sim <- read.csv(f)$RelLen
     
      dp <- data.frame( divergence=c(data,sim), y=c(rep("data",length(data)), rep("sim",length(sim)) ) )
      plots[[i]] <- ggplot(dp, aes(x = divergence, fill = y)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + xlim(-2, 5) + xlab("") + ylab("")
   }

   pdf("plot-fit-top-sims.pdf",height=6,width=15)
   multiplot(plotlist=plots, cols=5)
   dev.off()
}

if(0){
   #####
   ## Extract the parameters that gave rise to these
   #####
   top.sims <- top.sims[1:5]

   top.pars <- gsub("results","parameters", files[top.sims] )
   top.pars.files <-  paste(path,top.pars,sep='')

   cat("running cat of parameter files\n")
   comm1 <-  paste("cat", paste(top.pars.files,collapse=" "), " > posterior.dat", collapse="" )
   system(comm1)

   cat("running resimulation\n")
   comm2 <- paste("mkdir resimulated-data-sets; cd resimulated-data-sets; ../../dynamicCIN/simulation/CINmodel Particles=",n.resim," Model=2 Resim=../posterior.dat; cd ../",sep="")
   system(comm2)
}


if(0){
   #####
   ## Resimulations: posterior predictive distribution
   #####
   files.resim <- list.files(path=path.resim, pattern="results-resim-S*")

   dens.pred.x <- c()
   dens.pred.y <- c()

   all.sim <- c()

   for(i in 1:n.resim){
      f <- paste(path.resim,files.resim[i],sep="")
      sim <- read.csv(f)$RelLen

      # **** hack ****
      sim <- sim[ which(sim < 5 & sim > -2) ]
      #print( summary(sim) )

      all.sim <- cbind(all.sim,sim)

      dens <- density(sim, from=-2, to=10)
      dens.pred.x <- rbind(dens.pred.x, dens$x)
      dens.pred.y <- rbind(dens.pred.y, dens$y)
   }

   dens.data <- density(data)
   dens.pred.q <- apply(dens.pred.y,MARGIN=2,function(x){ quantile(x,p=c(0.025,0.5,0.975) ) })

   xvals <- dens.pred.x[1,] 
   dp <- data.frame( x=c(dens.data$x, xvals ), y=c(dens.data$y,dens.pred.q[2,]), z=c(rep(1,512),rep(2,512) ) )
   #dp2 <- data.frame( x=c(xvals,xvals), y=c(dens.pred.q[1,],dens.pred.q[3,]), z=c(rep("0.025",512),rep("0.975",512))  )
   dp3 <- data.frame( x=xvals, ymin=dens.pred.q[1,], ymax=dens.pred.q[3,] )

   p <- ggplot() + geom_area(aes(x=x,y=y, fill=factor(z)), data=dp, stat="identity", alpha = 0.5, colour="black") + guides(fill=FALSE) + xlim(-2, 5) + xlab("fractional chromosomal divergence") + ylab("")
   p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax), data=dp3, linetype="dotted", alpha=0.3)

   pdf("plot-fit-post-pred-density.pdf")
   print(p)
   dev.off()
}

if(1){
   
   top.sims <- order(dists)[1:40]

   top.pars <- gsub("results","parameters", files[top.sims] )
   top.pars.files <-  paste(path,top.pars,sep='')

   cat("running cat of parameter files\n")
   comm3 <-  paste("cat", paste(top.pars.files,collapse=" "), " > posterior.dat", collapse="" )
   system(comm3)

   post <- read.table("posterior.dat",header=F)
   names(post) <- c("Npop","ngen","g_d","mu_i","p_trans","svp2","trp2","SV_min","svgradient","svtransgrad","threshold_g","SV_max","mu_i_range","gnmdou","maxchr","minchr")

   plots <- vector("list",3)
   dp <- data.frame( p_g=post$gnmdou, p_t=post$p_trans, p_v=post$g_d )
   plots[[1]] <- ggplot(dp, aes(x = p_v)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + xlim(0, 1) + xlab("p(insertion/deletion)") + ylab("")
   plots[[2]] <- ggplot(dp, aes(x = p_t)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + xlim(0, 1) + xlab("p(translocation)") + ylab("") 
   plots[[3]] <- ggplot(dp, aes(x = p_g)) + geom_density(alpha = 0.5) + guides(fill=FALSE) + xlim(0, 1) + xlab("p(genome doubling)") + ylab("")
   
   pdf("plot-posterior-1D.pdf",height=3,width=9)
   multiplot(plotlist=plots, cols=3)
   dev.off()
}