library(gtools)
library(ggplot2)
source("multiplot.R")

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
   path <- "./simulated-data-sets/"

   files <- list.files(path=path, pattern="results-S*")
   files <- mixedsort(sort(files))
   nfiles <- length(files)
   cat("n particles total:", nfiles, "\n", sep=" ")

   dists <- rep(0,nfiles)
   for(i in 1:nfiles){
      f <- paste(path,files[i],sep="")
      sim <- read.csv(f)$RelLen
      dists[i] <- ks.test(sim,data)$statistic
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

   top.pars <- gsub("results","parameters", files[top.sims] )
   top.pars.files <-  paste(path,top.pars,sep='')

   cat("running cat of parameter files\n")
   comm1 <-  paste("cat", paste(top.pars.files,collapse=" "), " > posterior.dat", collapse="" )
   system(comm1)

   cat("running resimulation\n")
   comm2 <- "mkdir resimulated-data-sets; cd resimulated-data-sets; ../../dynamicCIN/simulation/CINmodel Particles=100 Model=2 Resim=../posterior.dat; cd ../"
   system(comm2)
}

path.resim <- "./resimulated-data-sets/"
files.resim <- list.files(path=path.resim, pattern="results-resim-S*")

dens.pred.x <- c()
dens.pred.y <- c()

all.sim <- c()

for(i in 1:100){
   f <- paste(path.resim,files.resim[i],sep="")
   sim <- read.csv(f)$RelLen

   # hack
   sim <- sim[ which(sim < 5 & sim > -2) ]
   print( summary(sim) )

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
#p <- p + geom_line(aes(x=x,y=y), linetype="dotted", data=dp2, alpha=0.2)
p <- p + geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax), data=dp3, linetype="dotted", alpha=0.1)
#p <- p + scale_fill_manual( values = c("red","blue", "green") )

pdf("plot-fit-post-pred-density.pdf")
print(p)
dev.off()
