library("readxl")
library(ggplot2)
library(reshape2)
library(parallel)
library(foreach)
library(doParallel)
library(fields)

Euclideandistance <- function(x1, y1, x2, y2){
  return sqrt((x1-x2)²+(y1-y2)²)
}


no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)


labels <- read_excel("training_data/training/labels.xlsx")

trainingdata = list.files('training_data/training', "*.csv") 

foreach(i = 1:(length(trainingdata)))%dopar%{
  assign(trainingdata[i], as.matrix(read.csv(paste('training_data/training/', trainingdata[i], sep=''))[130:2630,4:302]))
  }

sample(which(labels$candy == "halls"), 1)

Grad = colorRampPalette(c('grey',"blue","red","black"))
{
  dataHal = paste0(labels[sample(which(labels$candy == "halls"), 1),1])
  par(mar = c(0,0,1,0))
  image(get(dataHal), useRaster=TRUE, axes=FALSE, col = Grad(10), main = "Halls candy")
  image.plot(legend.only=TRUE, zlim=c(min(get(dataHal)),max(get(dataHal))),col = Grad(100))
}


Grad = colorRampPalette(c('grey',"blue","red","black"))
{
  dataHal = paste0(labels[sample(which(labels$candy == "citrus"), 1),1])
  par(mar = c(0,0,1,0))
  image(get(dataHal), useRaster=TRUE, axes=FALSE, col = Grad(10), main = "citrus candy")
  image.plot(legend.only=TRUE, zlim=c(min(get(dataHal)),max(get(dataHal))),col = Grad(100))
}

  
  trainingdata_peax = list.files('data/training', "*.csv")
  
  foreach(i = 1:(length(trainingdata_peax)))%dopar%{
    assign(trainingdata_peax[i], read.csv(paste('training_data/training/', trainingdata_peax[i], sep='')))
  }
  
  





