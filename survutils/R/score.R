#' @export

score <- function(filename, numtimes)
{
  
  data <- read.csv(filename, sep=",")
  
  cols<-ncol(data)
  features<-data[,3:cols]
  features<-as.matrix(features)
  
  if(numtimes !=0 )
  {
    cutoffs <- sort(unique(data[which(data$event==1),1]))[1:numtimes]
    output <- pseudosurv(time=data$time,event=data$event, tmax=cutoffs)
  }
  else
  {
    output <- pseudosurv(time=data$time,event=data$event)
  }
  
  multi_output_survival <- as.matrix(output$pseudo)
  
  mrcefit <- mrce(Y=multi_output_survival, X=features, lam1=1e-3, lam2=1e-3, method="single", cov.tol=0.1, tol.out=1e-10)
  
  finalbeta <- as.matrix(mrcefit$Bhat)
  
  mrcepredictions <- features %*% finalbeta;
  mrcepredictions <- as.matrix(mrcepredictions)
  
  time <- data$time
  status <- data$event
  surv <- Surv(time, status)
  
  auc <- sapply(1:ncol(mrcepredictions), function(i){
    event_prediction <- 1-mrcepredictions[,i]
    survConcordance(surv ~ event_prediction)$concordance
  })
  
  brier <- sapply(1:ncol(mrcepredictions), function(i){
    event_prediction <- 1-mrcepredictions[,i]
    sum((status-event_prediction)^2)/nrow(data)
  })
  
  return(c(auc=mean(auc), auc.sd=sd(auc), brier.score=mean(brier)))
  
}