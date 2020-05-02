library(parallel)


cl <- makeCluster(38)

parLapply(cl = cl,rep(list(7),38),function(x){
  
  save.file <- "galgo.search.parallel.Rdata"
  
  library(galgo)
  loadObject(save.file)
  assignParallelFile(galgo.search)
  blast(galgo.search)
  
  
})


stopCluster(cl)