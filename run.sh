#!/usr/bin/env Rscript

library(stringr)

read_params = function(filepath) {
  params=list()
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    spl=stringr::str_split(line,'=')
    print(spl[[1]])
    if (length(spl[[1]])>1){
        params=c(params,list(spl[[1]][length(spl[[1]])]))
    }
  }
  close(con)
  return(params)
}

params=read_params('./params.txt')

source('./climorec.r')
results=apply_rec(params[[1]],params[[3]],params[[4]],params[[5]],params[[6]],params[[8]],params[[7]],params[[9]],params[[10]],params[[11]],params[[12]])

setwd(params[[1]])
dir.create(params[[2]])
setwd(paste(params[[1]],params[[2]],sep='/'))

write.csv(results[[1]],paste(params[[2]],'_renormalized_reconstruction.csv',sep=''),quote=F,row.names=F)

write.csv(results[[2]],paste(params[[2]],'_original_reconstruction.csv',sep=''),quote=F,row.names=F)                                                                            

write.csv(results[[3]],paste(params[[2]],'_correlation_scores.csv',sep=''),quote=F,row.names=F)                                                                                

write.csv(results[[4]],paste(params[[2]],'_RMSE_scores.csv',sep=''),quote=F,row.names=F) 

write.csv(results[[5]],paste(params[[2]],'_nb_records.csv',sep=''),quote=F,row.names=F)

write.csv(results[[6]],paste(params[[2]],'_name_proxies.csv',sep=''),quote=F,row.names=F)

write.csv(results[[7]],paste(params[[2]],'_individual_reconstructions.csv',sep=''),quote=F,row.names=F)

write.csv(results[[8]],paste(params[[2]],'_NSCE_scores.csv',sep=''),quote=F,row.names=F)

write.csv(results[[9]],paste(params[[2]],'_individual_se.csv',sep=''),quote=F,row.names=F)

write.csv(results[[10]],paste(params[[2]],'_training_samples.csv',sep=''),quote=F,row.names=F)

write.csv(results[[11]],paste(params[[2]],'_testing_samples.csv',sep=''),quote=F,row.names=F)

write.csv(results[[12]],paste(params[[2]],'_ShapiroWilk_pvalues_residuals.csv',sep=''),quote=F,row.names=F)