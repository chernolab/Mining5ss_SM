sp = c('apl', 'bta', 'clu', 'eca', 'mdo', 'oan', 'ocu', 'sha', 'ssa', 'ssc', 'xtr')
logos = lapply(sp, FUN=function(x){return(readRDS(paste('data/seqs/logo_', x, '.rds', sep='')))})
aux = lapply(logos, FUN=function(x){return(c(nrow(x), nrow(unique(x)), nrow(unique(x[which(x[,4] == 'G' & x[,5] == 'T'),]))))})
print(lapply(aux, FUN=function(x){return(paste0(x, collapse=' & '))}))
