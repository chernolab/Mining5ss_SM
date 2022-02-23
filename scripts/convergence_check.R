directory = 'runs/'
data_dir = 'data/'
species = c('ath', 'mtr', 'osa', 'ppa', 'ptri', 'hvu', 'sly', 'vvi',
            'dme', 'dps', 'dya', 
            'cel', 'cbr',
            'cne', 'ncr', 'mor', 'cci', 'ani',
            'apl', 'bta', 'clu', 'dre', 'eca', 'ggo', 'hsa', 'mdo', 'mmu', 'oan', 'ocu', 'ptr', 'sha', 'ssa', 'ssc', 'xtr',
            'metazoa', 'plantae', 'fungi')
gammas = c('0.010000', '0.015000', '0.020000', '0.025000', '0.030000', '0.035000', '0.040000', '0.045000', '0.050000')
plot_conv = F

fis = matrix(data=numeric(), nrow=0, ncol=36)
for(sp in species){
  print(sp)
  files = lapply(gammas, FUN=function(x){return(list.files(paste(directory, sp, '/gamma', x, sep='')))})
  iterations = lapply(files, FUN=function(x){return(unique(gsub(x, pattern='energy|Param.|Pi|Pij|sequences|\\.gz', replacement='')))})
  iterations = lapply(iterations, FUN=function(x){return(as.character(sort(as.numeric(x))))})
  
  for(g in gammas){
    fis = rbind(fis, do.call('rbind', lapply(paste(directory, sp, '/gamma', g, '/Pi', iterations[[which(gammas == g)]], sep=''), FUN=function(x){return(unlist(read.table(x)))})))
  }
  
  if(plot_conv){
    plot(1:nrow(fis), fis[,1], col=1, type='l', ylim=c(min(unlist(fis)), max(unlist(fis))))
    for(i in 2:ncol(fis)){
      lines(1:nrow(fis), fis[,i], col=i)
    }
  }
    
  pis = unlist(read.table(paste(data_dir, 'freqs/', sp, '/fi.txt', sep='')))
  if(plot_conv){
    plot(1:nrow(fis), abs(fis[,1] - pis[1]), col=1, log='y', type='l', ylim=c(min(abs(unlist(fis) - pis[1])), abs(max(unlist(fis) - pis[1]))))
    for(i in 2:ncol(fis)){
      lines(1:nrow(fis), abs(fis[,i] - pis[i]), col=i)
    }
  }
  print(max(abs(fis[nrow(fis),] - pis)))
}

jijs = matrix(data=numeric(), nrow=0, ncol=36**2)
for(sp in species){
  files = lapply(gammas, FUN=function(x){return(list.files(paste(directory, sp, '/gamma', x, sep='')))})
  iterations = lapply(files, FUN=function(x){return(unique(gsub(x, pattern='energy|Param.|Pi|Pij|sequences', replacement='')))})
  iterations = lapply(iterations, FUN=function(x){return(as.character(sort(as.numeric(x))))})
  
  for(g in gammas){
    jijs = rbind(jijs, do.call('rbind', lapply(paste(directory, sp, '/gamma', g, '/ParamJ', iterations[[which(gammas == g)]], sep=''), FUN=function(x){return(unlist(read.table(x)))})))
  }
  
  nonzero_indices = which(colSums(jijs) != 0)
  if(plot_conv){
    plot(1:nrow(jijs), jijs[,nonzero_indices[1]], col=1, type='l', ylim=c(min(unlist(jijs)), max(unlist(jijs))))
    for(i in nonzero_indices[2:ncol(jijs)]){
      lines(1:nrow(jijs), jijs[,i], col=i)
    }
  }
}

fijs = matrix(data=numeric(), nrow=0, ncol=36*36)
for(sp in species){
  files = lapply(gammas, FUN=function(x){return(list.files(paste(directory, sp, '/gamma', x, sep='')))})
  iterations = lapply(files, FUN=function(x){return(unique(gsub(x, pattern='energy|Param.|Pi|Pij|sequences|\\.gz', replacement='')))})
  iterations = lapply(iterations, FUN=function(x){return(as.character(sort(as.numeric(x))))})
  
  for(g in gammas){
    fijs = rbind(fijs, do.call('rbind', lapply(paste(directory, sp, '/gamma', g, '/Pij', iterations[[which(gammas == g)]], sep=''), FUN=function(x){return(unlist(read.table(x)))})))
  }
  
  pijs = unlist(read.table(paste(data_dir, 'freqs/', sp, '/fij.txt', sep='')))
  if(plot_conv){
    plot(1:nrow(fis), abs(fijs[,1] - pijs[1]), col=1, type='l', ylim=c(min(abs(unlist(fijs) - pijs[1])), 0.1))
    for(i in 2:ncol(fijs)){
      lines(1:nrow(fijs), abs(fijs[,i] - pijs[i]), col=i)
    }
  }
  print(max(abs(fijs[nrow(fijs),] - pijs)))
}
