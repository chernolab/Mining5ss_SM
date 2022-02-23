run_dir = './runs/ani/'
data_dir = '../data/freqs/zma/'
gamma = 0.05
full_dir = paste(run_dir, 'gamma', format(gamma, nsmall=6), sep='')

iters = sort(as.numeric(gsub(dir(full_dir, pattern='ParamH'), pattern='ParamH', replacement='')))
iters = iters[which(sapply(paste(full_dir, '/Pi', iters, sep=''), FUN=file.size) > 0)]
iters = iters[which(iters < 100)]

#Lo hacemos con Pi.
pis = lapply(iters, FUN=function(i){return(read.table(paste(full_dir, '/Pi', i, sep='')))})
fis = read.table(paste(data_dir, 'fi.txt', sep=''))
pis_mat = array(as.numeric(unlist(pis)), dim=c(nrow(pis[[1]]), ncol(pis[[1]]), length(pis)))

for(j in 1:nrow(fis)){
  for(k in 1:ncol(fis)){
    if(j==1 & k==1){
      plot(iters, pis_mat[j,k,] - fis[j,k], type='l', col=j+nrow(fis)*(k-1), ylim=c(-1,1))
    } else{
      lines(iters, pis_mat[j,k,] - fis[j,k], col=j+nrow(fis)*(k-1))
    }
  }
}

#Lo hacemos con Pij.
pijs = lapply(iters, FUN=function(i){return(read.table(paste(full_dir, '/Pij', i, sep='')))})
fijs = read.table(paste(data_dir, 'fij.txt', sep=''))
pijs_mat = array(as.numeric(unlist(pijs)), dim=c(nrow(pijs[[1]]), ncol(pijs[[1]]), length(pijs)))

for(j in 1:nrow(fijs)){
  for(k in 1:ncol(fijs)){
    if(j==1 & k==1){
      plot(iters, pijs_mat[j,k,] - fijs[j,k], type='l', col=j+nrow(fijs)*(k-1), ylim=c(-1,1))
    } else{
      lines(iters, pijs_mat[j,k,] - fijs[j,k], col=j+nrow(fijs)*(k-1))
    }
  }
}

#Lo hacemos con Jij.
jijs = lapply(iters, FUN=function(i){return(read.table(paste(full_dir, '/ParamJ', i, sep='')))})
jijs_mat = array(as.numeric(unlist(jijs)), dim=c(nrow(jijs[[1]]), ncol(jijs[[1]]), length(jijs)))

for(j in 1:nrow(jijs[[1]])){
  for(k in 1:ncol(jijs[[1]])){
    if(j==1 & k==1){
      plot(iters, jijs_mat[j,k,], type='l', col=j+nrow(jijs[[1]])*(k-1), ylim=c(-10,10))
    } else{
      lines(iters, jijs_mat[j,k,], col=j+nrow(jijs[[1]])*(k-1))
    }
  }
}
