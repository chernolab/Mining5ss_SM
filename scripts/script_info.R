directory_fis = 'data/freqs/'
organisms = list.files(directory_fis)
organisms = organisms[which(nchar(organisms) == 3)]

fis_list = lapply(paste(directory_fis, organisms, '/fi.txt', sep=''), FUN=read.table)

sitewise.entropy = function(x){
  pilogpi = x * log(x)
  pilogpi = as.matrix(pilogpi)
  pilogpi[which(is.nan(pilogpi))] = 0
  return(-rowSums(pilogpi))
}

info_list = lapply(fis_list, FUN=sitewise.entropy)
info_mat = 2 - matrix(unlist(info_list), nrow=9)
rownames(info_mat) = c(-(3:1),1:6)
colnames(info_mat) = organisms
