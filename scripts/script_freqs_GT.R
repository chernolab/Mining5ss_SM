species_list = sp_list = c('metazoa', 'plantae', 'fungi')#c('ath', 'mtr', 'osa', 'ppa', 'ptri', 'hvu', 'sly', 'vvi', 'dme', 'cel', 'cne', 'ncr', 'mor', 'cci', 'ani', 'apl', 'bta', 'clu', 'dre', 'eca', 'ggo', 'hsa', 'mdo', 'mmu', 'oan', 'ocu', 'sha', 'ssa', 'ssc', 'xtr', 'metazoa', 'plantae', 'fungi')
directory_read = './data/seqs_GT/'
directory_write = './data/freqs/'
save = T#Save data or compare to existing?

for(species in species_list){
  logo = readRDS(paste(directory_read, 'logo_', species, '.rds', sep=''))
  logo = logo[which(apply(MARGIN=1, logo, FUN=function(x){return(sum(is.na(x) | x == ''))}) == 0),]#Remove nonsense rows.
  logo = logo[,1:9]#Keep needed rows.
  logo = logo[which(logo[,4] == 'G' & logo[,5] == 'T'),]
  fi = apply(logo, MARGIN=2, FUN=function(x){return(table(factor(x, levels=c('A', 'C', 'G', 'T')))/nrow(logo))})
  fij =  matrix(data=0, nrow=36, ncol=36)
  for(site_1 in 1:8){
    for(site_2 in (site_1+1):9){
      for(amino_1 in c('A', 'C', 'G', 'T')){
        for(amino_2 in c('A', 'C', 'G', 'T')){
          i_1 = 4 * (site_1 - 1) + which(c('A', 'C', 'G', 'T') == amino_1)
          i_2 = 4 * (site_2 - 1) + which(c('A', 'C', 'G', 'T') == amino_2)
          fij[i_1, i_2] = sum(logo[,site_1] == amino_1 & logo[,site_2] == amino_2)/nrow(logo)
          fij[i_2, i_1] = fij[i_1, i_2]
        }
      }
    }
  }
  
  if(save){#If I choose to save, do so.
    write.table(t(fi), file=paste(directory_write, species, '_gt/fi.txt', sep=''), row.names=F, col.names=F, quote=F)
    write.table(fij, file=paste(directory_write, species, '_gt/fij.txt', sep=''), row.names=F, col.names=F, quote=F)
  } else{#Else check if currently saved versions match.
    fi_now = read.table(paste(directory_write, species, '_gt/fi.txt', sep=''), header=F)
    fij_now = read.table(paste(directory_write, species, '_gt/fij.txt', sep=''), header=F)
    
    print(species)
    print(length(which(fi_now != t(fi), arr.ind=T)) != 0)
    print(length(which(fij_now != fij, arr.ind=T)) != 0)
  }
}
