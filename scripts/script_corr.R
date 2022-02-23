rm(list=ls())

folders = c('freqs', 'freqs_GT')
df_corrs = data.frame()
df_corrs_GT = data.frame()
threshhold = 0.05

for(folder in folders){
  data_dir = paste('./data/', folder, '/', sep='')
  sp_list = list.files(data_dir)
  sp_list = sp_list[which(nchar(sp_list) == 3)]
  
  sp_exclude = c("cbr", "dps", "dre", "dsi", "dya", "gga", "oan", "ptr", "zma", "hvu", "vra", "tae", 'ssa', 'xtr')
  sp_list = sp_list[which(!sp_list %in% sp_exclude)]
  
  plants_list = c('ath', 'hvu', 'mtr', 'osa', 'sly', 'ptri', 'ppa', 'vvi')
  animals_list = c('apl', 'bta', 'cel', 'cbr', 'clu', 'dme', 'dps', 'dre', 'dsi', 'dya', 'eca', 'ggo', 'hsa', 'mdo', 'mmu', 'ocu', 'sha', 'ssa', 'ssc', 'xtr')
  fungi_list = c('ani', 'cci', 'cne', 'mor', 'ncr')
  
  fi_list = lapply(paste(data_dir, sp_list, '/fi.txt', sep=''), FUN=read.table)
  fij_list = lapply(paste(data_dir, sp_list, '/fij.txt', sep=''), FUN=read.table)
  names(fi_list) = sp_list
  names(fij_list) = sp_list
  
  fij_list = lapply(fij_list, FUN=function(x){
    diag(x) = rowSums(x)/8
    return(x)
  })
  
  fi_list = lapply(fi_list, FUN=as.matrix)
  fij_list = lapply(fij_list, FUN=as.matrix)
  
  fcorr_list = lapply(1:length(sp_list), FUN=function(i){
    finum = as.numeric(t(fi_list[[i]]))
    fimat = finum %o% finum
    for(j in 1:9){
      inds = (4*(j-1)+1):(4*j)
      fimat[inds,inds] = 0
    }
    diag(fimat) = finum ** 2
    return(fij_list[[i]] - fimat)
  })
  
  sites_names = c(-3:-1,1:6)
  letter_names = c('A', 'C', 'G', 'T')
  
  corr_plotready = lapply(fcorr_list, FUN=function(x){
    ready = numeric()
    k = 1
    for(i in 1:8){
      for(j in (i+1):9){
        for(is in 1:4){
          for(js in 1:4){
            ready = c(ready, x[(4*(i-1)+is),(4*(j-1)+js)])
            names(ready)[k] = paste(sites_names[i], ':', letter_names[is], '_', sites_names[j], ':', letter_names[js], sep='')
            k = k + 1
          }
        }
      }
    }
    return(ready)
  })
  names(corr_plotready) = sp_list
  
  if(folder == 'freqs_GT'){
    corr_plotready_GT = corr_plotready
  } else{
    corr_plotready_all = corr_plotready
  }
  rm(corr_plotready)
}
 
 
equals_zero = abs(corr_plotready_all[[1]]) < threshhold & abs(corr_plotready_GT[[1]]) < threshhold
for(i in 2:length(corr_plotready_all)){
  equals_zero = equals_zero & ((abs(corr_plotready_all[[i]]) < threshhold) & (abs(corr_plotready_GT[[i]]) < threshhold))
}

for(folder in folders){
  corr_plotready = list()
  if(folder == 'freqs_GT'){
    corr_plotready = corr_plotready_GT
  } else{
    corr_plotready = corr_plotready_all
  }
  
  nonzero_corrs = lapply(corr_plotready, FUN=function(x){return(x[which(!equals_zero)])})
  names(nonzero_corrs) = sp_list
  
  #colors_sp = ifelse(sp_list %in% plants_list, 'green', ifelse(sp_list %in% animals_list, 'red', ifelse(sp_list %in% fungi_list, 'blue', 'black')))
  
  #mincorr = min(unlist(nonzero_corrs))
  #maxcorr = max(unlist(nonzero_corrs))
  #
  #plot(nonzero_corrs[[1]], type='p', xaxt = "n", xlab=NA, ylab='Correlation', col=colors_sp[1], ylim=c(mincorr, maxcorr))
  #axis(1, at=1:length(names(nonzero_corrs[[1]])), labels=names(nonzero_corrs[[1]]), las=2)
  #for(i in 2:length(nonzero_corrs)){
  #  points(nonzero_corrs[[i]], type='p', col=colors_sp[i])
  #}
  #legend(1, 0, pch=1, col=c('green', 'red', 'blue'), legend=c('Plantae', 'Metazoa', 'Fungi'))
  
  species_vec = sapply(names(unlist(nonzero_corrs)), FUN=function(x){return(strsplit(x, split='\\.')[[1]])})[1,]
  site_vec = sapply(names(unlist(nonzero_corrs)), FUN=function(x){return(strsplit(x, split='\\.')[[1]])})[2,]
  group_vec = ifelse(species_vec %in% plants_list, 'plantae', ifelse(species_vec %in% animals_list, 'metazoa', ifelse(species_vec %in% fungi_list, 'fungi', 'other')))
  
  if(folder == 'freqs_GT'){
    df_corrs_GT = data.frame(INT=site_vec, VAL=unlist(nonzero_corrs), SP=species_vec, GRP=group_vec)
  } else{
    df_corrs_all = data.frame(INT=site_vec, VAL=unlist(nonzero_corrs), SP=species_vec, GRP=group_vec)
  }
}

df_corrs_all = df_corrs_all[which(!df_corrs_all$INT == '1:G_2:T'),]#Remove trivial interaction.
df_corrs_GT = df_corrs_GT[which(!df_corrs_GT$INT == '1:G_2:T'),]#Remove trivial interaction.

df_corrs_all$OGT = 'All sequences'
df_corrs_GT$OGT = 'GT only'
df_corrs_full = rbind(df_corrs_all, df_corrs_GT)

library(ggplot2)

#Plotting all together.
#plot <- ggplot(df_corrs, aes(x=factor(INT), y=VAL, color=factor(GRP), fill=factor(GRP)))+
#  geom_boxplot()+
#  theme(legend.position='bottom', legend.title=element_blank())+
#  scale_color_manual(breaks = c('metazoa', 'fungi', 'plantae'), values=c('darkred', 'darkblue', 'darkgreen'))+
#  scale_fill_manual(breaks = c('metazoa', 'fungi', 'plantae'), values=c('pink', 'cyan', 'lightgreen'))+
#  xlab('Site interaction')+
#  ylab('Correlation value')

#plot

#Plot in facet grid.
WIDTH = 5

p = ggplot(df_corrs_full, aes(x=factor(OGT), y=VAL, color=factor(GRP), fill=factor(GRP)))+
    geom_boxplot()+
    theme(legend.position='bottom', legend.title=element_blank())+
    scale_color_manual(breaks = c('metazoa', 'fungi', 'plantae'), values=c('darkred', 'darkblue', 'darkgreen'))+
    scale_fill_manual(breaks = c('metazoa', 'fungi', 'plantae'), values=c('pink', 'cyan', 'lightgreen'))+
    xlab('')+
    ylab('Correlation value')+
    facet_wrap(vars(INT), ncol=WIDTH)+
    geom_point(color='black', position=position_jitterdodge(), size=0.8)
p

#Test de Kolmogorov-Smirnov para ver si correlaciones vienen de misma distribucion.
library(data.table)
dt_corrs = as.data.table(cbind(df_corrs_all[,which(!colnames(df_corrs_all) %in% c('OGT', 'SP'))], VALGT=df_corrs_GT$VAL[match(rownames(df_corrs_GT), rownames(df_corrs_all))]))
pvs = as.data.frame(dt_corrs[, KS:=wilcox.test(VAL, VALGT, paired=T, exact=T)$p.value, by=list(GRP, INT)])
pvs = pvs[,which(!colnames(pvs) %in% c('VAL', 'VALGT'))]
pvs = pvs[!duplicated(pvs),]
pvs$KS = p.adjust(pvs$KS)

#Ploteamos cuales puntos se extravian.
plot(dt_corrs$VAL, dt_corrs$VALGT, col=c('red', 'green', 'blue')[match(dt_corrs$GRP, c('metazoa', 'plantae', 'fungi'))])
abline(0, 1)
legend(0.02, -0.02, legend=c('Metazoa', 'Plantae', 'Fungi'), col=c('red', 'green', 'blue'), pch=1)
