#sp_list = c('cne', 'ath', 'mtr', 'osa', 'ptr', 'ggo', 'hsa', 'cle', 'cbr', 'aga', 'dme', 'dps', 'dsi', 'dya')
#sp_full_list = c('cryptococcus', 'arabidopsis', 'medicago', 'oryza', 'pan', 'gorilla', 'humano', 'celegans', 'briggsae', 'gambiae', 'drosophila', 'pseudoobscura', 'simulans', 'yacuba')

#sp_list = c('hsa_lowenergy', 'hsa_nonu12')
#sp_full_list = c('humano_lowenergy', 'humano_u12')

#sp_list = 'hsa_nonu12'
#sp_full_list = 'humano_u12_invertida'

sp_list      = c('ath', 'mtr', 'osa', 
                 'ggo', 'hsa', 
                 'dre','mmu',
                 'dme',
                 'cel',
                 'cne', 'apl', 'bta', 'clu', 'eca', 'mdo', 'oan', 'ocu', 'sha', 'ssa', 'hvu', 'sly', 'ssc', 'xtr', 'ani', 'cci', 'mor', 'ncr', 'ppa', 'ptri', 'vvi')
#                 'pan', 'gorilla', 'humano',
#                 'zebrafish','mouse',
#                 'drosophila','pseudoobscura', 'yacuba', 
#                 'celegans', 'briggsae',
#                 'cryptococcus')

plants = c('ath', 'mtr', 'osa', 'hvu', 'ppa', 'ptri', 'sly', 'vvi')
fungi = c('cne', 'ani', 'ncr', 'mor', 'cci')
vertebrates = c('ggo', 'hsa', 'dre', 'mmu', 'bta', 'clu', 'xtr', 'ssc', 'mdo', 'oan', 'ocu', 'sha', 'apl', 'ssa', 'eca')
invertebrates = c('dme', 'cel')

gamma <- '0.025000'
N     <- 200  #top colors


# Leo tabloas de gamma especificado ----
paramJ_list <- paramH_list <- Pi_list <- Pij_list <- fi_list <- fij_list <- list()
for(isp in seq_along(sp_list)){
  if(length(gamma)>1){
    ggamma <- gamma[isp]
  }else{
    ggamma <- gamma
  }
  ppath <- paste0('./runs/',sp_list[isp],'/gamma',format(ggamma))
  d    <- dir(ppath,pattern = 'ParamH')
  irun <- max(as.numeric(sub('ParamH','',d)))
  
  jij<- read.table(paste0(ppath,'/ParamJ',irun)) 
  hi <- read.table(paste0(ppath,'/ParamH',irun)) 
  colnames(hi) <- c("A","C","G","T")
  rownames(hi) <- c(-3:-1,1:6)
  colnames(jij)<- rownames(jij) <- apply(expand.grid(colnames(hi),rownames(hi))[,2:1],1,paste,collapse=":")
  
  paramJ_list[[sp_list[isp]]] <- jij
  paramH_list[[sp_list[isp]]] <- hi
  
  
  # . Eleccion de Gauge ----
  # 'A general pairwise interaction model provides accurate description
  # of in Vivo transcription factor binding sites´ Santolini2014 PLoSOne
  if(FALSE){
    hhi <- t(apply(hi,1,function(x){x-mean(x)}))
    aux <-limma::strsplit2(colnames(jij),':')
    pos <- aux[,1]
    base<- aux[,2]
    jjij <- as.matrix(jij)
    for(ipos in 1:length(unique(pos))){
      for(jpos in 1:length(unique(pos))){
        if(jpos==ipos) next
        jaux <- jjij[pos==unique(pos)[ipos],pos==unique(pos)[jpos]]
        jaux <- jaux + mean(jaux) - t(matrix(rep(apply(jaux,1,mean),4),ncol=4))- matrix(rep(apply(jaux,2,mean),4),ncol=4)
        
        jjij[pos==unique(pos)[ipos],pos==unique(pos)[jpos]] <- jaux
        
      }
    }
  }
  
  Pij_list[[sp_list[isp]]] <-read.table(paste0(ppath,'/Pij',irun)) 
  Pi_list[[sp_list[isp]]] <-read.table(paste0(ppath,'/Pi',irun)) 
  
  fi_list[[sp_list[isp]]] <-read.table(paste('./data/freqs/',sp_list[isp], '/fi.txt', sep=''))
  fij_list[[sp_list[isp]]] <-read.table(paste('./data/freqs/',sp_list[isp], '/fij.txt', sep=''))
  
  
}

paramJ_list = lapply(FUN=as.matrix, paramJ_list)
paramJ_list_normalized = lapply(paramJ_list, FUN=function(x){return(as.matrix(x/sum(abs(x))))})
list_of_Js = paramJ_list_normalized

get.pijmpi = function(pij_mat, pi_mat){
  pi_array = t(pi_mat)[1:36]
  return(pij_mat - pi_array %o% pi_array)
}

heatmap.image = function(xlist,metric=c('inner','euclidean','cos','cosB'),
                         type=c('dendro','dendrot','heatmap')[1],
                         main='',horiz=FALSE,bexp=FALSE){
  mat_J_vecs = do.call(rbind, lapply(FUN=function(x){return(as.numeric(t(as.matrix(x))))}, xlist))
  mat_J_vecs[!is.finite(mat_J_vecs)]<- -10
  if(bexp) mat_J_vecs <- exp(mat_J_vecs)
  if(FALSE){
    n<-nrow(mat_J_vecs)
    par(mar=c(5,4,2,2))
    matplot(t(mat_J_vecs),col=1:n,lty=1:n,typ='b',pch=20)
    abline(v=0.5+seq(4,32,4),lty=2,col='gray')
    legend('topleft',rownames(mat_J_vecs),col=1:n,lty=1:n,cex=0.6,inset=0.02,bty='n')
  }
  if(metric=='cos'){
    similarity = mat_J_vecs %*% t(mat_J_vecs)
    similarity <- similarity/sqrt(diag(similarity)%*%t(diag(similarity)))
    distancias = 1-similarity
  }
  if(metric=='cosB'){
    similarity = mat_J_vecs %*% t(mat_J_vecs)
    similarity <- similarity/sqrt(diag(similarity)%*%t(diag(similarity)))
    distancias = 1/similarity
  }
  if(metric=='inner'){
    similarity = mat_J_vecs %*% t(mat_J_vecs)
    distancias = 1/(similarity)
  }
  if(metric=='euclidean'){
    distancias = as.matrix(dist(mat_J_vecs, upper=T, method='euclidean'))    
  }

  colnames(distancias) = sp_list
  rownames(distancias) = sp_list
  #plot(heatmap(distancias, symm=T, scale='none'))
  bplotted <- FALSE
  if(type=='heatmap'){
   heatmap(distancias, symm=F, scale='none',main=main)
   bplotted <- TRUE
  }
  if(type=='dendro'){
    dendro1<-as.dendrogram(hclust(as.dist(distancias),method='complete'))
    plot(dendro1,main=main,horiz=horiz)
    bplotted <- TRUE
  }
  if(type=='dendrot'){
    dendro1<-as.dendrogram(hclust(as.dist(distancias),method='complete'))
    plot(dendro1, nodePar = list(pch = c(1,NA), cex = 0.8, lab.cex = 0.8),
         type = "t", center = TRUE,main=main,horiz=horiz)
    bplotted <- TRUE
  }
  if(!bplotted) cat('type argument should be one of: dendro, dendrot, heatmap \n')
}

im     <- 2
ttype  <- c('dendro','dendrot','heatmap')[1]
horiz  <- TRUE
bexp   <- F
heatmap.image(fi_list, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,horiz=TRUE,main='f_i')
heatmap.image(fij_list,metric=c('inner','euclidean','cos','cosB')[im],type=ttype,horiz=TRUE,main='f_ij')

heatmap.image(Pi_list, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,horiz=TRUE,main='P_i')
heatmap.image(Pij_list, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,horiz=TRUE,main='P_ij')

heatmap.image(paramH_list, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,bexp=bexp,main='H_i')
heatmap.image(paramJ_list, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,bexp=bexp,main='J_ij')

#Con J individuales.
J_nonulos = lapply(paramJ_list, FUN=function(x){return(which(x != 0, arr.ind=T))})
J_nonulos = unique(do.call('rbind', J_nonulos))
rownames(J_nonulos) = character()
J_nonulos[,1] = colnames(paramJ_list[[1]])[as.numeric(J_nonulos[,1])]
J_nonulos[,2] = colnames(paramJ_list[[1]])[as.numeric(J_nonulos[,2])]
J_nonulos = as.data.frame(J_nonulos, stringsAsFactors=F)

pval_list = numeric()
for(i in 1:nrow(J_nonulos)){
  site_1 = J_nonulos[i,1]
  site_2 = J_nonulos[i,2]
  J_particular = unlist(lapply(paramJ_list, FUN=function(x){return(x[site_1,site_2])}))
  names(J_particular) = names(paramJ_list)
  J_particular = sort(J_particular)

  J_particular_filtered = J_particular[which(!names(J_particular) %in% c('dps', 'dya', 'cbr', 'ptr', 'ggo'))]
  con_table = matrix(data=0, ncol=2, nrow=3)
  rownames(con_table) = c('Plants', 'Fungi', 'Metazoa')
  colnames(con_table) = c('Present', 'Not Present')
  
  con_table[1,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% plants)] != 0)
  con_table[2,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% fungi)] != 0)
  con_table[3,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% c(invertebrates, vertebrates))] != 0)
  con_table[1,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% plants)] == 0)
  con_table[2,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% fungi)] == 0)
  con_table[3,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% c(invertebrates, vertebrates))] == 0)
  pval_list = c(pval_list, fisher.test(con_table)$p.value)
}

J_nonulos = cbind(J_nonulos, PV=pval_list)
J_nonulos = J_nonulos[order(J_nonulos$PV),]
J_nonulos = J_nonulos[1+2*(1:(nrow(J_nonulos)/2)-1),]
J_nonulos$PV = p.adjust(J_nonulos$PV, method='fdr')

#Armo matriz de Js con especies.
mat_J_esp = matrix(nrow=nrow(J_nonulos), ncol=length(paramJ_list), data=0)
rownames(mat_J_esp) = apply(MARGIN=1, J_nonulos[,1:2], FUN=paste0, collapse='::')
colnames(mat_J_esp) = names(paramJ_list)
for(i in 1:nrow(mat_J_esp)){
  for(j in 1:ncol(mat_J_esp)){
    sp_name = colnames(mat_J_esp)[j]
    row_j = J_nonulos$row[i]
    col_j = J_nonulos$col[i]
    mat_J_esp[i,j] = paramJ_list[[which(names(paramJ_list) == sp_name)]][row_j,col_j]
  }
}
write.table(mat_J_esp, file='species_j_mat.csv', sep=',', row.names=T, col.names=T, quote=F)
write.table(J_nonulos, file='j_nonulos.csv', sep=',', row.names=F, col.names=T, quote=F)

#Para ver uno a la vez.
site_1_list = J_nonulos$col
site_2_list = J_nonulos$row
presence_df = matrix(nrow=length(site_1_list), ncol=3)
colnames(presence_df) = c('plants', 'animals', 'fungi')
for(i in 1:length(site_1_list)){
  J_particular = unlist(lapply(paramJ_list, FUN=function(x){return(x[site_1_list[i],site_2_list[i]])}))
  present_in = names(J_particular)[which(J_particular != 0)]
  presence_df[i,] = c(sum(present_in %in% plants), sum(present_in %in% c(vertebrates, invertebrates)), sum(present_in %in% fungi))
}
sig_ints = as.data.table(cbind(J_nonulos[1:length(site_1_list),-which(colnames(J_nonulos) == 'PV')], presence_df, PV=J_nonulos$PV[1:length(site_1_list)]))

#heatmap.image(J_particular, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,bexp=bexp)
#plot(J_particular, col=ifelse(names(J_particular) %in% plants, 'green', ifelse(names(J_particular) %in% fungi, 'red', ifelse(names(J_particular) %in% invertebrates, 'blue', 'black'))))
#legend(20, -1, legend=c('Plants', 'Fungi', 'Invertebrates', 'Vertebrates'), col=c('green', 'red', 'blue', 'black'), pch=1)

J_particular_filtered = J_particular[which(!names(J_particular) %in% c('dps', 'dya', 'cbr', 'ptr'))]
con_table = matrix(data=0, ncol=2, nrow=3)
rownames(con_table) = c('Plants', 'Fungi', 'Metazoa')
colnames(con_table) = c('Present', 'Not Present')

con_table[1,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% plants)] != 0)
con_table[2,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% fungi)] != 0)
con_table[3,1] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% c(invertebrates, vertebrates))] != 0)
con_table[1,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% plants)] == 0)
con_table[2,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% fungi)] == 0)
con_table[3,2] = sum(J_particular_filtered[which(names(J_particular_filtered) %in% c(invertebrates, vertebrates))] == 0)

orgs_order = c(plants, fungi, invertebrates, vertebrates)
ind_js = unlist(lapply(paramJ_list, FUN=function(x){return(x[site_1,site_2])}))
ind_js = ind_js[match(orgs_order, names(ind_js))]

#El grafico vector presencia.
sps_ach = c('sly', 'ptri', 'mtr', 'ath', 'vvi', 'hvu', 'osa', 'ppa', 'cel', 'dme', 'oan', 'hsa', 'ggo', 'mmu', 'ocu', 'ssc', 'bta', 'eca', 'clu', 'sha', 'mdo', 'apl', 'xtr', 'dre', 'ssa', 'ani', 'ncr', 'mor', 'cci', 'cne')
presence_mat = matrix(nrow=length(sps_ach), ncol=length(site_1_list))
rownames(presence_mat) = sps_ach
colnames(presence_mat) = paste(gsub(site_2_list, pattern=':', replacement=''), gsub(site_1_list, pattern=':', replacement=''), sep=':')
for(i in 1:length(site_1_list)){
  J_particular = unlist(lapply(paramJ_list, FUN=function(x){return(x[site_1_list[i],site_2_list[i]])}))
  J_particular = ifelse(J_particular == 0, 0, ifelse(J_particular > 0, '+', '-'))
  presence_mat[,i] = J_particular[match(rownames(presence_mat), names(J_particular))]
}
write.table(presence_mat, file='presence_mat.csv', sep=',', quote=F, row.names=T, col.names=T)

#El grafico vector presencia.
#sps_ach = c('sly', 'ptri', 'mtr', 'ath', 'vvi', 'hvu', 'osa', 'ppa', 'cel', 'dme', 'oan', 'hsa', 'ggo', 'mmu', 'ocu', 'ssc', 'bta', 'eca', 'clu', 'sha', 'mdo', 'apl', 'xtr', 'dre', 'ssa', 'ani', 'ncr', 'mor', 'cci', 'cne')
js_mat = matrix(nrow=length(sps_ach), ncol=length(site_1_list))
rownames(js_mat) = sps_ach
colnames(js_mat) = paste(gsub(site_2_list, pattern=':', replacement=''), gsub(site_1_list, pattern=':', replacement=''), sep=':')
for(i in 1:length(site_1_list)){
  J_particular = unlist(lapply(paramJ_list, FUN=function(x){return(x[site_1_list[i],site_2_list[i]])}))
  js_mat[,i] = J_particular[match(rownames(presence_mat), names(J_particular))]
}
write.table(t(js_mat), file='presence_mat.csv', sep=',', quote=F, row.names=T, col.names=T)


if(F){
  #Producto interno de Js para hacer el dendrograma.
  scaling_vector = 1 - J_nonulos$PV
  scaling_vector = scaling_vector / sqrt(sum(scaling_vector ^ 2))
  paramJ_list_adjusted = lapply(paramJ_list, FUN=function(x){
    result = x
    for(i in 1:nrow(J_nonulos)){
      result[J_nonulos$row[i], J_nonulos$col[i]] = x[J_nonulos$row[i], J_nonulos$col[i]] * sqrt(scaling_vector[i])
      result[J_nonulos$col[i], J_nonulos$row[i]] = result[J_nonulos$row[i], J_nonulos$col[i]]
    }
    return(result)
  })
  
  heatmap.image(paramJ_list_adjusted, metric=c('inner','euclidean','cos','cosB')[im],type=ttype,bexp=bexp)
  
  mat_similaridades = matrix(0, nrow=length(paramJ_list), ncol=length(paramJ_list))
  colnames(mat_similaridades) = names(paramJ_list)
  rownames(mat_similaridades) = names(paramJ_list)
  for(i in 1:nrow(mat_similaridades)){
    for(j in 1:i){
      mat_similaridades[i,j] = max(paramJ_list_adjusted[[i]] - paramJ_list_adjusted[[j]])
      mat_similaridades[j,i] = mat_similaridades[i,j]
    }
  }
  plot(as.dendrogram(hclust(as.dist(mat_similaridades)), method='complete'))
}