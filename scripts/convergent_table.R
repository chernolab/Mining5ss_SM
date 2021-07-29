rm(list = ls())
require(grDevices)

n_sites = 9#Number of sites within sequences analysed.
logos_directory = 'data/'#Folder where observed frequency data will be found.
runs_directory = 'runs/'#Folder where the fitting results will be found.
gamma = '0.025000'#Gamma result selected to analyse.

#The following lists show the available organisms utilized. The script will run including all of them.
sp_list = c('cne', 'ath', 'mtr', 'osa', 'ptr', 'ggo', 'hsa', 'dre', 'mmu', 'cel', 'cbr', 'dme', 'dps', 'dya')
sp_full = c('cryptococcus', 'arabidopsis', 'medicago', 'oryza', 'pan', 'gorilla', 'humano', 'zebrafish', 'mouse', 'celegans', 'briggsae', 'drosophila', 'pseudoobscura', 'yacuba')

load.logo = function(sp){
    #Cargo los datos del genoma
    i = match(sp, sp_list)
    load(paste(logos_directory, 'seqs/logo_', sp, '.RData', sep=''))
    if(sp_full[i] == 'humano'){#En el caso de homo sapiens, la variable no se llama logo.
        logo = logo_hsa
    }
    if(sp_full[i] == 'arabidopsis'){
        logo <- logo
    }
    if(sp_full[i] == 'celegans'){
        logo = logo_cle
    }
    if(sp_full[i] == 'briggsae'){
        logo = logo_cbr
    }
    if(sp_full[i] == 'cryptococcus'){
        logo <- logo_cne
    }
    if(sp_full[i] == 'drosophila'){
        logo = logo_dme
    }
    if(sp_full[i] == 'medicago'){
        logo <- logo_mtr
    }
    if(sp_full[i] == 'oryza'){
        logo = logo_osa
    }
    if(sp_full[i] == 'pan'){
        logo = logo_ptr
    }
    if(sp_full[i] == 'gambiae'){
        logo = logo_aga
    }
    if(sp_full[i] == 'pseudoobscura'){
        logo <- logo_dps
    }
    if(sp_full[i] == 'simulans'){
        logo = logo_dsi
    }
    if(sp_full[i] == 'yacuba'){
        #logo = logo#El logo ya es correcto en yacuba, se llama logo.
    }
    if(sp_full[i] == 'gorilla'){
        logo = logo_ggo
    }
    if(sp_full[i] == 'mouse'){
      logo = logo_org
    }
    if(sp_full[i] == 'zebrafish'){
      #El logo ya es correctp en dre.
    }
    return(logo)
}

fi.distribution = function(sp){
    logo = load.logo(sp)
    #estimo fi_exp y fij_exp
    auxt<-apply(logo,2,table)
    if(class(auxt)=="list"){
        laux<-lapply(auxt,function(x){
            res<-rep(0,4)
            names(res)<-c("A","C","G","T")
            res[names(x)]<-x
            return(res)
        })
        auxt<-t(matrix(unlist(laux),byrow=TRUE,ncol=4))
        rownames(auxt)<-c("A","C","G","T")
    }
    fi <- t(apply(auxt,2,function(x){x/sum(x)}))
    return(fi)
}

fi.distribution.2 = function(logo){
    #estimo fi_exp y fij_exp
    auxt<-apply(logo,2,table)
    if(class(auxt)=="list"){
        laux<-lapply(auxt,function(x){
            res<-rep(0,4)
            names(res)<-c("A","C","G","T")
            res[names(x)]<-x
            return(res)
        })
        auxt<-t(matrix(unlist(laux),byrow=TRUE,ncol=4))
        rownames(auxt)<-c("A","C","G","T")
    }
    fi <- t(apply(auxt,2,function(x){x/sum(x)}))
    return(fi)
}

fij.distribution = function(sp){
    #Cargo los datos del genoma
    logo = load.logo(sp)
    fij <- matrix(0,ncol=36,nrow=36)                  #frecuencias de coaparicion i-letra_i j-letra_j
    for(i in 1:9){
        for(j in (i+1):9){
            if(j>9) break
            t<-table(logo[,i],logo[,j])  
            if(any(dim(t)!=c(4,4))){
                t0 <- matrix(0,ncol=4,nrow=4)
                colnames(t0)<-rownames(t0)<-c("A","C","G","T")
                t0[rownames(t),colnames(t)]<-t
                t<-t0
            }
            pij <- t/nrow(logo)
            fij[(i-1)*4+(1:4),(j-1)*4+(1:4)]<-pij

            x<-table(logo[,i])/nrow(logo)
            y<-table(logo[,j])/nrow(logo)
            pipj <- x%*%t(y)
        }
    } 
    fij<-fij+t(fij)
    colnames(fij)<-rownames(fij)<-paste(rep(c(-3:-1,1:6),each=4),c("A","C","G","T"),sep=":")
    return(fij)
}

fij.distribution.2 = function(logo){
    #Cargo los datos del genoma
    fij <- matrix(0,ncol=36,nrow=36)                  #frecuencias de coaparicion i-letra_i j-letra_j
    for(i in 1:9){
        for(j in (i+1):9){
            if(j>9) break
            t<-table(logo[,i],logo[,j])  
            if(any(dim(t)!=c(4,4))){
                t0 <- matrix(0,ncol=4,nrow=4)
                colnames(t0)<-rownames(t0)<-c("A","C","G","T")
                t0[rownames(t),colnames(t)]<-t
                t<-t0
            }
            pij <- t/nrow(logo)
            fij[(i-1)*4+(1:4),(j-1)*4+(1:4)]<-pij

            x<-table(logo[,i])/nrow(logo)
            y<-table(logo[,j])/nrow(logo)
            pipj <- x%*%t(y)
        }
    } 
    fij<-fij+t(fij)
    colnames(fij)<-rownames(fij)<-paste(rep(c(-3:-1,1:6),each=4),c("A","C","G","T"),sep=":")
    return(fij)
}

#Calculos sobre fi.
fi_list = lapply(FUN=fi.distribution, sp_list)
fi_mat = matrix(data=unlist(lapply(FUN=t, fi_list)), nrow=4*n_sites, ncol=length(sp_list))
rownames(fi_mat) = paste(unlist(lapply(FUN=rep, c(-3:-1, 1:6), 4)), ':', rep(c('A', 'C', 'G', 'T'), n_sites), sep='')
colnames(fi_mat) = sp_list

ppath_list = paste(runs_directory, sp_list, '/gamma', gamma, sep='')
d = lapply(FUN=dir, ppath_list, pattern='ParamH')
irun = unlist(lapply(FUN=function(x){return(max(as.numeric(sub('ParamH','',x)))[1])}, d))#El numero de la maxima iteracion.

paramJ_list = lapply(FUN=function(x){return(read.table(x))}, paste('runs/', sp_list, '/gamma', gamma, '/ParamJ', irun, sep=''))#Las nuevas de Cherno.
#paramJ_list = lapply(FUN=read.table, paste(params_directory, 'ParamJ_', sp_list, '_ach_15', sep=''))
paramJ_list = lapply(FUN=as.matrix, paramJ_list)

intersection.intensities = function(pi, jij, sp='', consenso_lim=0.25, plot=F){#Me retorna un array con las intensidades promedio entre los 10 tipos posibles de interaccion.
    #Ahora clasifico los consensos, y no-consenso.
    unl_pi = t(pi)[1:36]

    #Ahora armo grupos.
    con_ex = which(unl_pi[1:12] > consenso_lim)
    con_in = 12 + which(unl_pi[13:length(unl_pi)] > consenso_lim)
    nocon_ex = which(unl_pi[1:12] <= consenso_lim)
    nocon_in = 12 + which(unl_pi[13:length(unl_pi)] <= consenso_lim)

    jij_used = matrix(rep(NA, nrow(jij)*ncol(jij)), nrow=nrow(jij), ncol=ncol(jij))
    for(i in 1:8){#Meto solo lo que no esta en los bloques diagonales nulos de 4x4.
        jij_used[(4*(i-1))+(1:4), (4*i+1):ncol(jij)] = jij[(4*(i-1))+(1:4), (4*i+1):ncol(jij)]
    }

    #Ahora armo todos los grupos de jijs.
    jij_con_in_con_ex = jij_used[con_ex, con_in]#Nunca tendra NAs.
    jij_con_in_nocon_ex = jij_used[nocon_ex, con_in]#Nunca tendra NAs.
    jij_con_ex_nocon_in = jij_used[con_ex, nocon_in]#Nunca tendra NAs.
    jij_nocon_in_nocon_ex = jij_used[nocon_ex, nocon_in]#Nunca tendra NAs.
    jij_con_in_con_in = jij_used[con_in, con_in]
    jij_con_in_nocon_in = jij_used[con_in, nocon_in]
    jij_nocon_in_nocon_in = jij_used[nocon_in, nocon_in]
    jij_con_ex_con_ex = jij_used[con_ex, con_ex]
    jij_con_ex_nocon_ex = jij_used[con_ex, nocon_ex]
    jij_nocon_ex_nocon_ex = jij_used[nocon_ex, nocon_ex]
    
    #A los bloques diagonales tengo que quitarles las componentes simetricas.
    jij_con_in_con_in = jij_con_in_con_in[!is.na(jij_con_in_con_in)]
    jij_con_in_nocon_in = jij_con_in_nocon_in[!is.na(jij_con_in_nocon_in)]
    jij_nocon_in_nocon_in = jij_nocon_in_nocon_in[!is.na(jij_nocon_in_nocon_in)]
    jij_con_ex_con_ex = jij_con_ex_con_ex[!is.na(jij_con_ex_con_ex)]
    jij_con_ex_nocon_ex = jij_con_ex_nocon_ex[!is.na(jij_con_ex_nocon_ex)]
    jij_nocon_ex_nocon_ex = jij_nocon_ex_nocon_ex[!is.na(jij_nocon_ex_nocon_ex)]
    
    if(plot){
        #png(paste('./logos/inter_', sp, '.png', sep=''))
    }
    grupos = c('IC-EC', 'IC-ENC', 'INC-EC', 'INC-ENC', 'IC-IC', 'IC-INC', 'INC-INC', 'EC-EC', 'EC-ENC', 'ENC-ENC')
    valores = c(mean(jij_con_in_con_ex), mean(jij_con_in_nocon_ex), mean(jij_con_ex_nocon_in), mean(jij_nocon_in_nocon_ex), mean(jij_con_in_con_in), mean(jij_con_in_nocon_in), mean(jij_nocon_in_nocon_in), mean(jij_con_ex_con_ex), mean(jij_con_ex_nocon_ex), mean(jij_nocon_ex_nocon_ex))
    if(plot){
        barplot(valores, ylab='Mean intensity', col=ifelse(valores>0, "green", "red"), names.arg=grupos, las=2, ylim=c(-0.2, 0.2))
        abline(mean(abs(jij_used[!is.na(jij_used)])), 0, col='blue', lty=2)
        abline(-mean(abs(jij_used[!is.na(jij_used)])), 0, col='blue', lty=2)
        #dev.off()
    }

    return(valores)
}

mean_values = unlist(lapply(paramJ_list, FUN=function(x){return(mean(abs(x)))}))
inter_intensities = mapply(FUN=intersection.intensities, fi_list, paramJ_list, sp_list)
interaction_names = c('IC-EC', 'IC-ENC', 'INC-EC', 'INC-ENC', 'IC-IC', 'IC-INC', 'INC-INC', 'EC-EC', 'EC-ENC', 'ENC-ENC')
colnames(inter_intensities) = sp_list
rownames(inter_intensities) = interaction_names
inter_intensities = inter_intensities / (rep(1, nrow(inter_intensities)) %o% sqrt(colSums(inter_intensities ** 2)))#Normalizo cada especie por norma 2.
inter_intensities_char = t(matrix(as.character(round(inter_intensities, 2)), nrow=nrow(inter_intensities), ncol=ncol(inter_intensities)))
colorit = F#Don't color entries by their sign.
for(i in 1:length(sp_list)){
    avg = 0
    for(j in 1:nrow(inter_intensities)){
      if(colorit){
        inter_intensities_char[i,j] = ifelse(abs(inter_intensities[j,i]) > avg, ifelse(inter_intensities[j,i] > 0, paste('\\color{green}', inter_intensities_char[i,j], sep=''), paste('\\color{red}', inter_intensities_char[i,j], sep='')),  paste('\\color{blue}', inter_intensities_char[i,j], sep=''))
      }
    }
}
rownames(inter_intensities_char) = paste('\\textbf{', colnames(inter_intensities), '}', sep='')
lines = paste0(c('', rownames(inter_intensities)), collapse='&')
for(i in 1:nrow(inter_intensities_char)){
  lines = c(lines, paste0(c(rownames(inter_intensities_char)[i], inter_intensities_char[i,]), collapse='&'))
}
#We can then output the results in a LaTeX-friendly format for substituting within a tabular environment.
write.table(file='inter_intensities.txt', lines, sep='', quote=F, col.names=F, row.names=F, eol='\\\\\n')
