dir_params = './runs/'#Directorio donde esten los parametros ajustados por la corrida.
dir_seqs = './data/seqs/'#Directorio donde esten las secuencias de entrada.
species = c('hsa')#Especies.
gamma = '0.025000'#Gamma elegido.

max.run = function(sp){
  ppath = paste0('./runs/', sp,'/gamma',format(gamma))
  d = dir(ppath,pattern = 'ParamH')
  return(max(as.numeric(sub('ParamH','',d))))
}

read.seq = function(sp){
  ppath = paste0('./data/seqs/')
  d = dir(ppath, pattern=sp)
  d = d[grep(d, pattern='\\.rds')]
  logo = readRDS(paste(ppath, d, sep=''))
  return(logo)
}

paramH_list = lapply(species, FUN=function(sp){return(read.table(paste(dir_params, sp, '/gamma', gamma, '/ParamH', max.run(sp), sep='')))})
paramJ_list = lapply(species, FUN=function(sp){return(read.table(paste(dir_params, sp, '/gamma', gamma, '/ParamJ', max.run(sp), sep='')))})
seq_list = lapply(species, FUN=read.seq)

energy.calc = function(seq, paramH, paramJ){#Para una sola seq.
  site_idx = 4*(0:8) + match(seq, c('A', 'C', 'G', 'T'))
  paramH_array = as.numeric(t(paramH))
  nrg_H = sum(paramH_array[site_idx])
  nrg_J = sum(paramJ[site_idx,site_idx])/2#Esto agarra todas las entradas de sitio presentes en la secuencia y suma las energias J. Es simetrica entonces divido por dos.
  return(nrg_H + nrg_J)
}

#Finalmente las energias. Una lista con una entrada por especie, cada entrada es un array con las energias de cada secuencia.
energies_list = lapply(1:length(species), FUN=function(i){return(apply(seq_list[[i]], MARGIN=1, FUN=energy.calc, paramH=paramH_list[[i]], paramJ=paramJ_list[[i]]))})
