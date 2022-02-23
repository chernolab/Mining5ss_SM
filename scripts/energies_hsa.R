dir_params = './'#Directorio donde esten los parametros ajustados por la corrida.
dir_seqs = './'#Directorio donde esten las secuencias de entrada.

max.run = function(sp){
  ppath = paste0('./runs/', sp,'/gamma',format(ggamma))
  d = dir(ppath,pattern = 'ParamH')
  return(max(as.numeric(sub('ParamH','',d))))
}

read.seq = function(sp){
  ppath = paste0(dir_seqs)
  d = dir(ppath, pattern=sp)
  logo = load(paste(ppath, d, sep=''))
  return(eval(as.name(logo)))
}

paramH = read.table(paste(dir_params, '/ParamH', sep=''))
paramJ = read.table(paste(dir_params, '/ParamJ', sep=''))
seqs = read.seq('logo')

energy.calc = function(seq, paramH, paramJ){#Para una sola seq.
  site_idx = 4*(0:8) + match(seq, c('A', 'C', 'G', 'T'))
  paramH_array = as.numeric(t(paramH))
  nrg_H = sum(paramH_array[site_idx])
  nrg_J = sum(paramJ[site_idx,site_idx])/2#Esto agarra todas las entradas de sitio presentes en la secuencia y suma las energias J. Es simetrica entonces divido por dos.
  return(nrg_H + nrg_J)
}

#Finalmente las energias. Una lista con una entrada por especie, cada entrada es un array con las energias de cada secuencia.
energies = apply(seqs, MARGIN=1, FUN=energy.calc, paramH=paramH, paramJ=paramJ)
