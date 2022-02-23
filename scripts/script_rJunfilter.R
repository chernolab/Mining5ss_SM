data_dir = './data/seqs/'

hist.loglog = function(data, num_pts=10, fit=T, xlab_fun='', ylab_fun='Density', title_fun=''){
  breaks_log = min(data) + (max(data) - min(data)) * c(0, 0.5**(0:num_pts))
  histo = hist(data, breaks=breaks_log, plot=F)
  plot(histo$mids, histo$density, log='xy', xlab=xlab_fun, ylab=ylab_fun)
  if(fit==T){
    non_zero = which(histo$density > 0)
    fit_plaw = lm(log(histo$density[non_zero]) ~ log(histo$mids[non_zero]))
    lines(histo$mids[non_zero], exp(fit_plaw$coefficients[1]) * (histo$mids[non_zero] ** fit_plaw$coefficients[2]), col='red')
    title(title_fun)
  }
}

hist.loglog(rjunBaseCoor$Normal.median, xlab_fun='Median activity', title_fun='Normal')
hist.loglog(rjunBaseCoor$Tumor.median, xlab_fun='Median activity', title_fun='Tumor')
hist.loglog(rjunBaseCoor$NT.median, xlab_fun='Median activity', title_fun='NT')

sep.sequences = function(x){
  result = as.data.frame(do.call(rbind, strsplit(x$seq5ss, split='')))
  rownames(result) = x$JunctionID
  colnames(result) = c(-3:-1,1:6)
  return(result)
}

seqs_gr0 = sep.sequences(rjunBaseCoor[which(rjunBaseCoor$Normal.median > 0),])
seqs_gr5 = sep.sequences(rjunBaseCoor[which(rjunBaseCoor$Normal.median > 5),])

saveRDS(seqs_gr0, file=paste(data_dir, 'logo_seqs_gr0.rds', sep=''))
saveRDS(seqs_gr5, file=paste(data_dir, 'logo_seqs_gr5.rds', sep=''))
