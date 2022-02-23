# Organisms to include in the analysis
sp_list = c('seqs_gr0', 'seqs_gr5')#c('plantae', 'metazoa')#c('ath', 'mtr', 'osa', 'ppa', 'ptri', 'hvu', 'sly', 'vvi',
            #'dme',
            #'cel',
            #'cne', 'ncr', 'mor', 'cci', 'ani',
            #'apl', 'bta', 'clu', 'dre', 'eca', 'ggo', 'hsa', 'mdo', 'mmu', 'oan', 'ocu', 'sha', 'ssa', 'ssc', 'xtr',
            #'metazoa', 'plantae', 'fungi')

gamma = '0.025000'         #Gamma result selected to analyse.
runs_directory = 'runs/'   #Folder where the fitting results will be found.



N = 100  #Number of top J params to include within the circos diagram (to avoid clutter and focus on the important ones).

#----------------------------------------------------------------------------#

ppath_list = paste(runs_directory, sp_list, '/gamma', gamma, sep='')
d = lapply(FUN=dir, ppath_list, pattern='ParamH')
irun = unlist(lapply(FUN=function(x){return(max(as.numeric(sub('ParamH','',x)))[1])}, d)) #El numero de la maxima iteracion.

paramJ_list = lapply(FUN=function(x){return(read.table(x))}, paste('runs/', sp_list, '/gamma', gamma, '/ParamJ', irun, sep=''))#Las nuevas de Cherno.
Pij_list = lapply(FUN=function(x){return(read.table(x))}, paste('runs/', sp_list, '/gamma', gamma, '/Pij', irun, sep=''))#Las nuevas de Cherno.
Pi_list = lapply(FUN=function(x){return(read.table(x))}, paste('runs/', sp_list, '/gamma', gamma, '/Pi', irun, sep=''))#Las nuevas de Cherno.
paramJ_list = lapply(FUN=as.matrix, paramJ_list)
paramJ_list_normalized = lapply(paramJ_list, FUN=function(x){return(as.matrix(x/sum(abs(x))))})
list_of_Js = paramJ_list_normalized

fi_list = lapply(FUN=function(x){return(read.table(paste('data/freqs/', x, '/fi.txt', sep='')))}, sp_list)

list_of_Js = paramJ_list

for(species in 1:length(sp_list)){
    indices = which(is.finite(list_of_Js[[species]]) & list_of_Js[[species]] != 0, arr.ind=T)
    values = list_of_Js[[species]][indices]
    
    top_indices = order(abs(values), decreasing=T)[1:N]#Selects top J parameters to keep most important.
    indices = indices[top_indices,]
    values = values[top_indices]

    index_names = paste(c(rep('', 12), rep('+', 24)), c(ceiling(-(15:4)/4), floor((4:27)/4)), ':', rep(c('A', 'C', 'G', 'T'), 9), sep='')
    library(BioCircos)

    custom_genome = as.list(t(fi_list[[species]])[1:36])
    names(custom_genome) = index_names
    not_consenso = which(!(as.numeric(custom_genome) > 0.25))

    links_chromosomes_1 = index_names[indices[,1]]#Start
    links_chromosomes_2 = index_names[indices[,2]] #End
    
    n_bins = 10
    link_colors = rep('', length(values))
    link_background_color = 'white'
    ramp_palette_positive = colorRampPalette(c(link_background_color, 'green'))(n_bins+1)[2:(n_bins+1)]
    link_colors[which(values > 0)] = ramp_palette_positive[cut(values[which(values > 0)], n_bins, include_lowest=T, labels=1:n_bins)]
    ramp_palette_negative = colorRampPalette(c(link_background_color, 'red'))(n_bins+1)[2:(n_bins+1)]
    if(length(which(values < 0)) > 0){
      link_colors[which(values < 0)] = ramp_palette_negative[cut(abs(values[which(values < 0)]), n_bins, include_lowest=T, labels=1:n_bins)]
    }
    
    posiciones = cumsum(as.numeric(custom_genome))
    posiciones_resta = abs(c(0, posiciones[1:(length(posiciones)-1)]) - posiciones)
    links_pos_1 = posiciones_resta[indices[,1]]
    links_pos_2 = posiciones_resta[indices[,2]]
    links_labels = rep('', length(links_chromosomes_1))

    null = rep(0, length(links_pos_1))

    tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 1, borderSize = 0, fillColors = link_background_color)  

    different_colors = unique(link_colors)

    for(color in different_colors){
        which_have_color = which(link_colors == color)
        tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1[which_have_color], null[which_have_color], links_pos_1[which_have_color], links_chromosomes_2[which_have_color], null[which_have_color], links_pos_2[which_have_color], maxRadius = 0.95, labels = links_labels[which_have_color], color=color)
    }

    exon_site_palette = colorRampPalette(c('gold', 'darkorange2'))(3)
    intron_site_palette = colorRampPalette(c('cyan4', 'darkorchid4'))(6)
    genome_colors = c(exon_site_palette, intron_site_palette)[floor((0:35)/4) + 1]
    genome_colors[not_consenso] = paste(genome_colors[not_consenso], '60', sep='')

    plot = BioCircos(tracklist, genome = custom_genome, genomeFillColor = genome_colors, chrPad = 0.035, displayGenomeBorder = FALSE, yChr =  FALSE, genomeTicksDisplay = FALSE, genomeLabelDisplay = TRUE, genomeLabelOrientation = 270, genomeLabelTextSize = '16pt', genomeLabelDy = 18, width='1000px', height='1000px')
    
    library(htmlwidgets)
    saveWidget(plot, paste('circos_', sp_list[species], '.html', sep=''), selfcontained=T)
}

