directory_fis = 'data/seqs/'
organisms = list.files(directory_fis)
organisms = organisms[grep(organisms, pattern='\\.rds')]
sp_list = gsub(organisms, pattern='logo_|\\.rds', replacement='')

seqs_list = lapply(paste(directory_fis, organisms, sep=''), FUN=readRDS)
names(seqs_list) = sp_list

num_GT = lapply(seqs_list, FUN=function(x){return(sum(x[,4] == 'G' & x[,5] == 'T'))})
num_total = lapply(seqs_list, FUN=nrow)
num_unique_GT = lapply(seqs_list, FUN=function(x){return(sum(unique(x)[,4] == 'G' & unique(x)[,5] == 'T'))})
num_unique = lapply(seqs_list, FUN=function(x){return(nrow(unique(x)))})

df_nums = data.frame(TOTAL=unlist(num_total), GT=unlist(num_GT), UNIQUE=unlist(num_unique), UNIQUE.GT=unlist(num_unique_GT), PROP=unlist(num_GT)/unlist(num_total))
rownames(df_nums) = names(seqs_list)

write.table(df_nums, file='seq_nums.csv', quote=F, row.names=T, col.names=T)