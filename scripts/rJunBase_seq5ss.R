#Obtener junturas expresadas en humanos para analizar sus 5'ss

#Libraries
library(limma)
library("BSgenome.Hsapiens.NCBI.GRCh38")

?#Descargo de rJunBase la tabla con todas las junturas lineales
#http://www.rjunbase.org/rJunBase/download/toIndex
rjunBase = read.delim("detail_LS_annotation.txt")
#nrow = 682017

#Algunos filtros para quedarme con junturas mas confiables: me quedo con junturas anotadas y que pertenecen a genes que codifican proteinas

rjunBase <- rjunBase[rjunBase$Annotation == "yes", ]
#nrow = 357366

rjunBase <- rjunBase[rjunBase$Gene.type == "protein coding gene", ]
#nrow = 309029
#N de genes = 18670

#Separo el dato de la localizacion de las junturas en distintas columnas
tmp <- strsplit2(rjunBase$Junction.location, split = ":")
tmp <- data.frame(tmp[, 1], strsplit2(tmp[, 2], split = "\\|"), tmp[, 3])

colnames(tmp) <- c("chr", "start", "end", "strand")

#Modifico el nombre de los chr
tmp$chr <- gsub(pattern = "chr", replacement = "", x = tmp$chr)

#Agrego la info a la tabla principal
rjunBaseCoor <- cbind(rjunBase, tmp)
rm(tmp)

#Elimino las junturas en las que coincidan el sitio 5'ss y obtengo las secuencias
pos     <- rjunBaseCoor[rjunBaseCoor$strand == "+", ]
pos     <- pos[!duplicated(paste(pos$chr, pos$start)), ]
seq_pos <- getSeq(Hsapiens, 
                  names = pos$chr, 
                  start = (as.numeric(pos$start) - 2), 
                  width = rep(9, length(as.numeric(pos$start))))

neg     <- rjunBaseCoor[rjunBaseCoor$strand == "-", ]
neg     <- neg[!duplicated(paste(neg$chr, neg$end)), ]
seq_neg <- reverseComplement(getSeq(Hsapiens, 
                                    names = neg$chr, 
                                    start = (as.numeric(neg$end) - 6), 
                                    width = rep(9, length(as.numeric(neg$start)))))

#Mergeo de nuevo las junturas y agrego las secuencias a la tabla
rjunBaseCoor        <- rbind(pos, neg)
rjunBaseCoor$seq5ss <- c(as.character(seq_pos), as.character(seq_neg))
#nrow=239611

save(rjunBaseCoor, file = "rJunBase_seq5ss_noanotate.RData")


