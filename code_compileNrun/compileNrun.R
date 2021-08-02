orgs    <- c('ath','mtr','osa',
             'hsa','pan','ggo','mmu','dre',
             'dme','dya','dps',
             'cel','cbr',
             'cne')
             
orgname <- orgs

for(iorg in seq_along(orgs)){
  # Header generation ----
  source(paste0('include/params_',orgs[iorg],'.txt'))
  {
    
    headerFileName <- 'parameters.h'
    
    cat(paste0('#define N_SITES 9\n'),file = headerFileName)  #Cantidad de sitios en las secuencias.
    cat(paste0('N termalizacion_sin_cambio =',nTermalizacion,';\n'),append=TRUE,file = headerFileName)
    
    cat(paste0('#define MAX_NO_CONVERGENCE ', format(maxNoConvergence,scientific=FALSE),'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define epsilonSingle ',epsilonSingle,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define epsilonJoint ',epsilonJoint,'\n'),append=TRUE,file=headerFileName)  
    
    cat(paste0('#define tolNum ',format(tolNum,scientific = FALSE),'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define TOLERANCE_HI ',TOLERANCE_HI,'\n'),append=TRUE,file=headerFileName) 
    cat(paste0('#define TOLERANCE_JIJ ', TOLERANCE_JIJ,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define tolabsHi ', format(tolabsHi,scientific = FALSE),'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define tolabsJij ', format(tolabsJij,scientific = FALSE),'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define memory_drag ',memory_drag,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define tolPi ', tolPi,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define tolPij ',tolPij,'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define SavingDataSteps ',SavingDataSteps,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define SaveSequencesRate ',SaveSequencesRate,'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define TOTAL_SEQUENCES ',TOTAL_SEQUENCES,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define MUTATION_STEPS ',MUTATION_STEPS,'\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('std::string filename_marginales = ', filename_marginales,';\n'),append=TRUE,file=headerFileName)
    cat(paste0('std::string filename_dobles = ', filename_dobles,';\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('std::string filename_Jij = ',filename_Jij,';\n'),append=TRUE,file=headerFileName)
    cat(paste0('std::string filename_hi = ',filename_hi,';\n'),append=TRUE,file=headerFileName)
    cat(paste0('std::string prefixDirectory = ',prefixDirectory,';\n'),append=TRUE,file=headerFileName)
    cat(paste0('std::string filename_errorfile= ',filename_errorfile,';\n'),append=TRUE,file=headerFileName)
    
    cat(paste0('#define gamma_start ',gamma_start,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define gamma_end ',gamma_end,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('#define gamma_step ',gamma_step,'\n'),append=TRUE,file=headerFileName)
    cat(paste0('bool through_zero = false;\n'),append=TRUE,file=headerFileName) #//Si bajo hasta el cero y vuelvo, o no.
    cat(orgs[iorg],' :header file generated.\n')
  }
  
  # Code compilation ----
  {
    saux <- paste0('g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall')
    system(saux)
    cat(orgs[iorg],': code compiled.\n')
  }
  
  # Runs ----
  {
    # . create output directories ----
    ppath   <- paste0('./corridas/',orgname[iorg],'_ach_h0/')
    ppath   <- sub('\"','',sub('\"','',prefixDirectory))
    gg       <- seq(gamma_start,gamma_end,gamma_step*sign(gamma_end-gamma_start))
    dirName <- paste0(ppath,'gamma',sprintf('%.6f',gg))
    saux    <- paste0("mkdir -p ",dirName)
    for(i in seq_along(saux)){
      system(saux[i])  
    }
    
    system(paste('cp',paste0('include/params_',orgs[iorg],'.txt'),ppath))
    
    # . start runs ----
    #saux <- paste0('./main_',orgs[iorg])
    saux <- './main'
    system(saux)
    
    # . gzip seqs and energy files ----
    cat(orgs[iorg],' :Compressing seq and energy files...\n')
    dirName <- paste0(ppath,'gamma',sprintf('%.6f',gg))
    sauxE   <- paste0("gzip ",dirName, "/ener*")
    sauxS   <- paste0("gzip ",dirName, "/seq*")
    for(i in seq_along(dirName)){
      system(sauxE[i])  
      system(sauxS[i])  
    }
  }
}
