# Mining5ss_SM
Supplementary material for Mining 5ss paper

## Directory Guide

**code_standalone/**: C++ code and program used for fitting observed single-site and two-site frequencies into an H-parameter and regularized J-parameter model. This folder contains all necessary code to compile and run for a single fitting session.

**code_compileNrun**: Same as above, but with parameter generation script for ease of use.

**data/**: Directory containing data used for the fitting algorithm, fed into C++ programs present in **code_standalone/** and **code_compileNrun/**. **data/freqs/** contains single-site and two-site observed frequencies given to the fitting program, under a subfolder pertinent to the species needed. **data/seqs/** contains the sample sequences used to construct those observed values, in **logo_[species]** files named after each species.

**runs/**: is the standard directory where fitting runs will be stored, although other locations may be specified within the C++ code. It contains subfolders for each species, and within each subfolder, **runs/[species]/error.txt** tracks the evolution of error as the run went on. Within each species' subfolder are many folders where the results for the run at different gammas is stored. Because of the decreasing order with which gammas are applied (as regularization increases), higher gammas (typically **gamma0.05000**) occur before lower gammas (typically **gamma0.010000**).

**scripts/**: contains the R scripts used to create the circos diagrams for exploring divergent signals in different organisms, and the table used to analyse convergent signals. These scripts are intended to be used using the base directory as their working directory. This is to say, it's recommended to use **Mining5ss_SM** as the working directory for these scripts.

##C++ code

The C++ code manages the task of inputting observed statistical data from sequences, and using it to fit a regularized parameter model.

Within the **code_standalone** directory, we can find the **main.cpp** file. This should be compiled using

```
g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall
```

It is recommended to use locally, but if one wishes to use it elsewhere one may, making sure to change the appropriate directories used within the code.

In order to do this, one needs to explore the **include** folder, which contains files necessary to run the code.

Within **include/parameters.h**, one can find a variety of parameters used, including the directories that one can choose to change.

The user is recommended not to edit the other files within the **include/** folder. These are:

**include/dirtools.h**: Used to check for directories used for output, and creates them should they not exist.

**include/lectura_escritura.h**: Used to read input files and write output files.

**include/generacion_secuencia**: Used to evolve sequences by a Markov Chain Montecarlo algorithm through rejection sampling, and later generate a sequence ensemble. Both of these are performed using current fitted parameters, and the results are used to compare the simulated statistics of the system with the observed statistics of real sequences. This is used for the fitting algorithm, so that fitted parameters slowly approximate those best suited to represent real-world data.

Most of the fitting work not covered in **include/generacion_secuencia** is covered within the **main.cpp** file, in the **fit_iteration()** function. Both **main.cpp** and **include/generacion_secuencia** are commented for ease of understanding.

##R Scripts

Within the directory **scripts/** we can find two R scripts.

**scripts/convergent_table.R** is used to produce the table comparing Intron-Exon signals for different species. This table is used within the article to illustrate similar behaviors among different species.

**scripts/circos.R** is used to produce the circos diagrams which illustrate the important J parameters fitted for different organisms, as well as the proportion of times a given nucleic acid was represented at a given site.

Both scripts should be run using the base directory as their working directory, but parameters within the files can be changed in order to use them elsewhere.