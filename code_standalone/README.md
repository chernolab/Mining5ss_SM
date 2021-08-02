# Compilation

The **main.cpp** file should be compiled using

```
g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall
```

# Dependencies

The **include** folder contains include files where parameters and auxiliary functions are declared.

* **parameters.h**: many input/output and running parameters are declared (e.g. input-output files and directories, regularization gamma values, 
output saving rates, and convergence criteria for error values and parameter changes). 

* **include/dirtools.h** (should not be modified): Functionality to check for directories used for output, and creates them should they not exist.

* **include/lectura_escritura.h** (should not be modified): Functionality to read input files and write output files.

* **include/generacion_secuencia** (should not be modified) Functionality to produce ensembles of sequences by a Markov Chain Montecarlo algorithm through rejection sampling.

# Input

* Empirical one-site and two-site frequencies (for instance f<sub>i</sub> and f<sub>ij</sub> provided at data/freqs/hsa for Homo Sapiens) should be specified at 
lines 32 and 33 of the **parameters.h** file (std::string filename_marginales and std::string filename_dobles respectively). 

* Previously estimated single-site fields and two-site couplings could be optionally specified at lines 35 and 36 of **parameters.h** 
(std::string filename_Jij and std::string filename_hi respectively) to resume previous runs.

# Output

Outputs will be stored in the appropriate runs/ folder (e.g. runs/hsa for Homo Sapiens). Within each species' folder, there exists a separate folder for each
gamma value used in the fitting process, and an error file. Within each gamma directory (example, runs/hsa/gamma0.050000/) are stored the different 
values for fitted H and J parameters (for example, /ParamH101, /ParamJ101 for the 101st iteration), the estimated single-site and two-site probabilities (/Pi101, /Pij101), 
and generated sequence and energy ensembles (/sequences101, /energies101) for each different fitting iteration within that gamma value.

There is a different storage rate for parameters and probabilities, compared to sequences and energies (for example, at current values the first are saved every 100 iterations, and the latter are not saved). 
This can be changed in parameters.h.

After the last fitting iteration of a given gamma value (suppose 0.05) the algorithm proceeds to use the same parameters to start fitting at a lower gamma value (0.045 in this case). 
Thus, we transition from the folder gamma0.050000/ to the folder gamma0.045000/.

The error file error.txt exists outside the gamma folders and stores error information for the entire run. It contains, in the following order in each line, the current gamma value, the number of 
nonconverged H and J parameters, the number of nonzero J parameters, and the maximum error for Pis and Pijs respectively (6 entries, tab separated format).

