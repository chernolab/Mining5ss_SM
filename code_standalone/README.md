C++ code to get a family of regularized maximum entropy models from observed single-site and two-site frequencies. In a typical run, models are obtained sequentially increasing the regularization level. In this way, after the last fitting iteration at a given gamma value the algorithm  authomatically continue the fitting procedure with the next specified regularization level.


# Compilation

The **main.cpp** file should be compiled using

```
g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall
```

# Dependencies

The **include** folder contains include files where parameters and auxiliary functions are declared.

* **parameters.h**: many input/output and running parameters are declared (e.g. input-output files and directories, regularization gamma values, 
output saving rates, and convergence criteria for error values and parameter changes).  For instance, fhe frequency at which parameters and probabilities are saved to disk can be specified at line 26 of **parameters.h** (SavingDataSteps). The frequency at which sequences and energies are saved to disk can be specified at line 27 of **parameters.h** (SaveSequencesRate). See **parameters.h** for further details.


* **include/dirtools.h** (should not be modified): Functionality to check for directories used for output, and creates them should they not exist.

* **include/lectura_escritura.h** (should not be modified): Functionality to read input files and write output files.

* **include/generacion_secuencia** (should not be modified) Functionality to produce ensembles of sequences by a Markov Chain Montecarlo algorithm through rejection sampling.

# Input

* Empirical one-site and two-site frequencies files (for instance f<sub>i</sub> and f<sub>ij</sub> provided at data/freqs/hsa for Homo Sapiens) should be specified at 
lines 32 and 33 of the **parameters.h** file (std::string filename_marginales and std::string filename_dobles respectively). 

* Previously estimated single-site fields and two-site couplings could be optionally specified at lines 35 and 36 of **parameters.h** 
(std::string filename_Jij and std::string filename_hi respectively) to resume previous runs.

# Output

Outputs will be stored in the appropriate runs/ subfolder (e.g. runs/hsa for Homo Sapiens). Within each species' folder, there exists a separate folder for each
gamma value run (e.g. runs/hsa/gamma0.050000/)

* ParamH : single site fitting parameters
* ParamJ : two-site fitting parameters
* Pi     : ensemble estimated single site probabilities
* Pij    : ensemble estimated two-site probabilities
* sequeneces : ensemble sequences
* energies   : data-driven energy of ensemble sequences
* **error.txt**  : error information for the entire run. It contains, in the following order in each line, the current gamma value, the number of 
nonconverged H and J parameters, the number of nonzero J parameters, and the maximum error for Pis and Pijs respectively (6 entries, tab separated format).



