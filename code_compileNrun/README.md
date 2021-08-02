We provide R-scripts and C++ code to infer multiple fitting models for different organisms at once.
Once **main.cpp** was compiled, you should source the **compileNrun.R** file to fit regularized maximum entropy models
for all the organisms specified in the first line of the R-script. For instance, the following line will consider all the organisms studied in the paper:

```
orgs    <- c('ath','mtr','osa',
             'hsa','pan','ggo','mmu','dre',
             'dme','dya','dps',
             'cel','cbr',
             'cne')
```



# Input
For each organism, ORG, fitting input parameters (e.g. input-output files and directories, regularization gamma values,
output saving rates, and convergence criteria for error values and parameter changes) should be specified in a dedicated **parameter_ORG.txt** file in the *include* 
folder. For instance **parameter_ath.txt** and **parameter_mtr.txt** correspond to input parameters for Arabidopsis thaliana and Medicago truncatula respectively. 
All the parametes_ORG.txt files used in the paper can be found in the *include* folder.

These files will be sourced from the **compileNrun.R** script and will autmoatically generate the appropiate include files for code compilation. 
The script will then launch the code compilation and the model estimation procedure for each considered organism.

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


