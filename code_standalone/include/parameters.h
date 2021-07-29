//Model parameters
#define N_SITES 9                 //Number of sites
N termalizacion_sin_cambio = 1000;//Number of thermalization steps without changes after finding the lowest energy

//Fitting parameters
#define MAX_NO_CONVERGENCE 100000000

#define epsilonSingle 0.2 // epsilon to renew parameters hi
#define epsilonJoint 0.2  // epsilon to renew parameters jij

#define tolNum 0.0000001  //Numeric tolerance

// Convergence criterium: change < TOL * value + tolabs.
#define TOLERANCE_HI 0.005 // Fractional tolerane for parameter changes
#define TOLERANCE_JIJ 0.01
#define tolabsHi 0.000001  // Absolute tolerance for parameter changes
#define tolabsJij 0.000001

// Probability convergence: error should not decrease
// We consider a wheighted average between current and previous error values
// memoria_error = memoria_error_previo * memory_drag + (1.0 - memory_drag) * error_actual
#define memory_drag 0.9  // wheight of the previous error value
#define tolPi 0.001      // tolerance to the memory of error change to claim convergence
#define tolPij 0.001

#define SavingDataSteps 10     // After how many steps data should be saved
#define SaveSequencesRate 10   // After how many steps sequences and energiers should be saved
#define TOTAL_SEQUENCES 130000 // Number of sequence in the ensemble
#define MUTATION_STEPS 10000

//Files and file-paths
std::string filename_marginales = "../data/freqs/hsa/fi.txt";   // empirical one-site frequencies
std::string filename_dobles = "../data/freqs/hsa/fij.txt";      // empirical two-site frequencies

std::string filename_Jij = "";                                // Previous Jij estimations to start with
std::string filename_hi =  "";                                // Previous hi estimations to start with

std::string prefixDirectory   = "../runs/hsa_new_test/";          // Output path
std::string filename_errorfile= "../runs/hsa_new_test/error.txt"; // Output error file

// Regularization parameters
#define gamma_start 0.05    // Starting gamma value
#define gamma_end 0.01      // Final gamma value
#define gamma_step 0.005    // gamma steps (absolute value)
bool through_zero = false;  // Cross zero gamma value
