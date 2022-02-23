#define N_SITES 9
N termalizacion_sin_cambio =1000;
#define MAX_NO_CONVERGENCE 100000000
#define epsilonSingle 0.2
#define epsilonJoint 0.2
#define tolNum 0.0000001
#define TOLERANCE_HI 0.005
#define TOLERANCE_JIJ 0.01
#define tolabsHi 0.000001
#define tolabsJij 0.000001
#define memory_drag 0.9
#define tolPi 0.001
#define tolPij 0.001
#define SavingDataSteps 10
#define SaveSequencesRate 10
#define TOTAL_SEQUENCES 130000
#define MUTATION_STEPS 10000
std::string filename_marginales = "../data/freqs/seqs_gr5/fi.txt";
std::string filename_dobles = "../data/freqs/seqs_gr5/fij.txt";
std::string filename_Jij = "";
std::string filename_hi = "";
std::string prefixDirectory = "runs/seqs_gr5/";
std::string filename_errorfile= "runs/seqs_gr5//error.txt";
#define gamma_start 0.05
#define gamma_end 0.01
#define gamma_step 0.005
bool through_zero = false;
