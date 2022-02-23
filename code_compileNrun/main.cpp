  // Program implementing fitting algorithm.
  // Hamiltonian (energy) operator:  sum h_i + sum sum J_{ij}, fitted as described in paper.
  // Compiled with command g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall
  // Flags indicate the following: Use appropriate c++ version, compile any possible algorithms in parallel (see https:  //gcc.gnu.org/onlinedocs/libstdc++/manual/parallel_mode.html), optimize for speed, show all compiler warnings.

  // Starting values for h_i and J_{ij} can be specified in files given in parametros.h.

  // Convergence criterion is based on rate of change of h_i and J_{ij}.

#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <ctime>
#include <random>
#include <functional>   //For std::bind
#include <chrono>
#include <map>
#include <algorithm>   //For execution compatible with parallelization.
#include <numeric>

#include <float.h> //For comparison to -inf.

#include "./include/dirtools.h"  //For preparing the directories needed.

using namespace std;

typedef int N;  //Can be substituted for longs if needed.
typedef double T;  //Can be substituded for floats or other floating point types if needed.

#include "./parameters.h"     //File containing parameters that can be changed between different fitting runs.

std::uniform_int_distribution<N> uni20(0,3);   // Generate random integer values for nucleic acid bases: (0, 1, 2, 3) -> (A, C, G, T).
std::uniform_int_distribution<N> uni65(0,N_SITES-1);   // Generate random integer values for nucleic acid locations: (0-8) -> (1-9). N_SITES = 9 located in parametros.h
std::uniform_real_distribution<T> dist(0.0, 1.0);   // Generate random floating-point values between 0 and 1.

const N naa = 4;  //Number of nucleic acid bases.
const N naanpos = naa * N_SITES;  //Columns and rows of the matrix Pij/Jij.

  //Model parameter matrices (Initially hi = log (fi) and Jij=0).
  //We use the vector class and a nested vector so as to compatibilize with operators within the algorithm library.
std::vector<T> hi(naanpos);
std::vector<std::vector<T>> jij(naanpos, std::vector<T>(naanpos));

N numHiNonConverged;
N numJijNonConverged;

const int ncomb = N_SITES * (N_SITES + 1)/2 - N_SITES;  //Number of possible combinations of positions (1-1, 1-2, ... , 1-naa, 2-2,2-3,... etc).
int PosPos[ncomb][2];  //All possible combinations of positions.

#include "./include/lectura_escritura.h"  //File containing functions used for reading and writing files.
#include "./include/generacion_secuencia.h"  //File containing the code used to simulate the evolution of genetic sequences and generate a number of samples from their distribution.

  //One fitting iteration. Returns true on convergence and false otherwise.
template<typename N, typename T>
bool fit_iteration(N counter, T gamma, std::vector<T> p1, std::vector<std::vector<T>> pij, std::vector<T> f1Data, std::vector<std::vector<T>> f2Data, T &maxPi, T &maxPij, T &memory_change_pi, T &memory_change_pij, std::chrono::high_resolution_clock::time_point t1, N &numHiNonConverged, N &numJijNonConverged, N numJijNonzero, FILE * errorFile){
    T change;

    bool writeFiles = (counter % SavingDataSteps == 1);

    std::ofstream sequencesFile;
    std::ofstream energyFile;

    if(counter % SaveSequencesRate == 1){  //Sequence ensemble and energy distribucion is saved in files this often.
        cout << counter << endl;
          //File to save sequences.
        std::string fullnameSequences = prefixDirectory;
		fullnameSequences.append("gamma");
		fullnameSequences.append(std::to_string(gamma));
		fullnameSequences.append("/sequences");
        fullnameSequences.append(std::to_string(counter));
		std::cout << fullnameSequences << std::endl;
        sequencesFile.open(fullnameSequences);

          //File to save sequence energies.
        std::string fullnameEnergy = prefixDirectory;
        fullnameEnergy.append("gamma");
        fullnameEnergy.append(std::to_string(gamma));
        fullnameEnergy.append("/energy");
        fullnameEnergy.append(std::to_string(counter));
        energyFile.open(fullnameEnergy);
    }

      //Reset joint probability (Pij) values to zero in order to later calculate them.
    std::fill(pij.begin(), pij.end(), std::vector<T>(naanpos, 0.0));

    int sign_change, sign_j;
    std::vector<T> energies(TOTAL_SEQUENCES);
    std::vector<std::vector<N>> sequences(TOTAL_SEQUENCES);

    generar_secuencias<N,T>(TOTAL_SEQUENCES, MUTATION_STEPS, sequences, energies);  //Generate sequence ensemble, saved within the variable "sequences". Their energies are saved in "energies".

    //Recalculate Pij values by summing occurrences.
    for(unsigned int i=0; i<sequences.size(); i++){
        vector<N> sequence = sequences[i];
        T energy = energies[i];

          // Actualizo Pij con los valores de esta secuencia
        for(int j=0; j<ncomb; j++){
            N index1 = PosPos[j][0] * naa + sequence[PosPos[j][0]];
            N index2 = PosPos[j][1] * naa + sequence[PosPos[j][1]];
            pij[index1][index2]++;
            if(index1!=index2){
                pij[index2][index1]++;
            }
        }

        if(sequencesFile.is_open()){  //Write sequence and its energy into files.
            for (int j=0; j<N_SITES; j++){
                sequencesFile << sequence[j] << " ";
            }
            sequencesFile << endl;
            energyFile << energy << endl;
        }
    }

    if(sequencesFile.is_open()){  // Close files that will no longer be in use.
        sequencesFile.close();
        energyFile.close();
    }

    // Normalization of Pij.
    for(int index1=0;index1<naanpos; index1++){
        for(int index2=0;index2<naanpos;index2++){
            pij[index1][index2] = pij[index1][index2]/TOTAL_SEQUENCES;
        }
    }

    // Calculate Pi values as partial sums of Pij.
    for(int indA=0;indA<naanpos;indA++){
        p1[indA]=0;
    }

    // For position 1.
    for(int indA=0;indA<naa;indA++){
        for(int indx=naa;indx<(2*naa);indx++){
            p1[indA] += pij[indA][indx];
        }
    }

    // For other positions.
    for(int indA=naa;indA<naanpos;indA++){
        for(int indx=0;indx<naa;indx++){
            p1[indA] += pij[indx][indA];
        }
    }

      // Compare Pij and Pi (from the model) to fij and fi (from data).
    double m = 0.0;
    double mp = 0.0;
    int n1, n2;
    for(int fila=0;fila<naanpos;fila++){
        n1 = fila/naa;
        for(int col=0;col<naanpos;col++){
            n2 = col/naa;
            if(n1 != n2){
                m = fabs(f2Data[fila][col] - pij[fila][col]);
                if(m > mp){
                    mp = m;
                }
            }
        }
    }

	  //This indicator takes into account the maximum error of Pij, but averages it with its previous values so as to have a smooth estimate of error.
    memory_change_pij = memory_drag * (memory_change_pij) + (1.0 - memory_drag) * fabs(mp - maxPij);
    maxPij = mp;

    m = 0.0;
    mp = 0.0;
    for(int i=0;i<naanpos;i++){
        m = fabs(f1Data[i]-p1[i]);
        if(m>mp){
            mp = m;
        }
    }

	  //Same indicator for Pi.
    memory_change_pi = memory_drag * (memory_change_pi) + (1.0 - memory_drag) * fabs(mp - maxPi);
    maxPi = mp;

    bool error = true;
	  //The convergence criterion requires all Hi and Jij parameters to have converged, but also for the last few error values for Pij and Pi to be within acceptable tolerance parameters.
    if(numHiNonConverged == 0 && numJijNonConverged == 0 && memory_change_pij < tolPij && memory_change_pi < tolPi){
        error = false;
    }

    if(writeFiles){  //Write output of parameters and probabilities.
        std::string fullnameHI = prefixDirectory;
        std::string fullnameJIJ = prefixDirectory;
        std::string fullnamePI = prefixDirectory;
        std::string fullnamePIJ = prefixDirectory;

        fullnameHI.append("gamma");
        fullnameJIJ.append("gamma");
        fullnamePI.append("gamma");
        fullnamePIJ.append("gamma");

        fullnameHI.append(std::to_string(gamma));
        fullnameJIJ.append(std::to_string(gamma));
        fullnamePI.append(std::to_string(gamma));
        fullnamePIJ.append(std::to_string(gamma));

	fullnameHI.append("/ParamH");
	fullnameJIJ.append("/ParamJ");
	fullnamePI.append("/Pi");
	fullnamePIJ.append("/Pij");

        fullnameHI.append(std::to_string(counter));
        fullnameJIJ.append(std::to_string(counter));
        fullnamePI.append(std::to_string(counter));
        fullnamePIJ.append(std::to_string(counter));

        write_1(fullnameHI, hi);
        write_2(fullnameJIJ, jij);
        write_1(fullnamePI, p1);
        write_2(fullnamePIJ, pij);
    }

    //Redefine parameters values (hi and Jij), if error is still large.
    //The Jij are defined according to L1 regularization.
    if(error){
	numJijNonzero = N_SITES * (N_SITES - 1) * 4 * 4;
        numHiNonConverged  = 0;
        numJijNonConverged = 0;

        for(int i=0;i<naanpos;i++){
            change = epsilonSingle * (f1Data[i] - p1[i]);
            if(abs(change) > abs(hi[i] * TOLERANCE_HI) + tolabsHi){
                numHiNonConverged++;
	    }
            hi[i]+=change;
        }

        double auxJIJ;
	std::vector<std::vector<T>> prevJij = jij;
        int n1, n2;

        for(int fila=0; fila<naanpos; fila++){
            n1 = fila / naa;
            for(int col=0; col<naanpos; col++){
                n2 = col / naa;
                if(n1 != n2){
		    if(f2Data[fila][col] < pij[fila][col]){  //This defines the sign of the difference between observed relative frequency and fitted probability.
                        sign_change = -1;
                    }
                    else{
                        sign_change = 1;
                    }

                    if(jij[fila][col] == 0){  //If the J parameter is zero.
                        if(abs(f2Data[fila][col] - pij[fila][col]) < gamma){  //In this case, the zero parameter will not become nonzero because the motive to do so was not strong enough.
                            jij[fila][col] = 0;
			    numJijNonzero--;
                        }
                        else{//Otherwise, it is changed and becomes nonzero.
			    change = epsilonJoint * ((f2Data[fila][col] - pij[fila][col]) - gamma * sign_change);
                            if(abs(change) > fabs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){  //If it changes more than a given amount, it is not considered converged.
                                numJijNonConverged++;
			    }
                            jij[fila][col] =+ change;
                        }
                    }
                    else{//If the J parameter is already nonzero.
                    	if(jij[fila][col] < 0){  //Calculate the current sign of the J parameter.
                            sign_j = -1;
                        }
                    	else{
                            sign_j = 1;
                   	}

			change = epsilonJoint * (f2Data[fila][col] - pij[fila][col] - gamma * sign_j);
                        auxJIJ = jij[fila][col] + change;
                        if(auxJIJ * jij[fila][col] < 0){  //If it would cross the zero line.
                            if(abs(jij[fila][col]) > abs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){
                                numJijNonConverged++;
			    }
			    jij[fila][col] = 0;  //It will become zero.
			    numJijNonzero--;
                        }
                        else{  //If it would not cross the zero line.
			    change = epsilonJoint * (f2Data[fila][col] - pij[fila][col]);  //It evolves in the usual manner.
			    if(abs(change) > abs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){
                                numJijNonConverged++;
			    }
                            jij[fila][col] += change;
                        }
                    }
                }
                else{
                    jij[fila][col] = 0;  //Elements within the same site should be null (for example, 1:A and 1:G should have no interaction).
                }
            }
        }

          //Time benchmarking.
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        cout << "... Tiempo: " << duration/1000000.0L << endl;
    }

    fprintf(errorFile, "%f\t%d\t%d\t%d\t%f\t%f\n", gamma, numHiNonConverged, numJijNonConverged, numJijNonzero, maxPi, maxPij);
	printf("Current gamma value: %f\t H parameters not converged: %d\t J parameters not converged: %d\n", gamma, numHiNonConverged, numJijNonConverged);
	printf("Nonzero J parameters: %d\t Maximum Pi error: %f\t Maximum Pij error: %f\n", numJijNonzero, maxPi, maxPij);  //Live printing of values.
    return error;
}

int main(int argc, char* argv[]){
    double maxPij = 0.00;
    double maxPi = 0.00;

    /* ----------------------------------------------------------------------------- */
    /* 1- LOAD DATA AND DEFINE PARAMETERS .......................................... */
    /* ----------------------------------------------------------------------------- */
      //Parametros generales  | General parameters.
    static std::vector<T> p1(naanpos);
    static std::vector<std::vector<T>> pij(naanpos, std::vector<T>(naanpos));

    int a = 0;  //Auxiliary counting variable.
      //Save every possible two-site combination.
    for(int i=0;i<N_SITES;i++){
        for(int j=i;j<N_SITES;j++){
            if(i != j){
                PosPos[a][0]=i;
                PosPos[a][1]=j;
                a++;
            }
        }
    }

      //Read experimental frequencies data.
    std::vector<T> f1Data(naanpos);
    std::vector<std::vector<T>> f2Data(naanpos, std::vector<T>(naanpos));
    load_f1<T>(filename_marginales, f1Data);
    load_f2<T>(filename_dobles, f2Data);

    /* ----------------------------------------------------------------------------- */
    /* 2- BEGINS BIG LOOP FOR PARAMETERS ACTUALIZATION ............................. */
    /* ----------------------------------------------------------------------------- */
    /* MC parameters */

    ifstream j0File(filename_Jij);
    if(j0File.good()){
        cout << "Initial Jij values provided." << endl;
        load_f2(filename_Jij, jij);
    }
    else{
	for(int i=0; i<naanpos; i++){
	    int n1 = i / naa;
	    for(int j=0; j<naanpos; j++){
		int n2 = i / naa;
		if(n1 != n2){
		    jij[i][j] = 0;
		}
		else{
		    jij[i][j] = 0;
		}
	    }
	}
    }
    j0File.close();

    ifstream h0File(filename_hi);
    if(h0File.good()){
        cout << "Initial hi values provided." << endl;
        load_f1(filename_hi, hi);
    }
    else{
	for(int i=0; i<naanpos; i++){
		hi[i] = log(f1Data[i]);
	}
    }
    h0File.close();

      //Makes sure the directory for saving data files exists.
    prepare_directory(prefixDirectory.c_str());
    printf("%s\n", prefixDirectory.c_str());
    float gmax = gamma_end;
    float gmin = gamma_start;
    if(gamma_end < gamma_start){
	gmax = gamma_start;
	gmin = gamma_end;
    }
    for(float g=gmin; g<=gmax; g+=gamma_step){
	char new_gdir[200];
	sprintf(new_gdir, "%sgamma%f", prefixDirectory.c_str(), g);
        printf("%s\n", new_gdir);
        prepare_directory(new_gdir);
    }

    bool error = true;
    int counter = 0;
    int min_iterations = 10;

    std::ofstream sequencesFile;
    std::ofstream energyFile;

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    int n_gamma_steps;
    if(through_zero){
        n_gamma_steps = (int)((gamma_start + gamma_end)/gamma_step);
    }
    else{
        n_gamma_steps = (int)(abs(gamma_start - gamma_end)/gamma_step);
    }
    FILE * errorFile = fopen(filename_errorfile.c_str(), "w");
    for(int g=0; g<=n_gamma_steps; g++){
        double local_gamma = 1.0;  //The first iteration is performed with a huge value of gamma in order to obtain a high-regularization scenario before starting to lower it.
        if(g >= 0){
            if(through_zero){  //Deprecated scenario use where gamma was lowered to zero and then increased (an annealing-quenching scenario).
                local_gamma = abs(gamma_start - ((double) g) * gamma_step);
            }
            else{  //Otherwise used for a monotonous increase in gamma.
                if(gamma_start > gamma_end){  //The use given in the paper, where gamma is decreased.
                    local_gamma = gamma_start - ((double) g) * gamma_step;
                }
                else{  //Deprecated scenario where gamma is increased.
                    local_gamma = gamma_start + ((double) g) * gamma_step;
                }
            }
        }
        printf("Gamma: %f\n", local_gamma);
        N numHiNonConverged = naanpos;
        N numJijNonConverged = naanpos*naanpos;

        double memory_change_pi = 2 * tolPi;
    	double memory_change_pij = 2 * tolPij;

        maxPi = memory_change_pi;
    	maxPij = memory_change_pij;
        counter = 0;
		N numJijNonzero = 0;
        while(error || counter < min_iterations){  //We demand a minimum number of iterations so that parameters are allowed to change.
            counter++;
            error = fit_iteration<N,T>(counter, local_gamma, p1, pij, f1Data, f2Data, maxPi, maxPij, memory_change_pi, memory_change_pij, t1, numHiNonConverged, numJijNonConverged, numJijNonzero, errorFile);
			if(counter > MAX_NO_CONVERGENCE && g != 0){  //If it hasn't converged in a long time, gamma is too low (regularization too lax).
				if(through_zero){  //Deprecated use for the annealing-quenching case.
					g = (int)((local_gamma + gamma_start)/gamma_step) + 1;
				}
				else{  //Actual use given where gamma only decreases.
					g = n_gamma_steps;
				}
			}
        }
    }
    fclose(errorFile);

    return 0;
}
