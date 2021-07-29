//check
// Aprende parametros de Hamiltoniano
// Hamiltoniano:  sum h + sum sum J : Agrego multiplicadores de Lagrange de pId para la evolucion de parametros
// Compilar con g++ -std=c++11 -D_GLIBCXX_PARALLEL -fopenmp -o main main.cpp -Ofast -Wall
// Significado: Usar todo lo disponible hasta el c++ mas avanzado, compilar paralelizando lo posible (ver https://gcc.gnu.org/onlinedocs/libstdc++/manual/parallel_mode.html), optimizar por rapidez, mostrar todas las advertencias del compilador.

// Valores iniciales para Hi y Jij pueden especificarse en files:
//   Final_h_v1.txt y Final_jij_v1.txt

// Verificar para cada corrida las lineas con *check*

// 20180201: Stop criterium based on rate of change for h's and j's
// 201711  : Gamma value is required from command-line

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
#include <functional> //Para std::bind
#include <chrono>
#include <map>
#include <algorithm>
#include <numeric>

using namespace std;

typedef int N;
typedef double T;

#include "./parameters.h"

std::uniform_int_distribution<N> uni20(0,3); // genera uniformemente numeros enteros entre esos dos valores (cant de bases -1)
std::uniform_int_distribution<N> uni65(0,N_SITES-1); // genera uniformemente numeros enteros entre esos dos valores (cant de posiciones-1)
std::uniform_real_distribution<T> dist(0.0, 1.0); //genera uniformemente numeros reales entre 0 y 1

const N naa = 4;//Numero de aminoacidos. | Number of aminoacids.
const N naanpos = naa * N_SITES;//Columnas y filas de la matriz Pij | Columns and rows of the matrix Pij.

//Matrices de parametros del modelo. | Model parameter matrices (Initially hi = log (fi) and Jij=0).
std::vector<T> hi(naanpos);
std::vector<std::vector<T>> jij(naanpos, std::vector<T>(naanpos));

N numHiNonConverged;
N numJijNonConverged;

const int ncomb = N_SITES * (N_SITES + 1)/2 - N_SITES;//number of possible combinations of positions (1-1, 1-2, ... , 1-naa, 2-2,2-3,... etc)
int PosPos[ncomb][2];//all possible combinations of positions

#include "./include/lectura_escritura.h"
#include "./include/generacion_secuencia.h"

//Una iteracion de ajuste. Retorna un booleano que informa si ya convergio o no.
template<typename N, typename T>
bool fit_iteration(N counter, T gamma, std::vector<T> p1, std::vector<std::vector<T>> pij, std::vector<T> f1Data, std::vector<std::vector<T>> f2Data, T &maxPi, T &maxPij, T &memory_change_pi, T &memory_change_pij, std::chrono::high_resolution_clock::time_point t1, N &numHiNonConverged, N &numJijNonConverged, N numJijNonzero, FILE * errorFile){
    T change;

    bool writeFiles = (counter % SavingDataSteps == 1);

    std::ofstream sequencesFile;
    std::ofstream energyFile;

    if(counter % SaveSequencesRate == 1){//Tasa de guardado de secuencias y energias.
        cout << counter << endl;
        //File to save sequences:
        std::string fullnameSequences = prefixDirectory;
	fullnameSequences.append("gamma");
	fullnameSequences.append(std::to_string(gamma));
	fullnameSequences.append("/sequences");
        fullnameSequences.append(std::to_string(counter));
	std::cout << fullnameSequences << std::endl;
        sequencesFile.open(fullnameSequences);

        //File to save sequences' energies
        std::string fullnameEnergy = prefixDirectory;
        fullnameEnergy.append("gamma");
        fullnameEnergy.append(std::to_string(gamma));
        fullnameEnergy.append("/energy");
        fullnameEnergy.append(std::to_string(counter));
        energyFile.open(fullnameEnergy);
    }

    /* a) Recalculate Pij and Pi with the new parameters */
    // To do this I generate a large set of sequences with MC

    //Reinicio valores de pij a cero.
    std::fill(pij.begin(), pij.end(), std::vector<T>(naanpos, 0.0));

    int sign_change, sign_j;
    std::vector<T> energies(TOTAL_SEQUENCES);
    std::vector<std::vector<N>> sequences(TOTAL_SEQUENCES);

    generar_secuencias<N,T>(TOTAL_SEQUENCES, MUTATION_STEPS, sequences, energies);

    //Once every SaveSequenceRate, I refresh Pij values |
    //Renuevo los valores de Pij con los valores de una secuencia cada SaveSequencesRate
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

        if(sequencesFile.is_open()){//Write sequence and its energy into files.
            for (int j=0; j<N_SITES; j++){
                sequencesFile << sequence[j] << " ";
            }
            sequencesFile << endl;
            energyFile << energy << endl;
        }
    }

    if(sequencesFile.is_open()){
        sequencesFile.close();
        energyFile.close();
    }

    //Normalization of Pij | Normalizo los Pij.
    for(int index1=0;index1<naanpos; index1++){
        for(int index2=0;index2<naanpos;index2++){
            pij[index1][index2] = pij[index1][index2]/TOTAL_SEQUENCES;
        }
    }

    //Calculate Pi values as partial sums of Pij | Calculo los nuevos Pi como suma parcial de Pij.
    //Reinicio valores de p1 a cero.
    for(int indA=0;indA<naanpos;indA++){
        p1[indA]=0;
    }
    // For position 1 | Los de la posicion 1
    for(int indA=0;indA<naa;indA++){
        for(int indx=naa;indx<(2*naa);indx++){
            p1[indA] += pij[indA][indx];
        }
    }
    // For other positions | Los de las demas posiciones
    for(int indA=naa;indA<naanpos;indA++){
        for(int indx=0;indx<naa;indx++){
            p1[indA] += pij[indx][indA];
        }
    }

    //b) Compare Pij and Pi (from the model) to fij and fi (from data).
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
    memory_change_pi = memory_drag * (memory_change_pi) + (1.0 - memory_drag) * fabs(mp - maxPi);
    maxPi = mp;

    bool error = true;
    if(numHiNonConverged == 0 && numJijNonConverged == 0 && memory_change_pij < tolPij && memory_change_pi < tolPi){
        error = false;
    }//ACH 20180201 + AR20181101

    if(writeFiles){
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

    //c) Redefine parameters values (hi and Jij), if error is still large.
    //The Jij are defined according to L1 regularization.
    if(error){
	numJijNonzero = N_SITES * (N_SITES - 1) * 4 * 4;
        numHiNonConverged  = 0;                               //ACH 20180201
        numJijNonConverged = 0;

        for(int i=0;i<naanpos;i++){
            change = epsilonSingle * (f1Data[i] - p1[i]);// / (tolNum + f1Data[i] * (1.0 - f1Data[i]));//Segunda parte (p * (1-p)) agregada BK20190411. Agrego tolNum al denominador para que el maximo cambio absoluto no sea infinito.	     //ACH 20180201
            if(abs(change) > abs(hi[i] * TOLERANCE_HI) + tolabsHi){//Cambie de dividir por hi[i] el lado izquierdo, a multiplicar el lado derecho para evitar infinitos que ocurrian. BK20190411.
                //printf("%d %f %f %f %f", i, f1Data[i], p1[i], fabs(change), TOLERANCE_HI * hi[i]);//Para debuggear.
                numHiNonConverged++; //ACH 20180201
	    }
            hi[i]+=change;
            //printf("\n");//Para debuggear.
        }

        double auxJIJ;
	std::vector<std::vector<T>> prevJij = jij;
        int n1, n2;

        for(int fila=0; fila<naanpos; fila++){
            n1 = fila / naa;
            for(int col=0; col<naanpos; col++){
                n2 = col / naa;
                if(n1 != n2){
		    if(f2Data[fila][col] < pij[fila][col]){//En este bloque if se define el signo de la diferencia.
                        sign_change = -1;
                    }
                    else{
                        sign_change = 1;
                    }

                    if(jij[fila][col] == 0){//Si es cero.
                        if(abs(f2Data[fila][col] - pij[fila][col]) < gamma){//Caso en que no se despega del cero.
                            jij[fila][col] = 0;
			    numJijNonzero--;
                        }
                        else{
			    change = epsilonJoint * ((f2Data[fila][col] - pij[fila][col]) - gamma * sign_change);
                            if(abs(change) > fabs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){//Caso en que si se despega del cero. Mismos cambios que en H a la ecuacion de tolerancia. BK20190412
                                numJijNonConverged++; //ACH 20180201
			    }
                            jij[fila][col] =+ change;
                            //printf("jij %f pij %f fij %f\n", jij[fila][col], pij[fila][col], f2Data[fila][col]);//Para debuggeo.
                        }
                    }
                    else{//Si no es cero.
                    	if(jij[fila][col] < 0){//En este bloque if se define el signo del jij.
                            sign_j = -1;
                        }
                    	else{
                            sign_j = 1;
                   	}

			change = epsilonJoint * (f2Data[fila][col] - pij[fila][col] - gamma * sign_j);
                        auxJIJ = jij[fila][col] + change;
                        if(auxJIJ * jij[fila][col] < 0){//Si cruzaria la linea del cero.
                            if(abs(jij[fila][col]) > abs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){//Cambia en cantidad jij[][] al irse desde ese valor a 0.
                                numJijNonConverged++; //ACH 20180201
			    }
			    jij[fila][col] = 0;//Se pega al cero.
			    numJijNonzero--;
                            //printf("jij %f pij %f fij %f\n", jij[fila][col], pij[fila][col], f2Data[fila][col]);//Para debuggeo.
                        }
                        else{//Si no cruzaria la linea del cero.
			    if(sign_change * sign_j > 0){
				change = epsilonJoint * (f2Data[fila][col] - pij[fila][col]);//Crece de forma normal.
                            }
			    else{
				change = epsilonJoint * (f2Data[fila][col] - pij[fila][col]);//Amplifico los cambios negativos con la intensidad del jij.
			    }
			    if(abs(change) > abs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){
                                numJijNonConverged++; //ACH 20180201
			    }
                            jij[fila][col] += change;
                            printf("jij %f pij %f fij %f exp %f change %f\n", jij[fila][col], pij[fila][col], f2Data[fila][col], pij[fila][col] - f2Data[fila][col], change);//Para debuggeo.
                        }
                    }
                }
                else{
                    jij[fila][col] = 0;//Me aseguro que los de un mismo sitio sean nulos.
                }
            }
        }

	//Regularizacion de cantidad total de Jij.
        /*T totalJijIntensity = 0;//Calculo que tanta intensidad total tienen los Jij como medida de cuanto deberia estar controlandolos.
        for(int i=0; i<naanpos; i++){
            for(int j=0; j<naanpos; j++){
                totalJijIntensity += abs(jij[i][j]);
            }
        }
	printf("%f\n", totalJijIntensity);

	float total_crit = (1.0 - gamma) * ((double) numJijNonzero);
	if(totalJijIntensity > total_crit){//Asi controlo que no se pasen de cierta cantidad. Cuanto mas permisivo el gamma, mas puede ser esto. El 0.01 es arbitrario.
		numJijNonConverged = 0;
		for(int fila=0; fila<naanpos; fila++){
			for(int col=0; col<naanpos; col++){
				jij[fila][col] *= total_crit / totalJijIntensity;
				if(abs(jij[fila][col] - prevJij[fila][col]) > abs(jij[fila][col]) * TOLERANCE_JIJ + tolabsJij){
					numJijNonConverged++;
				}
			}
		}
	}*/

        // Mido cuanto tardo
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        cout << "... Tiempo: " << duration/1000000.0L << endl;
    }

    fprintf(errorFile, "%f\t%d\t%d\t%d\t%f\t%f\n", gamma, numHiNonConverged, numJijNonConverged, numJijNonzero, maxPi, maxPij);
    printf("%f\t%d\t%d\t%d\t%f\t%f\n", gamma, numHiNonConverged, numJijNonConverged, numJijNonzero, maxPi, maxPij);//Para debuggeo.
    return error;
}

int main(int argc, char* argv[]){
    double maxPij = 0.00;
    double maxPi = 0.00;
    cout << "Here we go!" << endl;

    /* ----------------------------------------------------------------------------- */
    /* 1- LOAD DATA AND DEFINE PARAMETERS .......................................... */
    /* ----------------------------------------------------------------------------- */
    //Parametros generales  | General parameters.
    static std::vector<T> p1(naanpos);
    static std::vector<std::vector<T>> pij(naanpos, std::vector<T>(naanpos));

    int a = 0;//Auxiliary counting variable.
    //Guardo todas las combinaciones posibles de dos posiciones.
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

    ifstream j0File(filename_Jij); /*check*/
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
		    jij[i][j] = 0;//log(f2Data[i][j] / (f1Data[i] * f1Data[j]));
		}
		else{
		    jij[i][j] = 0;
		}
	    }
	}
    }
    j0File.close();

    ifstream h0File(filename_hi); /*check*/
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

    bool error=true;
    int counter = 0;
    //bool writeFiles = true;
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
        double local_gamma = 1.0;//La primera iteracion la hago con gamma sobreexigente para ajustar las hi.
        if(g >= 0){
            if(through_zero){//Si voy a traves del cero, bajo y luego subo.
                local_gamma = abs(gamma_start - ((double) g) * gamma_step);
            }
            else{//De otra forma es monotona la subida o bajada.
                if(gamma_start > gamma_end){
                    local_gamma = gamma_start - ((double) g) * gamma_step;
                }
                else{
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
        while(error || counter < min_iterations){//Exijo una cantidad minima de iteraciones para darle una oportunidad a los parametros a cambiar.
            counter++;
            error = fit_iteration<N,T>(counter, local_gamma, p1, pij, f1Data, f2Data, maxPi, maxPij, memory_change_pi, memory_change_pij, t1, numHiNonConverged, numJijNonConverged, numJijNonzero, errorFile);
	    if(counter > MAX_NO_CONVERGENCE && g != 0){//Si lleva demasiado sin converger, el gamma era muy bajo.
		if(through_zero){//Si el objetivo era bajar hasta un punto y luego volver.
		    g = (int)((local_gamma + gamma_start)/gamma_step) + 1;//Lo devuelvo al paso anterior, en el momento en que toca subir en vez de bajar.
		}
		else{//Si el objetivo era solamente bajar.
		    g = n_gamma_steps;//Termino la secuencia.
		}
	    }
        }

        /* ----------------------------------------------------------------------------- */
        /* 3- SAVE FINAL DATA .......................................................... */
        /* ----------------------------------------------------------------------------- */

        /*if(writeFiles){
	    std::string appendstring = "gamma";
            appendstring.append(std::to_string(local_gamma));
            appendstring.append("_Final.txt");

            std::string fullnameHI = prefixHI;
            std::string fullnameJIJ = prefixJij;
            std::string fullnamePI = prefixPI;
            std::string fullnamePIJ = prefixPIJ;
            fullnameHI.append(appendstring);
            fullnameJIJ.append(appendstring);
            fullnamePI.append(appendstring);
            fullnamePIJ.append(appendstring);

            write_1(fullnameHI, hi);
            write_2(fullnameJIJ, jij);
            write_1(fullnamePI, p1);
            write_2(fullnamePIJ, pij);
        }*/

        fprintf(errorFile, "%f\t%d\t%d\t%d\t%f\t%f\n", local_gamma, numHiNonConverged, numJijNonConverged, numJijNonzero, maxPi, maxPij);
    }
    fclose(errorFile);

    return 0;
}
