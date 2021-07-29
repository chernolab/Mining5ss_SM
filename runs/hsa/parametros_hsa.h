//ATENCION: Si se va a hacer una corrida pero no se desea borrar la anterior, CAMBIAR el directorio de guardado prefix_directory y filename_errorfile.
//PARAMETROS DE MODELO
#define N_SITES 9//Cantidad de sitios en las secuencias.
N termalizacion_sin_cambio = 1000;//Numero de pasos luego de encontrar la menor energia que se debe cumplir para considerar termalizado al sistema.

//PARAMETROS DE AJUSTE
#define MAX_NO_CONVERGENCE 100000000

#define epsilonSingle 0.2// epsilon to renew parameters hi
#define epsilonJoint 0.2// epsilon to renew parameters jij

#define tolNum 0.0000001//Tolerancia numerica maxima para no tener ceros en denominadores.

//Criterio de convergencia parametros: cambio < TOL * valor_actual + tolabs.
#define TOLERANCE_HI 0.005//Tolerancias fraccionales para el cambio de los parametros.
#define TOLERANCE_JIJ 0.01
#define tolabsHi 0.000001//Tolerancias absolutas para el cambio de los parametros.
#define tolabsJij 0.000001

//Criterio de convergencia en probabilidades: que el error respecto a las observaciones no disminuya.
//Esto se hace con un promedio pesado entre el error actual y el error anterior, para tomar en cuenta como viene disminuyendo.
//memoria_error = memoria_error_previo * memory_drag + (1.0 - memory_drag) * error_actual
#define memory_drag 0.9//Cuanto peso se le da al pasado.
#define tolPi 0.001//Cuan chico debe ser la memoria del cambio en el error para que haya convergido.
#define tolPij 0.001

#define SavingDataSteps 10//Cada cuanto se guardan todos los datos.
#define SaveSequencesRate 10//Cada cuanto se guardan ademas las secuencias y sus energias.

#define TOTAL_SEQUENCES 100000//Cantidad de secuencias a generar cuando se usa el modelo.
#define MUTATION_STEPS 10000

//PARAMETROS DE ARCHIVOS Y DIRECTORIOS
std::string filename_marginales = "archivos/humano/fi.txt";//Probabilidades de un sitio.
std::string filename_dobles = "archivos/humano/fij.txt";//Probabilidades de dos sitios.

std::string filename_Jij = ""; // corridas/humano_u12_invertida/gamma3.000000/ParamJ1231";//Si existen, estimaciones anteriores de parametros J.
std::string filename_hi =  ""; // corridas/humano_u12_invertida/gamma3.000000/ParamH1231";//Lo mismo para los h.

std::string prefixDirectory = "corridas/humanoU2_ach/";//El directorio de las corridas.

/*std::string prefixHI = "corridas/arabidopsis/ParamH";//Como comenzar el nombre de cada archivo de guardado.
std::string prefixJij="corridas/arabidopsis/ParamJ";
std::string prefixPI="corridas/arabidopsis/Pi";
std::string prefixPIJ="corridas/arabidopsis/Pij";
std::string prefixSequenceFile="corridas/arabidopsis/sequences";
std::string prefixEnergyFile="corridas/arabidopsis/Energy";*/
std::string filename_errorfile="corridas/humanoU2_ach/error.txt";

//PARAMETROS DE BARRIDO
#define gamma_start 0.05//En que gamma empiezo.
#define gamma_end 0.01//En que gamma termino.
#define gamma_step 0.005//Por cuanto voy barriendo.
bool through_zero = false;//Si bajo hasta el cero y vuelvo, o no.
