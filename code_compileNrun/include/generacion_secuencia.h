//Uso un template para poder elegir cualquier tipo de entero N y cualquier tipo de flotante T.
template<typename N, typename T>
struct GenerarSecuencia{
    private:
    N cantidad_de_mutaciones;
    std::vector<N> secuencia_inicial;

    public:
    GenerarSecuencia(std::vector<N> input_secuencia_inicial, N input_cantidad_de_mutaciones){
        secuencia_inicial = input_secuencia_inicial;
        cantidad_de_mutaciones = input_cantidad_de_mutaciones;
    }

    void termalizarSecuencia(T &energy, std::mt19937 &rng){
        std::vector<N> sequence(secuencia_inicial);

        T dE;
        T r, expE;
        N mpos, mres, oldRes;
        N pasos_sin_disminuir = 0;
        T min_energy = energy;
        //Elijo mutacion al azar | random mutation.
        while(pasos_sin_disminuir < termalizacion_sin_cambio){
            mpos = uni65(rng);//La posicion.
            mres = uni20(rng);//La mutacion candidato.
            oldRes = sequence[mpos];//El valor viejo en la posicion.

            if(mres != oldRes){//Si hay un cambio.
                dE = -hi[naa * mpos + mres] + hi[naa * mpos + oldRes];//Primero hago el cambio en la energia local.

		for(int i=0; i<N_SITES; i++){
                    if(i != mpos){
			dE -= jij[naa * mpos + mres][naa * i + sequence[i]];
                	dE += jij[naa * mpos + oldRes][naa * i + sequence[i]];
		    }
		}

                //Elijo si me quedo o no con la mutacion | Decide whether I keep the mutation or not.
                if(dE < 0){
                    sequence[mpos] = mres;
                    energy += dE;
		    //printf("        E: %f   |   dE: %f\n", energy, dE);
                }
                else{
                    r = dist(rng);//Esto es más uniforme pero más lento.
                    expE = exp(-dE);
                    if(r < expE){
                        sequence[mpos] = mres;
                        energy += dE;
			//printf("	E: %f   |   dE: %f\n", energy, dE);
			/*printf("	");
			for(int c=0; c<9; c++){
			   printf("%d", sequence[c]);
			}
			printf("\n");*/
                    }
                }
            }

            if(min_energy > energy){
                min_energy = energy;
                pasos_sin_disminuir = 0;
		//printf("	AAA: %f\n", min_energy);
            }
            else{
                pasos_sin_disminuir++;
            }
            //printf("%f, %d\n", min_energy, pasos_sin_disminuir);//Test para ver cuanto tarda en termalizar.
        }

        secuencia_inicial = sequence;
    }

    std::vector<N> operator()(T &energy, int rng_seed){
        std::vector<N> sequence(secuencia_inicial);
        std::mt19937 rng(rng_seed);//Random-number engine used (Mersenne-Twister in this case).

        T dE;
        T r, expE;
        N mpos, mres, oldRes;
        //Elijo mutacion al azar | random mutation.
        for(int iteracion=0; iteracion<cantidad_de_mutaciones; iteracion++){
            mpos = uni65(rng);//La posicion.
            mres = uni20(rng);//La mutacion candidato.
            oldRes = sequence[mpos];//El valor viejo en la posicion.

            if(mres != oldRes){//Si hay un cambio.
                dE = -hi[naa * mpos + mres] + hi[naa * mpos + oldRes];//Primero hago el cambio en la energia local.

                for(int i=0; i<N_SITES; i++){
		    if(i != mpos){
                	dE -= jij[naa * mpos + mres][naa * i + sequence[i]];
                	dE += jij[naa * mpos + oldRes][naa * i + sequence[i]];
		    }
                }


                //Elijo si me quedo o no con la mutacion | Decide whether I keep the mutation or not.
                if(dE < 0){
                    sequence[mpos] = mres;
                    energy += dE;
                }
                else{
                    r = dist(rng);//Esto es más uniforme pero más lento.
                    expE = exp(-dE);
                    if(r < expE){
                        sequence[mpos] = mres;
                        energy += dE;
                    }
                }
            }
        }

        return sequence;
    }
};

template<typename N, typename T>
struct calcEnergies{
	T operator()(std::vector<N> &sequence){
		T energy = 0.0;
    		for(N i=0; i<N_SITES; i++){
        		energy += -hi[naa*i+sequence[i]];
        		if(i!=(N_SITES-1)){
            			for(int j=(i+1); j<N_SITES; j++){
                			energy += -jij[naa*i+sequence[i]][naa*j+sequence[j]];
            			}
        		}
		}
		return energy;
	}
};

template<typename N, typename T>
void generar_secuencias(N cantidad_de_secuencias, N cantidad_de_mutaciones, std::vector<std::vector<N>> &sequences, std::vector<T> &energies){
    std::random_device rd;
    std::mt19937 rng(rd());//Random-number engine used (Mersenne-Twister in this case).

    //Random initial sequence.
    vector<N> sequence;
    for(N i=0; i<N_SITES; i++){
        N random_integer = uni20(rng);
        sequence.push_back(random_integer);
    }

    //Su energia inicial. | Its initial energy.
    T energy = 0.0;
    for(N i=0; i<N_SITES; i++){
        energy += -hi[naa*i+sequence[i]];
        if(i!=(N_SITES-1)){
            for(int j=(i+1); j<N_SITES; j++){
                energy += -jij[naa*i+sequence[i]][naa*j+sequence[j]];
            }
        }
    }
    //printf("ENERGY: %f\n", energy);

    GenerarSecuencia<N,T> seq_generator(sequence, cantidad_de_mutaciones);//Armo el functor de avance de secuencias.
    seq_generator.termalizarSecuencia(energy, rng);//Termalizo la secuencia principal.

    std::vector<int> rng_seeds(cantidad_de_secuencias, rd());//Inicializo el vector lleno de un seed que saco de rd.
    std::vector<int> counting_vector(cantidad_de_secuencias);
    std::iota(counting_vector.begin(), counting_vector.end(), 0);//Un vector que cuenta entre 0 y cantidad_de_secuencias - 1.
    std::transform(counting_vector.begin(), counting_vector.end(), rng_seeds.begin(), rng_seeds.begin(), std::plus<int>());//Sumando los dos, obtengo un vector aleatorio por rd() donde todas las entradas son diferentes.

    std::fill(energies.begin(), energies.end(), energy);
    std::transform(energies.begin(), energies.end(), rng_seeds.begin(), sequences.begin(), seq_generator);//Hago mutaciones.

    std::transform(sequences.begin(), sequences.end(), energies.begin(), calcEnergies<N,T>());
}
