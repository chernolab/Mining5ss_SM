//The template operator allows us to later choose any integer type and floating-point type for N and T respectively.
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

    //This function takes the initial sequence and, using the available fitted parameters, evolves it until it converges to the least-energy candidate, with some thermal noise.
    void termalizarSecuencia(T &energy, std::mt19937 &rng){//This makes sure that it's near the lowest possible energy, which is important for later producing a sequence ensemble.
        std::vector<N> sequence(secuencia_inicial);

        T dE;
        T r, expE;
        N mpos, mres, oldRes;
        N pasos_sin_disminuir = 0;
        T min_energy = energy;
        //Elijo mutacion al azar | random mutation.
        while(pasos_sin_disminuir < termalizacion_sin_cambio){
            mpos = uni65(rng);//The position for mutating.
            mres = uni20(rng);//The possible mutation (to be accepted or rejected).
            oldRes = sequence[mpos];//The old (current) value in that position.

            if(mres != oldRes){//If the site is indeed proposing a change.
                dE = -hi[naa * mpos + mres] + hi[naa * mpos + oldRes];//I first change local energy.

		            for(int i=0; i<N_SITES; i++){
                    if(i != mpos){
			                dE -= jij[naa * mpos + mres][naa * i + sequence[i]];
                	    dE += jij[naa * mpos + oldRes][naa * i + sequence[i]];
		                }
		            }

                //Decide whether I keep the mutation or not (rejection sampling).
                if(dE < 0){
                    sequence[mpos] = mres;
                    energy += dE;
                }
                else{
                    r = dist(rng);//This is a more uniform alternative to the typical ((float)rand())/((float)MAX_RAND) approach.
                    expE = exp(-dE);
                    if(r < expE){
                        sequence[mpos] = mres;
                        energy += dE;
                    }
                }
            }

            if(min_energy > energy){
                min_energy = energy;
                pasos_sin_disminuir = 0;
            }
            else{
                pasos_sin_disminuir++;
            }
        }

        secuencia_inicial = sequence;//The sequence is finally changed.
    }

    //This is the operator used to generate a given sequence by randomly evolving from the current one.
    //It will implement a similar procedure to the above function.
    std::vector<N> operator()(T &energy, int rng_seed){//It is set up as an operator to be used in either parallel or sequentially as determined by the compiler.
        std::vector<N> sequence(secuencia_inicial);
        std::mt19937 rng(rng_seed);//Random-number engine used (Mersenne-Twister in this case).

        T dE;
        T r, expE;
        N mpos, mres, oldRes;
        //Random mutation, similar to above.
        for(int iteracion=0; iteracion<cantidad_de_mutaciones; iteracion++){
            mpos = uni65(rng);
            mres = uni20(rng);
            oldRes = sequence[mpos];

            if(mres != oldRes){
                dE = -hi[naa * mpos + mres] + hi[naa * mpos + oldRes];
              
                for(int i=0; i<N_SITES; i++){
		              if(i != mpos){
                	  dE -= jij[naa * mpos + mres][naa * i + sequence[i]];
                	  dE += jij[naa * mpos + oldRes][naa * i + sequence[i]];
		              }
                }

                if(dE < 0){
                    sequence[mpos] = mres;
                    energy += dE;
                }
                else{
                    r = dist(rng);
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

    //Its initial energy.
    T energy = 0.0;
    for(N i=0; i<N_SITES; i++){
        energy += -hi[naa*i+sequence[i]];
        if(i!=(N_SITES-1)){
            for(int j=(i+1); j<N_SITES; j++){
                energy += -jij[naa*i+sequence[i]][naa*j+sequence[j]];
            }
        }
    }

    GenerarSecuencia<N,T> seq_generator(sequence, cantidad_de_mutaciones);//Sequence evolution functor.
    seq_generator.termalizarSecuencia(energy, rng);//Converge the main sequence close to minimum energy.

    std::vector<int> rng_seeds(cantidad_de_secuencias, rd());//Initialize a random seed vector used for different sequences to evolve according to different randomness.
    std::vector<int> counting_vector(cantidad_de_secuencias);
    std::iota(counting_vector.begin(), counting_vector.end(), 0);//A vector that counts between 0 and the number of sequences minus 1.
    std::transform(counting_vector.begin(), counting_vector.end(), rng_seeds.begin(), rng_seeds.begin(), std::plus<int>());//Adding both, I get a random vector where all entries are different.

    std::fill(energies.begin(), energies.end(), energy);//Initialize all sequence energies at the last sequence energy recorded.
    std::transform(energies.begin(), energies.end(), rng_seeds.begin(), sequences.begin(), seq_generator);//Mutate different sequences from this initial sequence, using the random seed vector previously created.

    std::transform(sequences.begin(), sequences.end(), energies.begin(), calcEnergies<N,T>());//Calculate their individual energies.
}
