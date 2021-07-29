//Escribo un vector a un ostream.
template<typename T>
ostream& operator << (ostream& out, const vector<T>& v){
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
    }
    return out;
}

//Cargo probabilidades (o parametros) de un sitio.
template<typename T>
void load_f1(std::string filename, std::vector<T> &f1Data){
    ifstream f1File(filename);  /*check*/
    for(int j=0; j<naanpos; j++){
        f1File >> f1Data[j];
	//std::cout << f1Data[j] << std::endl;
    }
    f1File.close();
}

//Cargo probabilidades (o parametros) de dos sitios.
template<typename T>
void load_f2(std::string filename, std::vector<std::vector<T>> &f2Data){
    ifstream f2File(filename);  /*check*/
    for(int fila=0;fila<naanpos;fila++){
        for(int col=0;col<naanpos;col++){
            f2File >> f2Data[fila][col];
	    //std::cout << f2Data[fila][col] << " ";
        }
	//std::cout << std::endl;
    }
    f2File.close();
}

//Escribo probabilidades (o parametros) de un sitio.
template<typename T>
void write_1(std::string fullname, std::vector<T> data){
    ofstream file(fullname);
    for (int j=0; j<naanpos; j++){
        file << data[j] << "\t";
        if((j+1) % naa == 0){
            file << endl;
        }
    }
    file.close();
}

//Escribo probabilidades (o parametros) de dos sitios.
template<typename T>
void write_2(std::string fullname, std::vector<std::vector<T>> data){
    ofstream file(fullname);
    for(int fila=0; fila<naanpos; fila++){
        for(int col=0; col<naanpos; col++){
            file << data[fila][col] << "\t";
        }
        file << endl;
    }
    file.close();
}
