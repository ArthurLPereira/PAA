#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>

using namespace std;

ofstream file;

/*
Calcula a distancia entre dois pontos
@author		Ricardo Xavier Sena
@author     Arthur Pereira

@param	x	Coordenada x do ponto1
@param	y	Coordenada y do ponto1
@param	a	Coordenada x do ponto2
@param	b	Coordenada y do ponto2

@return		Distancia, entre os dois pontos
*/
double distanciaEuclidiana(int x, int y, int a, int b){
    double resposta;

    //distancia euclidiana
    resposta = sqrt(pow(x - a, 2) + pow(y - b, 2));

    return resposta;
}

Graph construirGrafo(int n){
    random_device rd;
    mt19937 eng(rd());
    uniform_int_distribution<> distr(0, 9999);

    int x, y; //coordenadas x e y
    double aux;  //variavel aux

    Graph grafo(n);

    for (int i = 0; i < n; i++)
    {
        x = distr(eng);
        y = distr(eng);
        //file << x << " " << y << endl;
        //adiciona as coordenadas em uma matriz.
        grafo.addCoor(i, x, y);
    }

    //calcula distancia entre os verticies, partindo das coordenadas.
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i != j){
                //distancia euclidiana do grafo,
                //nas coordenadas (i,0), (i,1), (j,0) e (j,1).
                aux = distanciaEuclidiana(grafo.coor[i][0], grafo.coor[i][1], grafo.coor[j][0], grafo.coor[j][1]);

                // adiciona aresta no grafo, na posica i, j com o valor de aux.
                grafo.addEdge(i, j, aux);
            }
        }
    }
    return grafo;
}

/* metodo main aonde é feito chamada do metodo de branch&bound.
calcula tambem o caminho total percorrido.

*/
int main(){
    //cria arquivo com o output do branch and borund.
    file.open("branch-bound.txt");

    int n;
    cin >> n;
   

    for (int j = 0; j < 15; j++){
        vector<int> resposta;    //vetor que vai armazenar a ordem percorrida.
        double caminhoTotal = 0; //contem caminho total percorrido.

        //criando o grafo.
        Graph grafo = construirGrafo(n);

        //fazendo o branchBound.
        resposta = grafo.branchBound();
        string caminho = ""; //melhor caminho percorrido

        //for para cacular distancia total percorrida,
        //fazer a impressao na tela da ordem dos vertices percorridos.
        for (int i = 0; i < resposta.size(); i++){
            caminhoTotal += grafo.adj[resposta[i]][resposta[i + 1]]; //faz a soma da distancia do i-vertcie visitado e o seguinte.
            caminho += to_string(resposta[i] + 1) + " ";
        }
        file << caminhoTotal << endl; //print da distancia total percorrida
        file << caminho << endl << endl;      //printa o caminho percorrido.
    }
    file.close(); //fecha arq
}