#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <ctime>
#include <list>
#include <utility>
#include <map>
#include <set>
#include <stdlib.h>

using namespace std;


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

    Graph graph(n);

    for (int i = 0; i < n; i++){
       // x = distr(eng);
       // y = distr(eng);
       cin >> x;
       cin >> y;
       // file << x << " " << y << endl;
       graph.addCoor(i, x, y);
    }

    //calcula distancia entre os verticies, partindo das coordenadas.
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i != j){
                //distancia euclidiana do grafo,
                //nas coordenadas (i,0), (i,1), (j,0) e (j,1).
                aux = distanciaEuclidiana(graph.coor[i][0], graph.coor[i][1], graph.coor[j][0], graph.coor[j][1]);

                // adiciona aresta no grafo, na posica i, j com o valor de aux.
                graph.addEdge(i, j, aux);
            }
        }
    }
    return graph;
}

/* metodo main aonde é feito chamada do metodo genetico.
calcula tambem o caminho total percorrido.

*/
int main(){
   //cria arquivo com o output do genetico.
   ofstream file;
   file.open("genetic.txt");
   int n;
   cin >> n;

   //for (int j = 0; j < 15; j++){
        Graph graph = construirGrafo(n); //criando o grafo.
        vector<int> resposta;    //vetor que vai armazenar a ordem percorrida.
        double caminhoTotal = 0; //contem caminho total percorrido.

        algoritmoGenetico ag(&graph, 40, 100, 20, 0);

        ag.gerarMenorCusto();
        string caminho = ""; //melhor caminho percorrido

         /*
        //for para cacular distancia total percorrida,
        //fazer a impressao na tela da ordem dos vertices percorridos.
        for (int i = 0; i < resposta.size(); i++){
            caminhoTotal += grafo.adj[resposta[i]][resposta[i + 1]]; //faz a soma da distancia do i-vertcie visitado e o seguinte.
            caminho += to_string(resposta[i] + 1) + " ";
        }
        file << caminhoTotal << endl; //print da distancia total percorrida
        file << caminho << endl << endl;      //printa o caminho percorrido. */
 //print da distancia total percorrida e do menor caminho

        //cout << "Tempo BF: " << graph.timeT << endl;
    //}
    file.close(); //fecha arq
}