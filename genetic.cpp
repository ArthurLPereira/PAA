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

ofstream file;

double distanciaEuclidiana(int x, int y, int a, int b) {
	return sqrt(pow(x - a, 2) + pow(y - b, 2));
}

Graph construirGrafo(int n) {

	random_device rd;
	mt19937 eng(rd());
	uniform_int_distribution<> distr(0, 9999);

	double x, y; //coordenadas
	double aux; //variavel auxiliar 

	Graph graph(n);

	for (int i = 0; i < n; i++) {

		x = distr(eng);
		y = distr(eng);
        // file << x << " " << y << endl;
		graph.addCoor(i, x, y);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) {

				aux = distanciaEuclidiana(graph.coor[i][0], graph.coor[i][1], graph.coor[j][0], graph.coor[j][1]);
				graph.addEdge(i, j, aux);
			}
		}
    }

	return graph;
}

void algoritmoGenetico::iniciarPopulacao(){
      vector<double> pai;
      //Insere o vertice inicial no vetor pai
      pai.push_back(verticeInicial);
      //Cria o vetor pai
      for (int i = 0; i < grafo->n; i++){
         if(i != verticeInicial){
            pai.push_back(i);
         }
      }

      double custoTotal =  custoCaminho(pai);
      //Verifica se a rota contem ligacoes e uma rota valida, se for insere na populacao e incremeta o contador
      if(custoTotal != -1){
         populacao.push_back(make_pair(pai,custoTotal));
         tamanhoRealPopulacao++;
      }

      //Cria rotas aleatorias, posteriormente verifica se tem custo maior que zero
      //e um crossomo ja inserido e faz a insercao na populacao
      for (int i = 0; i < geracoes; i++){
         random_shuffle(pai.begin() + 1, pai.begin() + (rand() % (grafo->n - 1) + 1));
         double custoTotal = custoCaminho(pai);
         if(custoTotal != -1 && !existeCromossomos(pai)){
            populacao.push_back(make_pair(pai,custoTotal));
            tamanhoRealPopulacao++;
         }
         if(tamanhoPopulacao == tamanhoRealPopulacao)
            break;
      }

      if(tamanhoRealPopulacao == 0)
         cout << "\nPopulacao inicial vazia"<<endl;
      else
         sort(populacao.begin(),populacao.end(),ordenarPredecessor());
   }




int main() {
	file.open("branchbound.txt");
	int n;
	cin >> n;

	for(int j = 0; j < 15; j++) {

		Graph graph = construirGrafo(n);
		vector<int> resposta;
		double caminhoTotal = 0;

        algoritmoGenetico ag(graph,40,100,20,0);

		ag.gerarMenorCusto();
		string caminho = ""; // melhor caminho percorrido

		// for (int i = 0; i < resposta.size(); i++) {
			
		// 	caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
		// 	caminho+= to_string(resposta[i] + 1) + " ";
		// }
        caminhoTotal = ag.menorCusto();

		file << caminhoTotal << endl << caminho << endl << endl;

		//cout << "Tempo BF: " << graph.timeT << endl;
	}
	file.close();
}