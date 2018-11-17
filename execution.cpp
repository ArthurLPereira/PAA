#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <ctime>
#include <time.h>
#include <list>
#include <utility>
#include <map>
#include <set>
#include <stdlib.h>

#define MAX_NODES 12
#define TEST_CASES 15

using namespace std;

ifstream input;	
ofstream output;

double getTime(clock_t x) {
	return ((double)x)/CLOCKS_PER_SEC;
}

double distanciaEuclidiana(int x, int y, int a, int b) {

	return sqrt(pow(x - a, 2) + pow(y - b, 2));

}

void gerar_amostra(){
	random_device rd;
	mt19937 eng(rd());
	uniform_int_distribution<> distr(0, 9999);
	double x, y;
	output.open("input.txt");
	for(int n = 4; n <= MAX_NODES; n++){
		for(int i = 0; i < TEST_CASES; i++) {
			output << n << endl;
			for(int k = 0; k < n; k++){
				x = distr(eng);
				y = distr(eng);
				output << x << " " << y << endl;		
			}
		}
	}
	output.close();
}

Graph construirGrafo(int n) {

	double x, y; //coordenadas
	double aux;  //variavel auxiliar 

	Graph graph(n);

	for (int i = 0; i < n; i++) {

		input >> x;
		input >> y;
        // output << x << " " << y << endl;
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


int brute_force() {
	input.open("input.txt");
	for(int n = 4; n <= MAX_NODES; n++) {
		output.open("bruteforce/N"+to_string(n)+".txt");
		for(int j = 0; j < TEST_CASES; j++){
			input >> n;
			Graph graph = construirGrafo(n);
			vector<int> resposta;
			double caminhoTotal = 0;

			resposta = graph.bruteForce();
			//string caminho = "";

			// for (int i = 0; i < resposta.size(); i++) {
				
			// 	caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
			// 	caminho+= to_string(resposta[i] + 1) + " ";
			// }

			// cout << caminhoTotal << endl << caminho << endl << endl;
			output << getTime(graph.timeT) << endl;

		}
        output.close();

	}
	input.close();
}

int branch_bound() {
	input.open("input.txt");
	for(int n = 4; n <= MAX_NODES; n++) {
		output.open("branchbound/N"+to_string(n)+".txt");
		for(int j = 0; j < TEST_CASES; j++){
			input >> n;
			Graph graph = construirGrafo(n);
			vector<int> resposta;
			double caminhoTotal = 0;

			resposta = graph.branchBound();
			//string caminho = "";

			// for (int i = 0; i < resposta.size(); i++) {
				
			// 	caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
			// 	caminho+= to_string(resposta[i] + 1) + " ";
			// }

			// cout << caminhoTotal << endl << caminho << endl << endl;
			output << getTime(graph.timeT) << endl;
		}
        output.close();

	}
	input.close();
}

int dynamic() {
	input.open("input.txt");
	for(int n = 4; n <= MAX_NODES; n++) {
		output.open("dynamic/N"+to_string(n)+".txt");
		for(int j = 0; j < TEST_CASES; j++){
			input >> n;
			Graph graph = construirGrafo(n);
			vector<int> resposta;
			double caminhoTotal = 0;

			resposta = graph.dynamic();
			//string caminho = "";

			// for (int i = 0; i < resposta.size(); i++) {
				
			// 	caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
			// 	caminho+= to_string(resposta[i] + 1) + " ";
			// }

			// cout << caminhoTotal << endl << caminho << endl << endl;
			output << getTime(graph.timeT) << endl;
		}
        output.close();

	}
	input.close();
}

int genetic() {
	input.open("input.txt");
	for (int n = 4; n <= MAX_NODES; n++){
		output.open("genetic/N"+to_string(n)+".txt");

		for(int j = 0; j < TEST_CASES; j++) {
			input >> n;
			Graph graph = construirGrafo(n);
			vector<int> resposta;
			double caminhoTotal = 0;

			algoritmoGenetico ag(&graph,40,100,20,0);

			clock_t inicio, fim;
			inicio = clock();
			ag.gerarMenorCusto();
			fim = clock();
			clock_t timeT = fim - inicio;
			// melhor caminho percorrido

			// for (int i = 0; i < resposta.size(); i++) {
				
			// 	caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
			// 	caminho+= to_string(resposta[i] + 1) + " ";
			// }

			output << getTime(timeT) << endl;
			//cout << "Tempo BF: " << getTime(graph.timeT) << endl;
		}
		output.close();
	}
	input.close();
}

int main() {
	gerar_amostra();
	brute_force();
    branch_bound();
    dynamic();
	genetic();
}