#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>

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


int brute_force() {
	
	for(int n = 4; n < 15; n++) {
		file.open("bruteforce/N"+to_string(n)+".txt");
		for(int j = 0; j < 15; j++){
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
			file << graph.timeT << endl;

		}
        file.close();

	}
}

int branch_bound() {
	
	for(int n = 4; n < 15; n++) {
		file.open("branchbound/N"+to_string(n)+".txt");
		for(int j = 0; j < 15; j++){
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
			file << graph.timeT << endl;

		}
        file.close();

	}
}

int dynamic() {
	
	for(int n = 4; n < 15; n++) {
		file.open("dynamic/N"+to_string(n)+".txt");
		for(int j = 0; j < 15; j++){
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
			file << graph.timeT << endl;

		}
        file.close();

	}
}

int main() {
	
	brute_force();
    branch_bound();
    dynamic();

}