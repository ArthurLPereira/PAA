#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>

using namespace std;

ofstream file

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
        file << x << " " << y << endl;
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


int main() {
    int n;
    cin >> n;
    file.open("bruteforce.txt");
    for(int j = 0; j < 15; j++){
        Graph graph = construirGrafo(n);
        vector<int> resposta;
        double caminhoTotal = 0;

        resposta = graph.bruteForce();
        string caminho = "";

        for (int i = 0; i < resposta.size(); i++) {
            
            caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
            caminho+= to_string(resposta[i] + 1) + " ";
        }

        cout << caminhoTotal << endl << caminho << endl << endl;

    }
    file.close()
	//cout << "Tempo BF: " << graph.timeT << endl;
}