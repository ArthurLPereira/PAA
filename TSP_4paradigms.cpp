// #include "stdafx.h"
#include "Graph.h"
#include "math.h"
#include <iostream>
#include <vector>

using namespace std;

double distanciaEuclidiana(int x, int y, int a, int b) {

	return sqrt(pow(x - a, 2) + pow(y - b, 2));

}//fim distancia euclidiana
/*
4
100 100
900 100
900 900
100 900
*/
Graph construirGrafo() {

	int n; //qtd cidades
	double x, y; //coordenadas
	double aux; //variavel auxiliar 

	cin >> n;

	Graph graph(n);

	for (int i = 0; i < n; i++) {

		cin >> x;
		cin >> y;

		graph.addCoor(i, x, y);

	}//fim for

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			if (i != j) {

				aux = distanciaEuclidiana(graph.coor[i][0], graph.coor[i][1], graph.coor[j][0], graph.coor[j][1]);

				graph.addEdge(i, j, aux);

			}//fim if

		}//fim for

	}//fim for

	return graph;


}//fim construirGrafo


int main() {

	Graph graph = construirGrafo();
	vector<int> resposta;
	double caminhoTotal;

	/*Brute Force*/

	resposta = graph.bruteForce();

	//cout << resposta.size() << endl;

	caminhoTotal = graph.adj[resposta[0]][resposta[resposta.size() - 1]];

	for (int i = 0; i < resposta.size() - 1; i++) {

		caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
		cout << resposta[i] + 1 << " ";

	}//fim for

	cout << resposta[resposta.size() - 1] + 1 << endl;

	cout << caminhoTotal << endl;

	//cout << "Tempo BF: " << graph.timeT << endl;

	/*Branch and Bound*/

	resposta = graph.branchBound();

	//cout << resposta.size() << endl;

	caminhoTotal = graph.adj[resposta[0]][resposta[resposta.size() - 1]];

	for (int i = 0; i < resposta.size() - 1; i++) {

		caminhoTotal += graph.adj[resposta[i]][resposta[i + 1]];
		cout << resposta[i] + 1 << " ";

	}//fim for

	cout << resposta[resposta.size() - 1] + 1 << endl;

	cout << caminhoTotal << endl;


	cin.get();
	cin.get();


	/*debbug

	cout << " adj: " << endl;
	graph.printAdj();

	cout << "coor: " << endl;
	graph.printCoor();



	cin.get();
	cin.get();
	*/
}//fim main


