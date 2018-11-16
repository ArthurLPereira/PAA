#include <iostream>
#include <vector>
#include <ctime>
#define INFINITE 999999/0.001

using namespace std;

class Graph {

private:

	int n;	// Number of vertices
	vector<int> bruteForceR(int a, double res, double bestR,  vector<int>cidades, vector<int> bestC); //determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de for�a bruta
	vector<int> branchBoundR(int a, double res, double bestR, vector<int>cidades, vector<int> bestC); //determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de branch and bound 
	bool isIn(int x, vector<int> v); //retorna se um elemento esta ou nao no vetor de inteiros
	double best; //

public:

	time_t timeT;
	double **adj;	// matriz de adjacência/distancia. Se o valor for -1, os vértices não estão conectados
	double **coor;  // coordenadas das cidades
	Graph();  // construtor nao utilizado
	Graph(int size);	// Constructor
	~Graph();	// Destructor
	int getSize(); // recupera o tamanho do grafo
	void addCoor(int pos, double x, double y); //adiciona coordendas da cidade
	void addEdge(int x, int y, double wt);  // Add an edge to the graph
	bool isConnected(int x, int y);	// Check if two vertices are connected
	vector<int> bruteForce(); //determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de for�a bruta (chama a fun��o recursiva)
	vector<int> branchBound(); //determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de branch and bound (chama a fun��o recursiva)
	vector<int> dynamic(); //determina o menor caminho Hamiltoniano do grado utilizando programação dinamica
	void printAdj(); //imprime a matriz de adjacencia
	void printCoor(); //imprime matriz de coordenadas

};


// Constructor

Graph::Graph(int size) {

	if (size > 1) {

		n = size;

	}//fim if

	else {

		n = 2;

	}//fim else

	adj = new double*[n];
	coor = new double*[n];

	for (int i = 0; i < n; i++) {

		adj[i] = new double[n];	// Allocate memory for adjacency matrix
		coor[i] = new double[n];

		for (int j = 0; j < n; j++) {

			adj[i][j] = -1;	// Initialize the vertices to -1

		}//fim for

		for (int j = 0; j < 2; j++) {

			coor[i][j] = -1;	// Initialize the vertices to -1		

		}//fim for

	}//fim for

	this->best = INFINITE;
	this->timeT = 0;

}

//Constructor not utilized
Graph::Graph() {
	n = 0;
	this->best = INFINITE;
}


// Destructor not utilized
Graph::~Graph() {

}


int Graph::getSize() {

	return n;

}//fim getSize



void Graph::addCoor(int pos, double x, double y) {

	coor[pos][0] = x;
	coor[pos][1] = y;

}

void Graph::addEdge(int x, int y, double wt) {

	adj[x][y] = adj[y][x] = wt;

}

bool Graph::isConnected(int x, int y) {

	bool res = false;

	if (adj[x][y] >= 0) {

		res = true;

	}

	return res;

}

/*
Determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de forea bruta (chama a funcao recursiva)
@return	resposta	Vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho
@author				Henrique Schiess Pertussati
*/
vector<int> Graph::bruteForce() {

	time_t inicio, fim;
	vector<int> cidades;
	vector<int> resposta;
	cidades.push_back(0);
	this->best = INFINITE;

	inicio = clock();
	resposta = bruteForceR(0, 0, 0, cidades, cidades);
	fim = clock();
	timeT = fim - inicio;

	return resposta;
}

/*
Determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de forca bruta
@param	a		Vertice que estamos visitando
@param	res		Caminho total atual
@param	cidades	Vetor de inteiros com o caminho atual percorrido
@param	bestC	Atual vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho
@return	bestC	Vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho 	
@author			Arthur Ladislau Pereira
@author			Henrique Schiess Pertussati
*/
vector<int> Graph::bruteForceR(int a, double res, double bestR, vector<int>cidades,  vector<int> bestC) {

	for (int i = 0; i < n; i++) { // executa a permutacao as arestas do vertice corrente

		if (!isIn(i, cidades)) { // testa se o vertice destino ja esta presente no caminho

			res += adj[a][i]; // soma a distancia percorrida
			cidades.push_back(i); // insere o vertice destino no caminho

			bestC = bruteForceR(i, res, bestR, cidades, bestC); // chamada recursiva
			bestR = res; 

			if (cidades.size() == n && (bestR + adj[0][i]) < this->best) { // se o caminho é um circuito e é melhor que o melhor caso atual
				
				this->best = bestR + adj[0][i]; // atribui a distancia encontrada a melhor distancia 
				bestC = cidades; // atribui o caminho encontrado ao melhor caminho 

			}

			res -= adj[a][i]; // subtrai a distancia ate o vertice destino da distancia total
			cidades.pop_back(); // remove o vertice destino do caminho
		
		}//fim if

	}//fim for

	return bestC;

}//fim bruteForceR

/*
 Determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de branch and bound(chama a funcaoo recursiva)
 @return	resposta	Vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho
 @author				Henrique Schiess Pertussati
 */
vector<int> Graph::branchBound() {

	time_t inicio, fim;
	vector<int> cidades;
	vector<int> resposta;
	cidades.push_back(0);
	this->best = INFINITE;

	inicio = clock();
	resposta = branchBoundR(0, 0, 0, cidades, cidades);
	fim = clock();
	timeT = fim - inicio;

	return resposta;

}//fim branchBound

 /*
 Determina o menor caminho Hamiltoniano do grafo utilizando o paradigma de branch and bound
 @param		a		Vertice que estamos visitando
 @param		res		Caminho total atual
 @param		cidades	Vetor de inteiros com o caminho atual percorrido
 @param		bestC	Atual vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho
 @return	bestC	Atual vetor de inteiros que possui a ordem que os vertices devem ser visitados para obter o menor caminho
 @author			Arthur Ladislau Pereira
 @author			Henrique Schiess Pertussati
 */
vector<int> Graph::branchBoundR(int a, double res, double bestR, vector<int>cidades, vector<int> bestC) {

	for (int i = 0; i < n; i++) { // executa a permutacao as arestas do vertice corrente

		if (!isIn(i, cidades)) { // testa se o vertice destino ja esta presente no caminho
			
			res += adj[a][i]; // soma a distancia percorrida
			cidades.push_back(i); // insere o vertice destino no caminho

			if(res < this->best) { // verifica se o caminho em geracao ainda é válido

				bestC = branchBoundR(i, res, bestR, cidades, bestC);// chamada recursiva
				bestR = res;

				if (cidades.size() == n && (bestR + adj[0][i]) < this->best) { // se o caminho é um circuito e é melhor que o melhor caso atual
					this->best = bestR + adj[0][i]; // atribui a distancia encontrada a melhor distancia 
					bestC = cidades; // atribui o caminho encontrado ao melhor caminho 
				}
			}

			res -= adj[a][i]; // subtrai a distancia ate o vertice destino da distancia total
			cidades.pop_back(); // remove o vertice destino do caminho
		}
	}

	return bestC;
}

vector<int> Graph::dynamic() {
	long nsub = 1 << n;
	vector<vector<float>> opt(nsub, vector<float>(n));

	for (long s = 1; s < nsub; s += 2)
		for (int i = 1; i < n; ++i) {
			vector<int> subset;
			for (int u = 0; u < n; ++u)
				if (s & (1 << u))
					subset.push_back(u);

			if (subset.size() == 2)
				opt[s][i] = adj[0][i];

			else if (subset.size() > 2) {
				float min_subpath = INFINITE;
				long t = s & ~(1 << i);
				for (vector<int>::iterator j = subset.begin(); j != subset.end(); ++j)
					if (*j != i && opt[t][*j] + adj[*j][i] < min_subpath)
					min_subpath = opt[t][*j] + adj[*j][i];
				opt[s][i] = min_subpath;
			}
		}

	vector<int> tour;
	tour.push_back(0);

	bool selected[n];
	fill(selected, selected + n, false);
	selected[0] = true;

	long s = nsub - 1;

	for (int i = 0; i < n - 1; ++i) {
		int j = tour.back();
		float min_subpath = INFINITE;
		int best_k;
		for (int k = 0; k < n; ++k)
			if (!selected[k] && opt[s][k] + adj[k][j] < min_subpath) {
				min_subpath = opt[s][k] + adj[k][j];
				best_k = k;
			}
		tour.push_back(best_k);
		selected[best_k] = true;
		s -= 1 << best_k;
	}
	tour.push_back(0);

	return tour;
}

/*
Retorna se um elemento esta ou nao no vetor de inteiros
@param x	Elemento que deseja checar se existe
@param v	Vetor de inteiros usado para a checagem
@return	res	True(1) se o elemento esta no vetor, False(0) caso contrario
@author		Henrique Schiess Pertussati
*/
bool Graph::isIn(int x, vector<int> v) {

	bool res = false;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] == x) {
		
			res = true;
			i = v.size();
		}
	}
	return res;
}


/*
Funcao que imprime a matriz de adjacencia. Utilizada para debug
@author		Henrique Schiess Pertussati
*/
void Graph::printAdj() {
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < n; j++) 
			if (isConnected(i, j)) 
				cout << i << " " << j << " " << adj[i][j] << endl;
}

 /*
 Funcao que imprime a matriz de coordendas. Utilizada para debug
 @author		Henrique Schiess Pertussati
 */
void Graph::printCoor() {
	for (int i = 0; i < n; i++) {
		cout << i << " " << coor[i][0] << " " << coor[i][1] << endl;
	}
}