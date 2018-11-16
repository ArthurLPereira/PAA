#ifndef GRAFO_H_INCLUDED__
#define GRAFO_H_INCLUDED__

#include <iostream>
#include <list>
#include <vector>
#include <utility>
#include <map>
#include <set>
#include <stdlib.h>

using namespace std;

namespace Graph{
	class Grafo{
	
	private:	
		int num_vertices;
		bool direcionado;

	public:
		double mat[100][100];
		Grafo();
		Grafo(int,bool);
		bool inserirAresta(int,int,double);
		bool removerAresta(int ,int);
		bool temAresta(int,int);
		int contArestas();
		int contVertices();
		bool eCompleto();
		void getAdjacentes(list<double>&,int);
		void getAntecessor(list<double>&,int);
		Grafo* getComplementar(Grafo*);
		double getPeso(int,int);
		bool eComplementar(Grafo*);
		bool eDirecionado();
		void imprimir();

	   ~Grafo();
	   
	   friend class algoritmoGenetico;
	};

	typedef pair<vector<double>, double> popular;

struct ordenarPredecessor{
	bool operator()(const popular& primeiro, const popular& segundo)
		{
		return primeiro.second < segundo.second;
	}
};

class algoritmoGenetico
{

private:
	Grafo* grafo;
	std::vector<popular> populacao;
	int tamanhoPopulacao;
	int tamanhoRealPopulacao;
	int geracoes;
	int taxaMutacao;
	int verticeInicial;
	bool mostrarPopulacao;
private:
	void iniciarPopulacao();
public:
	algoritmoGenetico(Grafo*,int,int,int,int);
	double custoCaminho(vector<double>&);
	void imprimirPopulacao();
	void crossOver(vector<double>&,vector<double>&);
	void insercaoArvorePesquisa(vector<double>&,double);
	bool existeCromossomos(const vector<double>&);
	double menorCusto();
	void gerarMenorCusto();
	
	//~algoritmoGenetico();
	
};

}
#endif