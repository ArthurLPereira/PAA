#include <iostream>
#include <vector>
#include <ctime>
#include <list>
#include <utility>
#include <map>
#include <set>
#include <stdlib.h>
#include <algorithm>
#define INFINITE 999999 / 0.001
using namespace std;

/*
classe grafo para estruturar dados e
resolver o problema do caixeiro viajante
variaveis:
	n - numero de vertices
	adj - matriz de adjacencia. valor = -1 vertices nao conectados
	coor - coordenadas das cidades
@author		Henrique Schiess Pertussati
@author		Arthur Pereira
@author		Ricardo Sena
*/
class Graph{

private:

	int n;	//numero de vertices.
	vector<int> bruteForceR(int a, double res, double bestR,  vector<int>cidades, vector<int> bestC); //determina o menor caminho usando força bruta.
	vector<int> branchBoundR(int a, double res, double bestR, vector<int>cidades, vector<int> bestC); //determina o menor caminho usando branch and bound.
	vector<int> dynamicR(); //determina o menor caminho usando programação dinamica.
	bool isIn(int x, vector<int> v); //retorna se um elemento esta no vetor de inteiros.
	double best; //best vai ser usado pra calcular o PCV utilizando força bruta e branch and bound.

public:

	clock_t timeT;
	double **adj;	// matriz de adjacencia. se valor = -1, os vertices nao estao conectados.
	double **coor;	// coordenadas das cidades.
	Graph();		//construtor.
	Graph(int size);//construtor com paramentro size.
	~Graph();		//destrutor.
	int getSize();	//retorna valor tamanho n.
	void addCoor(int pos, double x, double y); //adiciona coordendas da cidade.
	void addEdge(int x, int y, double wt);	   //adiciona uma aresta ao grafo
	bool isConnected(int x, int y);			   //duas variaveis estao conectadas?
	vector<int> bruteForce();	//determina o menor caminho usando força bruta.
	vector<int> branchBound();  //determina o menor caminho usando branch and bound
	vector<int> dynamic();		//determina o menor caminho usando programação dinamica
	void printAdj();	//imprime a matriz de adjacencia
	void printCoor();	//imprime matriz de coordenadas
	double getPeso(int v1, int v2); //pegar peso entre vertices

	friend class algoritmoGenetico; //pode acessar private
};

typedef pair<vector<double>, double> popular;

struct ordenarPredecessor{
	bool operator()(const popular& primeiro, const popular& segundo){
		return primeiro.second < segundo.second;
	}
};

class algoritmoGenetico{
  private:
	Graph *grafo;
	std::vector<popular> populacao;
	int tamanhoPopulacao; //tamanho da populacao
	int tamanhoRealPopulacao;
	int geracoes; //quantidade de geracoes
	int taxaMutacao; //taxa usada pra mutacao
	int verticeInicial;
	bool mostrarPopulacao;

  private:
	void iniciarPopulacao();

  public:
	algoritmoGenetico(Graph *, int, int, int, int);
	double custoCaminho(vector<double> &);
	void imprimirPopulacao();
	void crossOver(vector<double> &, vector<double> &);
	void insercaoArvorePesquisa(vector<double> &, double);
	bool existeCromossomos(const vector<double> &);
	double menorCusto();
	void gerarMenorCusto();

	//~algoritmoGenetico();
};

/* 
Construtor
@param	size	Quantidade de vertices

@author			Ricardo Sena
@author			Desconhecido (checar Bibliografia[III])
@author			Arthur Pereira
*/
Graph::Graph(int size){

	if (size > 1){
		n = size;
	}
	else{
		n = 2;
	}

	adj = new double*[n];
	coor = new double*[n];

	for (int i = 0; i < n; i++){
		adj[i] = new double[n];	//matriz de adjacencia
		coor[i] = new double[n];//matriz de coordenadas

		for (int j = 0; j < n; j++){
			adj[i][j] = -1;	// Inicializa os vertices com -1
		}
		for (int j = 0; j < 2; j++){
			coor[i][j] = -1; //inicializa vertice com		
		}
	}
	this->best = INFINITE; //inicializa com "infinito"
	this->timeT = 0; 

}

/*
Construtor
*/
Graph::Graph(){
	n = 0;
	this->best = INFINITE;
}

/*
Destrutor
*/
Graph::~Graph() {
	for(int i = 0; i < this->n; i++){
		delete [] adj[i];
		delete [] coor[i];
	}
	delete [] adj;
	delete [] coor;
}


int Graph::getSize(){
	return n;
}


void Graph::addCoor(int pos, double x, double y){
	coor[pos][0] = x;
	coor[pos][1] = y;
}


void Graph::addEdge(int x, int y, double wt) {
	adj[x][y] = adj[y][x] = wt;
}


bool Graph::isConnected(int x, int y){
	bool res = false;

	if (adj[x][y] >= 0){
		res = true;
	}
	return res;
}

double Graph::getPeso(int v1, int v2) {
	return adj[v1][v2];
}

/*
algoritimo de força bruta que determina menor caminho
@return	resposta = ordem que vértices devem ser visitados para o menor caminho.
@author	Ricardo Sena
@author Henrique Schiess Pertussati
@author Arthur Pereira
*/
vector<int> Graph::bruteForce(){

	clock_t inicio, fim;

	vector<int> cidades; //vetor com as cidades
	vector<int> resposta; //vetor com a ordem das cidades da reposta
	cidades.push_back(0); 
  
	this->best = INFINITE;

	inicio = clock(); //começa contagem do tempo
	resposta = bruteForceR(0, 0, 0, cidades, cidades); //utiliza força bruta recursivamente 
	fim = clock(); //para o tempo
	timeT = fim - inicio; //checa tempo total gasto

	return resposta;
}

/*
algoritimo de força bruta que determina menor caminho
@param	a		vertice sendo visitado
@param	res		caminho atual total
@param	bestR	melhor Resposta com o caminho atual
@param	cidades	vetor com o caminho percorrido
@param	bestC	vetor com a ordem atual dos vértices

@return	bestC	retorna vetor dos inteiros com a ordem a ser seguida 	
@author			Ricardo Sena
@author     	Henrique Schiess Pertussati
@author			Arthur Pereira
*/
vector<int> Graph::bruteForceR(int a, double res, double bestR, vector<int>cidades, vector<int> bestC){
	for (int i = 0; i < n; i++){ //executa permutacao das arestas do vertice
		if(!isIn(i, cidades)){ //teste se vertice destino esta no caminho
			res += adj[a][i]; //soma distancia percorrida
			cidades.push_back(i); //insere vertice destino no caminho, operaçao relevante

			bestC = bruteForceR(i, res, bestR, cidades, bestC); //chamada recursiva
			bestR = res;

			if (cidades.size() == n && (bestR + adj[0][i]) < this->best){ //se caminho for circuito, e for melhor que o caminho atual
				this->best = bestR + adj[0][i]; //atribui como melhor distancia encontrada
				bestC = cidades; //atribui como melhor caminho enncontrado
			}
			res -= adj[a][i]; //subtrai distancia ate o vertice destino da distancia total
			cidades.pop_back(); //remove vertice destino do caminho percorrido
		}
	}
	return bestC;
}

/*
algoritimo de força bruta que determina menor caminho
@return	resposta = ordem que vértices devem ser visitados para o menor caminho.
@author      Henrique Schiess Pertussati
@author      Ricardo Sena
@author      Arthur Pereira
*/
vector<int> Graph::branchBound(){
	clock_t inicio, fim;
	vector<int> cidades;
	vector<int> resposta;
	cidades.push_back(0);
	this->best = INFINITE;

	inicio = clock(); //começa o tempo
	resposta = branchBoundR(0, 0, 0, cidades, cidades);
	fim = clock(); //para o tempo
	timeT = fim - inicio;

	return resposta; 
}

/*
algoritimo de branch and bround que determina menor caminho
 @param	 a		 vertice sendo visitado.
 @param	 res	 caminho atual total.
 @param	 bestR	 melhor Resposta com o caminho atual.
 @param	 cidades vetor com o caminho percorrido.
 @param	 bestC	 vetor com a ordem atual dos vértices.
 @return bestC	 retorna vetor dos inteiros com a ordem a ser seguida.
@author      Henrique Schiess Pertussati
@author      Ricardo Sena
@author      Arthur Pereira
 */
vector<int> Graph::branchBoundR(int a, double res, double bestR, vector<int>cidades, vector<int> bestC){
	for(int i = 0; i < n; i++){ //executa permutacao das arestas do vertice
        if (!isIn(i, cidades)){//teste se verice destino esta no caminho
            res += adj[a][i];  //soma distancia percorrida

            cidades.push_back(i); //insere vertice destino no caminho, operaçao relevante

            if(res < this->best){ //verifica se caminho gerado ainda é valido
				bestC = branchBoundR(i, res, bestR, cidades, bestC); //chamada recursiva
				bestR = res;

				if(cidades.size() == n && (bestR + adj[0][i]) < this->best){//se caminho for circuito, e for melhor que o caminho atual

                    this->best = bestR + adj[0][i]; //atribui como melhor distancia encontrada
                    bestC = cidades;                //atribui como melhor caminho enncontrado
                }
			}
            res -= adj[a][i];   //subtrai distancia ate o vertice destino da distancia total
            cidades.pop_back(); //remove vertice destino do caminho percorrido
        }
	}
	return bestC;
}

/*
algoritimo dinamico que determina menor caminho
@author      Ricardo Sena
@author      Arthur Pereira
 */
vector<int> Graph::dynamic(){
	vector<int> resposta;
	clock_t inicio, fim;
	inicio = clock();
	resposta = dynamicR();
	fim = clock();
	timeT = fim - inicio;

	return resposta;
}

/*
algoritimo dinamico que determina menor caminho
@author      Ricardo Sena
@author      Arthur Pereira
 */
vector<int> Graph::dynamicR(){
	long nsub = 1 << n;
	vector<vector<float>> opt(nsub, vector<float>(n));

	for (long s = 1; s < nsub; s += 2)
		for (int i = 1; i < n; ++i){
			vector<int> subset;
			for(int u = 0; u < n; ++u)
				if (s & (1 << u))
					subset.push_back(u);

			if (subset.size() == 2)
				opt[s][i] = adj[0][i];

			else if (subset.size() > 2){
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

	for (int i = 0; i < n - 1; ++i){
		int j = tour.back();
		float min_subpath = INFINITE;
		int best_k;
		for (int k = 0; k < n; ++k)
			if (!selected[k] && opt[s][k] + adj[k][j] < min_subpath){
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
retorna se elemento esta no vetor de inteiros
@param x	elemento que quer checar
@param v	vetor com os valores
@return		rep = true se estiver no vetor e false caso contrário.
@author      Henrique Schiess Pertussati
@author      Ricardo Sena
@author      Arthur Pereira
*/
bool Graph::isIn(int x, vector<int> v){
	bool res = false;

    //for para percorrer vetor
	for(int i = 0; i < v.size(); i++){
		if (v[i] == x){
			res = true;
			i = v.size(); //para sair do for
		}
	}
	return res;

}

/*
imprime a matriz de adjacencia.
@author      Henrique Schiess Pertussati
@author      Ricardo Sena
@author      Arthur Pereira
*/
void Graph::printAdj(){
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (isConnected(i, j)){
				cout << i << " " << j << " " << adj[i][j] << endl;
			}
		}
	}
}

/*
imprime a matriz de coordendas. 
@author      Henrique Schiess Pertussati
@author      Ricardo Sena
@author      Arthur Pereira
 */
void Graph::printCoor(){
	for (int i = 0; i < n; i++){
		cout << i << " " << coor[i][0] << " " << coor[i][1] << endl;
	}
}

algoritmoGenetico::algoritmoGenetico(Graph *grafo, int tamanhoPopulacao, int geracoes, int taxaMutacao, int verticeInicial)
{
	if (tamanhoPopulacao < 1)
	{
		cout << "Erro, populacao pequena" << endl;
		exit(1);
	}
	else if (taxaMutacao < 0 || taxaMutacao > 100)
	{
		cout << " Taxa de mutacao invalida";
		exit(1);
	}
	else if (verticeInicial > grafo->n)
	{
		cout << "Tamanho maior que numero de vertices";
		exit(1);
	}

	this->grafo = grafo;
	this->tamanhoPopulacao = tamanhoPopulacao;
	this->geracoes = geracoes;
	this->tamanhoRealPopulacao = 0;
	this->taxaMutacao = taxaMutacao;
	this->verticeInicial = verticeInicial;
}

/*
*Considera o custo de um caminho e sua solucao e valida
*@param solucao contendo a rota para se calcular o custo
*/

double algoritmoGenetico::custoCaminho(vector<double> &solucao)
{
	double custoTotal = 0;
	set<double> solucoes;

	for (int i = 0; i < grafo->n; i++)
		solucoes.insert(solucao[i]);

	if (solucoes.size() != (unsigned)grafo->n) //Verifica se os elementos não se repetem
		return -1;

	for (int i = 0; i < grafo->n; i++)
	{
		if (i + 1 < grafo->n)
		{
			double custo = grafo->getPeso(solucao[i], solucao[i + 1]);
			//Determina se tem ligacoes entre os vertices
			if (custo < 0)
			{
				return -1;
			}
			else
			{
				custoTotal += custo;
			}
		}
		else
		{
			double custo = grafo->getPeso(solucao[i], solucao[0]);
			if (custo < 0)
			{
				return -1;
			}
			else
			{
				custoTotal += custo;
				break;
			}
		}
	}

	return custoTotal;
}

/*
*Metodo responsavel por determinar se ja se encontra a rota dada na populacao
*@param a cromossomo contendo um rota
*/

bool algoritmoGenetico::existeCromossomos(const vector<double> &a)
{
	for (vector<pair<vector<double>, double>>::iterator it = populacao.begin(); it != populacao.end(); it++)
	{
		const vector<double> &b = (*it).first;
		if (equal(a.begin(), a.end(), b.begin()))
			return true;
	}

	return false;
}

/*
*Metodo responsavel por gerar a populacao inicial
*/

void algoritmoGenetico::iniciarPopulacao()
{
	vector<double> pai;
	//Insere o vertice inicial no vetor pai
	pai.push_back(verticeInicial);
	//Cria o vetor pai
	for (int i = 0; i < grafo->n; i++)
	{
		if (i != verticeInicial)
		{
			pai.push_back(i);
		}
	}

	double custoTotal = custoCaminho(pai);
	//Verifica se a rota contem ligacoes e uma rota valida, se for insere na populacao e incremeta o contador
	if (custoTotal != -1)
	{
		populacao.push_back(make_pair(pai, custoTotal));
		tamanhoRealPopulacao++;
	}

	//Cria rotas aleatorias, posteriormente verifica se tem custo maior que zero
	//e um crossomo ja inserido e faz a insercao na populacao
	for (int i = 0; i < geracoes; i++)
	{
		std::random_shuffle(pai.begin() + 1, pai.begin() + (rand() % (grafo->n - 1) + 1));
		double custoTotal = custoCaminho(pai);
		if (custoTotal != -1 && !existeCromossomos(pai))
		{
			populacao.push_back(make_pair(pai, custoTotal));
			tamanhoRealPopulacao++;
		}
		if (tamanhoPopulacao == tamanhoRealPopulacao)
			break;
	}

	if (tamanhoRealPopulacao == 0)
		cout << "\nPopulacao inicial vazia" << endl;
	else
		std::sort(populacao.begin(), populacao.end(), ordenarPredecessor());
}

/*
*Metodo para imprimir a populacao(rota contendo os custos referente a elas)
*/

void algoritmoGenetico::imprimirPopulacao()
{
	for (vector<pair<vector<double>, double>>::iterator it = populacao.begin(); it != populacao.end(); it++)
	{
		const vector<double> &b = (*it).first;

		for (int i = 0; i < grafo->n; i++)
		{
			cout << b[i] << " ";
		}
		cout << verticeInicial;
		cout << " | Custo: " << (*it).second << endl;
	}
	//cout << "\nTamanho da Populacao: " << tamanhoRealPopulacao << endl;
}

/*
*Realiza a insercao dos filhos, atraves de insercao de pesquisa
*@param filho gerado no algoritmo de cruzamento para ser inserido na populacao
*@param custoTotal o custo da rota gerada
*/
void algoritmoGenetico::insercaoArvorePesquisa(vector<double> &filho, double custoTotal)
{
	int min = 0;
	int max = tamanhoRealPopulacao - 1;

	while (max >= min)
	{
		int meio = min + (max - min) / 2;

		if (custoTotal == populacao[meio].second)
		{
			populacao.insert(populacao.begin() + meio, make_pair(filho, custoTotal));
			return;
		}
		else if (custoTotal > populacao[meio].second)
			min = meio + 1;
		else
			max = meio - 1;
	}
	populacao.insert(populacao.begin() + min, make_pair(filho, custoTotal));
}

/*
*Algoritmo de cruzamento que ira seleciona dois pontos aleatorio, que sera responsavel
*por gerar a filhos baseados nos pais
*A substring originaria dos pontos aleatorio do primeiro vetor sera inserida no segundo
*A substring do segundo sera inserido no primeiro
*Se torna um "filho" valido quando os genes ainda não foram usados
*/

void algoritmoGenetico::crossOver(vector<double> &pai1, vector<double> &pai2)
{

	vector<double> filho1, filho2;

	map<double, double> genes1, genes2; //Responsavel por o mapa de genes

	for (int i = 0; i < grafo->n; i++)
	{
		genes1[pai1[i]] = 0;
		genes2[pai2[i]] = 0;
	}

	//Pontos aleatorios a serem selecionados na substring
	int ponto1 = rand() % (grafo->n - 1) + 1;
	int ponto2 = rand() % (grafo->n - ponto1) + ponto1;

	if (ponto1 == ponto2)
	{
		if (ponto1 - 1 > 1)
			ponto1--;
		else if (ponto2 + 1 < grafo->n)
			ponto2++;
		else
		{
			int ponto = rand() % 10 + 1; // number in the range 1 to 10
			if (ponto <= 5)
				ponto1--;
			else
				ponto2++;
		}
	}

	//Gerar os "filhos"
	//Ate o ponto 1 ocorre a troca de genes entre os "filhos"

	for (int i = 0; i < ponto1; i++)
	{
		//Insere os genes
		filho1.push_back(pai1[i]);
		filho2.push_back(pai2[i]);
		//Marca se foram usados
		genes1[pai1[i]] = 1;
		genes2[pai2[i]] = 1;
	}
	//Marca os genes restantes
	for (int i = ponto2 + 1; i < grafo->n; i++)
	{
		genes1[pai1[i]] = 1;
		genes2[pai2[i]] = 1;
	}

	for (int i = ponto2; i >= ponto1; i--)
	{
		if (genes1[pai2[i]] == 0)
		{ //Caso nao usado
			filho1.push_back(pai2[i]);
			genes1[pai2[i]] = 1; //Marca o gene
		}
		else
		{
			for (map<double, double>::iterator it = genes1.begin(); it != genes1.end(); it++)
			{
				if (it->second == 0)
				{ //checa se o gene nao foi usado
					filho1.push_back(it->first);
					genes1[it->first] = 1; //marca o gene
					break;
				}
			}
		}

		if (genes2[pai1[i]] == 0)
		{
			filho2.push_back(pai1[i]);
			genes2[pai1[i]] = 1;
		}
		else
		{
			for (map<double, double>::iterator it = genes2.begin(); it != genes2.end(); it++)
			{
				if (it->second == 0)
				{
					filho2.push_back(it->first);
					genes2[it->first] = 1;
					break;
				}
			}
		}
	}

	for (int i = ponto2 + 1; i < grafo->n; i++)
	{
		filho1.push_back(pai1[i]);
		filho2.push_back(pai2[i]);
	}

	//Efetua o processo de mutacao dos genes
	int mutacao = rand() % 100 + 1;
	if (mutacao <= taxaMutacao)
	{ //Checa se o numero gerado nao e menor que a taxa de mutacao
		//Efetua a mutacao
		int indexGene1, indexGene2;
		indexGene1 = rand() % (grafo->n - 1) + 1;
		indexGene2 = rand() % (grafo->n - 1) + 1;
		//Insere os vertices nos "filhos"
		int aux = filho1[indexGene1];
		filho1[indexGene1] = filho1[indexGene2];
		filho1[indexGene2] = aux;

		aux = filho2[indexGene1];
		filho2[indexGene1] = filho2[indexGene2];
		filho2[indexGene2] = aux;
	}

	//Checa o custo dos filhos gerados
	double custoTotalFilho1 = custoCaminho(filho1);
	double custoTotalFilho2 = custoCaminho(filho2);
	//Caso nao tenha na populacao inserir
	if (custoTotalFilho1 != -1 && !existeCromossomos(filho1))
	{
		insercaoArvorePesquisa(filho1, custoTotalFilho1);
		tamanhoRealPopulacao++;
	}

	if (custoTotalFilho2 != -1 && !existeCromossomos(filho2))
	{
		insercaoArvorePesquisa(filho2, custoTotalFilho2);
		tamanhoRealPopulacao++;
	}
}

/*
*Metodo para gerar o menor custo
*/

void algoritmoGenetico::gerarMenorCusto()
{

	iniciarPopulacao(); //Gera a populacao inicial

	if (tamanhoRealPopulacao == 0)
		return;

	for (int i = 0; i < geracoes; i++)
	{
		int tempPopulacao = tamanhoRealPopulacao; //Tamanho real da populacao
		if (tamanhoRealPopulacao >= 2)
		{ //Seleciona dois pais para participar do processo de reproducao
			if (tamanhoRealPopulacao == 2)
			{
				crossOver(populacao[0].first, populacao[1].first); //Aplica o algoritmo de cruzamento
			}
			else
			{
				int pai1, pai2;
				do
				{
					pai1 = rand() % tamanhoRealPopulacao; //Seleciona pais aleatoriamente
					pai2 = rand() % tamanhoRealPopulacao;
				} while (pai1 == pai2);
				crossOver(populacao[pai1].first, populacao[pai2].first); //Aplica a reproducao
			}

			int diffPopulacao = tamanhoRealPopulacao - tempPopulacao; //Checa se a populacao aumentou

			if (diffPopulacao == 2)
			{
				if (tamanhoRealPopulacao > tamanhoPopulacao)
				{
					populacao.pop_back(); //Remove os piores pais da populacao
					populacao.pop_back();
					tamanhoRealPopulacao -= 2;
				}
				else if (diffPopulacao == 1)
				{
					if (tamanhoRealPopulacao > tamanhoPopulacao)
					{
						populacao.pop_back();
						tamanhoRealPopulacao--;
					}
				}
			}
			else
			{
				crossOver(populacao[0].first, populacao[0].first);

				if (tamanhoRealPopulacao > tamanhoPopulacao)
				{
					populacao.pop_back();
					tamanhoRealPopulacao--;
				}
			}
		}
	}

	// imprimirPopulacao();

	cout << "Custo: " << populacao[0].second<<endl;
	cout << "Caminho: ";
	const vector<double>& b = populacao[0].first;
	for(int i = 0; i < grafo->n; i++){
		cout << b[i]+1 << " ";
	}
	//cout << verticeInicial+1<<endl;
	cout<<endl;
}

/*
*Retorna a menor rota apos os processos
*/

double algoritmoGenetico::menorCusto()
{
	if (tamanhoRealPopulacao > 0)
		return populacao[0].second;
	else
		return -1;
}
