#include "Grafo.h"
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <stdlib.h>

using namespace std;

namespace Graph{

   /*
   *Construtor padrão para classe grafo
   */

   Grafo::Grafo(){
      this->num_vertices = 5;
      this->mat[5][5];
   }

   /*
   *Construtor padrão para classe grafo parametizado
   *@param num_vertices os vertices dos grafo
   *@param direc define se o grafo e ou nao direcionado
   */

   Grafo::Grafo(int num_vertices,bool direc){
      this->num_vertices = num_vertices;
      this->mat[num_vertices][num_vertices];
      this->direcionado = direc;

      for (int i = 0; i < num_vertices; i++){
         for (int j = 0; j < num_vertices; j++){
            this->mat[i][j] = 0.0;
         }
      }
   }

   /*
   *Metodo para inserir peso entre as arestas
   *@param v1 para primeiro vertice
   *@param v2 para segundo vertice
   *@param peso o valor de distancia entre os vertices
   */

   bool Grafo::inserirAresta(int v1,int v2,double peso){
      bool insert = true;
      if(!direcionado){
         this->mat[v1][v2] = peso;
         this->mat[v2][v1] = peso;
      }
      else{
         this->mat[v1][v2] = peso;
         this->mat[v2][v1] = -peso;
      }
      return insert;
   }

   /*
   *Metodo para remover aresta entre dois vertices
   *@v1 primeiro vertice
   *@v2 segundo vertice
   */

   bool Grafo::removerAresta(int v1, int v2){
      bool remove = true;
      this->mat[v1][v2];
      this->mat[v2][v1];
   }

   /*
   *Metodo para verificar se tem arestas entre os vertices
   *@v1 primeiro vertice
   *@v2 segundo vertice
   */

   bool Grafo::temAresta(int v1,int v2){
      bool resp = false;
      if(this->mat[v1][v2] > 0){
         resp = true;
      }

      return resp;
   }

   /*
   *Metodo para retornar o numero de arestas que contem o vertice
   */
   int Grafo::contArestas(){
      int cont = 0;
      for(int i = 0; i < this->num_vertices; i++){
         for(int j = 0; j < this->num_vertices; j++){
            if(this->mat[i][j] != 0 || this->mat[i][j] != -1){
               cont++;
            }
         }
      }
   }

   /*
   *Metodo para retonrar o numero de vertices
   */
   int Grafo::contVertices(){
      return this->num_vertices;
   }

   /*
   *Metodo para determinar se um grafo e completo
   */
   bool Grafo::eCompleto(){
      bool resp = false;
      int soma = 0;
      if(direcionado){
         return false;
      }
      else{
         for(int i = 0; i < this->num_vertices; i++){
            for(int j = 0; j < this->num_vertices; j++){
               if(this->mat[i][j]!=0){
                  soma +=1;
               }

            }
         }
      }

      if(soma == (this->num_vertices * this->num_vertices-1)/2){
         resp = true;
      }
      else{
         resp = false;
      }

      return resp;
   }

   /*
   *Metodo para retornar o numero de vertices adjacentes de um vertice
   *@param list para lista com os vertices
   *@param v o vertice
   */
   void Grafo::getAdjacentes(list<double>&a,int v){
      for(int i = 0; i < this->num_vertices; i++){
         if(this->mat[v][i] > 0){
            a.push_back(mat[v][i]);
         }
      }
   }

   /*
   *Metodo para retornar o antecessor do grafo
   *@param list para lista de vertices
   *@param v o vertices
   */

   void Grafo::getAntecessor(list<double>&a,int v){
      for(int i = 0; i < this->num_vertices; i++){
         if(this->mat[v][i] < 0){
            a.push_back(mat[v][i]);
         }
      }
   }

   /*
   *Metodo para retornar o complementar do grafo
   *@param Grafo que deseja o complementar
   */

   Grafo* Grafo::getComplementar(Grafo* g){
      Grafo* ret = new Grafo(g->num_vertices,false);

      for(int i = 0; i < this->num_vertices; i++){
         for(int j = 0; j < this->num_vertices; j++){
            if(temAresta(i,j) != true && i!=j){
               ret->inserirAresta(i,j,1);
            }
         }

      }
      return ret;
   }

   /*
   *Metodo para retornar o peso entre dos grafos
   *@param v1 vertice de origem
   *@param v2 vertice de destino
   */

   double Grafo::getPeso(int v1,int v2){
      return this->mat[v1][v2];
   }

   /*
   *Metodo para imprimir o grafo
   */

   void Grafo::imprimir(){
      cout<<"\t";
      for(int i = 0; i < this->num_vertices;i++)
         cout<<i<<"\t";
      cout<<"\n";
      for(int i = 0; i < this->num_vertices;i++){
         cout<<i<<"\t";
         for(int j = 0; j < this->num_vertices;j++)
            cout<<this->mat[i][j]<<"\t";
         cout<<"\n";
      }
   }


   /*
   *Construtor padrao para iniciar as condicoes para se rodar o algoritmo genetico
   *@param grafo o grafo representado
   *@param tamanhoPopulacao indica quantas rotas serao criadas inicialmente para determinar o menor custo
   *@param geracoes responsavel pelo numero de vezes que a populacao passara por mutacao
   *@param taxaMutacao o quanto a populacao ira sofrer o processo de mutacao
   *@param verticeInicial o vertice inicial do grafo
   */

   algoritmoGenetico::algoritmoGenetico(Grafo* grafo,int tamanhoPopulacao,int geracoes,int taxaMutacao,int verticeInicial){
      if(tamanhoPopulacao < 1){
         cout << "Erro, populacao pequena"<<endl;
         exit(1);
      }
      else if(taxaMutacao < 0 || taxaMutacao > 100){
         cout << " Taxa de mutacao invalida";
         exit(1);
      }else if(verticeInicial > grafo->num_vertices){
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

   double algoritmoGenetico::custoCaminho(vector<double>& solucao){
      double custoTotal = 0;
      set<double> solucoes;

      for(int i = 0; i < grafo->num_vertices;i++)
         solucoes.insert(solucao[i]);

      if(solucoes.size() != (unsigned) grafo->num_vertices) //Verifica se os elementos não se repetem
         return -1;

      for (int i = 0; i < grafo->num_vertices; i++){
         if(i+1 < grafo->num_vertices){
            double custo = grafo->getPeso(solucao[i],solucao[i+1]);
            //Determina se tem ligacoes entre os vertices
            if(custo < 0){
               return -1;
            }
            else{
               custoTotal+=custo;
            }
         }
         else{
            double custo = grafo->getPeso(solucao[i],solucao[0]);
            if(custo < 0){
               return -1;
            }
            else{
               custoTotal+=custo;
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

   bool algoritmoGenetico::existeCromossomos(const vector<double>& a){
      for(vector<pair<vector<double>, double> >::iterator it=populacao.begin(); it!=populacao.end(); it++){
         const vector<double>&b = (*it).first;
         if(equal(a.begin(),a.end(),b.begin()))
            return true;
      }

      return false;
   }

   /*
   *Metodo responsavel por gerar a populacao inicial
   */

   void algoritmoGenetico::iniciarPopulacao(){
      vector<double> pai;
      //Insere o vertice inicial no vetor pai
      pai.push_back(verticeInicial);
      //Cria o vetor pai
      for (int i = 0; i < grafo->num_vertices; i++){
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
         random_shuffle(pai.begin() + 1, pai.begin() + (rand() % (grafo->num_vertices - 1) + 1));
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

   /*
   *Metodo para imprimir a populacao(rota contendo os custos referente a elas)
   */

   void algoritmoGenetico::imprimirPopulacao(){
      for(vector<pair<vector<double>, double> >::iterator it=populacao.begin(); it!=populacao.end(); it++){
         const vector<double>& b = (*it).first;

         for(int i = 0; i < grafo->num_vertices;i++){
            cout << b[i] << " ";
         }
            cout << verticeInicial;
            cout <<" | Custo: " << (*it).second <<endl;

      }
      //cout << "\nTamanho da Populacao: " << tamanhoRealPopulacao << endl;
   }

   /*
   *Realiza a insercao dos filhos, atraves de insercao de pesquisa
   *@param filho gerado no algoritmo de cruzamento para ser inserido na populacao
   *@param custoTotal o custo da rota gerada
   */
   void algoritmoGenetico::insercaoArvorePesquisa(vector<double>& filho, double custoTotal){
      int min = 0;
      int max = tamanhoRealPopulacao - 1;

      while(max >= min)
      {
         int meio = min + (max - min) / 2;

         if(custoTotal == populacao[meio].second)
         {
            populacao.insert(populacao.begin() + meio, make_pair(filho, custoTotal));
            return;
         }
         else if(custoTotal > populacao[meio].second)
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

   void algoritmoGenetico::crossOver(vector<double>& pai1, vector<double>& pai2){

      vector<double> filho1, filho2;

      map<double, double> genes1, genes2;//Responsavel por o mapa de genes

      for(int i = 0; i < grafo->num_vertices; i++){
         genes1[pai1[i]] = 0;
         genes2[pai2[i]] = 0;
      }

      //Pontos aleatorios a serem selecionados na substring
      int ponto1 = rand() % (grafo->num_vertices - 1) + 1;
      int ponto2 = rand() % (grafo->num_vertices - ponto1) + ponto1;

      if(ponto1 == ponto2){
         if(ponto1 - 1 > 1)
            ponto1--;
         else if(ponto2 + 1 < grafo->num_vertices)
            ponto2++;
         else{
            int ponto = rand() % 10 + 1; // number in the range 1 to 10
            if(ponto <= 5)
               ponto1--;
            else
               ponto2++;
         }
      }

      //Gerar os "filhos"
      //Ate o ponto 1 ocorre a troca de genes entre os "filhos"

      for(int i = 0; i < ponto1; i++){
         //Insere os genes
         filho1.push_back(pai1[i]);
         filho2.push_back(pai2[i]);
         //Marca se foram usados
         genes1[pai1[i]] = 1;
         genes2[pai2[i]] = 1;
      }
      //Marca os genes restantes
      for(int i = ponto2 + 1; i < grafo->num_vertices; i++){
         genes1[pai1[i]] = 1;
         genes2[pai2[i]] = 1;
      }

      for(int i = ponto2; i >= ponto1; i--){
         if(genes1[pai2[i]] == 0){//Caso nao usado
            filho1.push_back(pai2[i]);
            genes1[pai2[i]] = 1;  //Marca o gene
         }
         else{
            for(map<double, double>::iterator it = genes1.begin(); it != genes1.end(); it++){
               if(it->second == 0){//checa se o gene nao foi usado
                  filho1.push_back(it->first);
                  genes1[it->first] = 1; //marca o gene
                  break;
               }
            }
         }

         if(genes2[pai1[i]] == 0){
            filho2.push_back(pai1[i]);
            genes2[pai1[i]] = 1;
         }
         else{
            for(map<double, double>::iterator it = genes2.begin(); it != genes2.end(); it++){
               if(it->second == 0){
                  filho2.push_back(it->first);
                  genes2[it->first] = 1;
                  break;
               }
            }
         }
      }



      for(int i = ponto2 + 1; i < grafo->num_vertices; i++){
         filho1.push_back(pai1[i]);
         filho2.push_back(pai2[i]);
      }

      //Efetua o processo de mutacao dos genes
      int mutacao = rand() % 100 + 1;
      if(mutacao <= taxaMutacao){//Checa se o numero gerado nao e menor que a taxa de mutacao
         //Efetua a mutacao
         int indexGene1, indexGene2;
         indexGene1 = rand() % (grafo->num_vertices - 1) + 1;
         indexGene2 = rand() % (grafo->num_vertices - 1) + 1;
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
      if(custoTotalFilho1 != -1 && !existeCromossomos(filho1)){
         insercaoArvorePesquisa(filho1, custoTotalFilho1);
         tamanhoRealPopulacao++;
      }

      if(custoTotalFilho2 != -1 && !existeCromossomos(filho2)){
         insercaoArvorePesquisa(filho2, custoTotalFilho2);
         tamanhoRealPopulacao++;
      }
   }

   /*
   *Metodo para gerar o menor custo
   */

   void algoritmoGenetico::gerarMenorCusto(){

      iniciarPopulacao();//Gera a populacao inicial

      if(tamanhoRealPopulacao == 0)
         return;

      for (int i = 0; i < geracoes; i++){
         int tempPopulacao = tamanhoRealPopulacao;//Tamanho real da populacao
         if(tamanhoRealPopulacao >= 2){//Seleciona dois pais para participar do processo de reproducao
            if(tamanhoRealPopulacao == 2){
               crossOver(populacao[0].first,populacao[1].first);//Aplica o algoritmo de cruzamento
            }
            else{
               int pai1,pai2;
               do{
                  pai1 = rand() % tamanhoRealPopulacao;//Seleciona pais aleatoriamente
                  pai2 = rand() % tamanhoRealPopulacao;
               }while(pai1 == pai2);
               crossOver(populacao[pai1].first,populacao[pai2].first);//Aplica a reproducao
            }

            int diffPopulacao = tamanhoRealPopulacao - tempPopulacao;//Checa se a populacao aumentou

            if(diffPopulacao == 2){
               if(tamanhoRealPopulacao > tamanhoPopulacao){
                  populacao.pop_back();//Remove os piores pais da populacao
                  populacao.pop_back();
                  tamanhoRealPopulacao -= 2;
               }
               else if(diffPopulacao == 1){
                  if(tamanhoRealPopulacao > tamanhoPopulacao){
                     populacao.pop_back();
                     tamanhoRealPopulacao--;
                  }
               }
            }
            else{
               crossOver(populacao[0].first, populacao[0].first);

               if(tamanhoPopulacao > tamanhoPopulacao){
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
      for(int i = 0; i < grafo->num_vertices; i++){
         cout << b[i]+1 << " ";
      }
      //cout << verticeInicial+1<<endl;
      cout<<endl;
   }

   /*
   *Retorna a menor rota apos os processos
   */

   double algoritmoGenetico::menorCusto(){
      if(tamanhoRealPopulacao > 0)
         return populacao[0].second;
      else
         return -1;
   }

}