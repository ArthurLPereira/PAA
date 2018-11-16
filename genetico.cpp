/************************************************
 * Projeto de Algoritmos - PUC Minas 15/11/2018
 * @author Matheus Kraisfeld
 * Solucoes para o Problema do Caixeiro Viajante.
 * Metodo Algortimo Genetico.
 ************************************************
 */
#include "Grafo.h"
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <time.h>

using namespace std;
using namespace Graph;

typedef
   struct coodernadas{
      int cidade;
      double x;
      double y;
   }coord;

/**
 * @brief Calcula a distancia entre duas cidades.
 * @param x1 coordenada x do primeiro vertice.
 * @param x2 coordenada x do segundo vertice.
 * @param y1 coordenada y do primeiro vertice.
 * @param y2 coordenada y do segundo vertice.
 * @return custo da distancia entre duas cidades.
 */
double distanciaDoisPontos(double x1, double x2, double y1, double y2){
   double resp = 0.0;
   double x,y,a_x,a_y;
   x = x2 - x1;
   y = y2 - y1;
   a_x = pow(x,2);
   a_y = pow(y,2);
   resp = sqrt(a_x+a_y);
   return resp;
}

/*Main*/
int main(){
   int num_vertices,tamanhoPopulacao,geracoes,taxaMutacao,verticeInicial;
   double x,y,peso,tempo;
   cin >> num_vertices;
   coord cidades[num_vertices];
   Grafo* g = new Grafo(num_vertices,false);

   //ler da entrada e criar um grafo com n vertices
   for(int i = 0; i < num_vertices; i++){
      cin>>x>>y;
      cidades[i].cidade = i;
      cidades[i].x = x;
      cidades[i].y = y;
   }

/* Calcular distancias entre as cidades. */
   for (int i = 0; i < num_vertices; i++){
      for (int j = 0; j < num_vertices; j++){
         peso = distanciaDoisPontos(cidades[i].x,cidades[j].x,cidades[i].y,cidades[j].y);
         g->inserirAresta(i,j,peso);
      }
   }

   //g->imprimir();

   algoritmoGenetico ag(g,40,100,20,0);

  /* Calcular o tempo de execucao. */
   //struct timespec tstart={0,0}, tend={0,0};
   //clock_gettime(CLOCK_MONOTONIC, &tstart);

 /* Obter solucao para o problema. */
    ag.gerarMenorCusto();

   //clock_gettime(CLOCK_MONOTONIC, &tend);
   //tempo = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) -((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
   //printf("\nTempo %.6f ms\n", tempo) ;
   return 0;
}