/* AnÃ¡lise nodal modificada de:
-  Ponto de OperaÃ§Ã£o
-  AnÃ¡lise no Estado Permanente
Elementos aceitos:
Resistor: R<nome> <nÃ³ +> <nÃ³ -> <ResistÃªncia>
Indutor: L<nome> <nÃ³ +> <nÃ³ -> <IndutÃ¢ncia>
Acoplamento entre indutores: K<nome> <La> <Lb> <k> (La e Lb nomes de indutores jÃ¡ declarados.)
Capacitor: C<nome> <nÃ³ +> <nÃ³ -> <CapacitÃ¢ncia>
Fonte de tensÃ£o controlada a tensÃ£o: E<nome> <nÃ³ V+> <nÃ³ V-> <nÃ³ v+> <nÃ³ v-> <Av>
Fonte de corrente controlada a corrente: F<nome> <nÃ³ I+> <nÃ³ I-> <nÃ³ i+> <nÃ³ i-> <Ai>
Fonte de corrente controlada a tensÃ£o: G<nome> <nÃ³ I+> <nÃ³ I-> <nÃ³ v+> <nÃ³ v-> <Gm>
Fonte de tensÃ£o controlada a corrente: H<nome> <nÃ³ V+> <nÃ³ V-> <nÃ³ i+> <nÃ³ i-> <Rm>
Fonte de corrente: I<nome> <nÃ³ +> <nÃ³ -> <mÃ³dulo> <fase (graus)> <valor contÃ­nuo>
Fonte de tensÃ£o: V<nome> <nÃ³ +> <nÃ³ -> <mÃ³dulo> <fase (graus)> <valor contÃ­nuo>
Amplificador operacional ideal: O<nome> <nÃ³ saÃ­da +> <nÃ³ saÃ­da -> <nÃ³ entrada +> <nÃ³ entrada ->
Transistor MOS: M<nome> <nÃ³ drain> <nÃ³ gate> <nÃ³ source> <nÃ³ base> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt 0> <lambda> <gamma> <theta> <Ld>
*/

/*
#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <complex.h>
#include <math.h>
#include <iostream>
*/

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <complex>
#include <iostream>
#include <string>
#include <time.h>
#include "fstream"

#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define LIMITE_LINHAS 5
#define TOLG_DC 1e-9
#define TOLG_AC 1e-30
#define DEBUG
#define FATORDC 1e9
#define MAX_ERRO 1e-9			//valor maximo de erro para iteracoes
#define MAX_ITER 50 			//maximo de iteracoes
#define REF_VAL 1 				//valor de referencia utilizado nos calculos de convergencias
#define MIN_ITER_CONV 2 		//numero minimo de iteracoes para considerar como estavel a solucao
#define TRY_CONV 5 				//o numero de tentativas que o algoritmo faz para calcular uma solucao, mesmo nao convergindo
#define J doubleComplex(0.0,1.0)
#define PI 3.141592653589793

#define REFVAL 1 //valor de referencia utilizado nos calculos de convergencias
#define MINCONV 2 //minimo de iteracoes para considerar como estavel a solucao
#define NCONV 5 // o numero de vezes que o algoritmo pode tentar calcular uma solucao
                //mesmo que esta nao esteja convergindo
#define MAX_IT 50


#define VALOR_UM 	 0.9999999999990
#define VALOR_ZERO 	 0.0000000000001

typedef std::complex<double> doubleComplex;

using namespace std;

enum pontoOperacao{
	corte, triodo, saturacao
};

enum tipoMOS{
	nmos,pmos};


typedef struct elemento { /* Defini??o de Elemento */
  char nome[MAX_NOME];
  double valor,modulo,fase;
  double cgB,cgS,cgD,gDS,gm,gmB;
  int a,b,c,d,x,y,tD,tG,tS,tB;  // Nos dos elementos

  double L,W,K,VT,LAMBDA,GAMMA,THETA,LD, ALPHA;
  char nomeLa[MAX_NOME], nomeLb[MAX_NOME], NPMOS[MAX_NOME];
  pontoOperacao operacaoTransistor;
  tipoMOS pnMOS;
} elemento;

elemento netlist[MAX_ELEM]; /* Lista de Elementos -> Netlist */

int
  ne, /* N?mero de Elementos */
  nv, /* N?mero de Variaveis */
  nn, /* N?mero de N?s */
  i,j,k;

char escala_frequencia[3] ;	// default para escolha de escala para frequencias
int freqInicialHz;		 // default de frequencia inicial para anÃ¡lise de pequenos sinais
int freqFinalHz;	// default de frequencia final para anÃ¡lise de pequenos sinais
int ptsFreq;			// default para quantidade de pontos obtidos entre limites de frequencia_variavel


char

/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  temp_char[MAX_LINHA+1],
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  ntD[MAX_NOME], ntG[MAX_NOME], ntS[MAX_NOME], ntB[MAX_NOME],
  ntTipo[MAX_NOME],nL[MAX_NOME],nW[MAX_NOME],nK[MAX_NOME],nVT[MAX_NOME],nLAMBDA[MAX_NOME],nGAMMA[MAX_NOME],nTHETA[MAX_NOME],nLD[MAX_NOME], //Vari?veis extras para o transistor
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;

FILE *arquivo;
FILE *tabelaFreq;


complex <double> gComplex;
complex <double> YnComplex[MAX_NOS+1][MAX_NOS+2];

double
  g,
  gm = 0,     //
  gDS = 0,   // // Variaveis para o transistor.
  gmB = 0,  //
  iO = 0,  //
  frequencia_variavel,
  tmp_err [MAX_NOS+1],
  tempVar [MAX_NOS+1],

  Yn[MAX_NOS+1][MAX_NOS+2],
  vAtual[MAX_NOS+1],
  vProximo[MAX_NOS+1];

  pontoOperacao operacaoTransistorAtual [MAX_NOS +1];
  pontoOperacao operacaoTransistorProximo [MAX_NOS +1];


  bool first_run = true;   // Primeira vez que roda a analise de pequenos sinais (escrever a primeira linha de informacoes da tabela)
  bool frequenciaHz = true; //Utiliza a frequencia_variavel em Hz por default
  bool linear = true;
  bool correct_model = true;
  bool flagDC = true; // Para comecar pela analise do ponto de operacao
  bool convergiu = false; // Comecar com false para entrar na primeira checagem


  int count_conv = 0; // Contador para convergencia
  int iteracoes = 0; // Numero de iteracoes
  int count_NOT_conv = 0; // Quantas vezes o algoritmo nao converge
  int vezes = 0;

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */

bool mantemModelo = true;
int contadorConv = 0; //conta quantas vezes os calculos deram erros menores
int vezNConvergiu = 0; //Conta a quantidade de vezes que o algoritmo nao convergiu
int numMaxIteracoes = 0;
int numRandIteracoes = 0;
double exc = 1.0;


 void calculoCapacitanciasParasitas(elemento netlist[]);
 void montarEstampa(double);
 int procuraIndutorTransformador(char *nomeElemento);



inline double sinDouble (double angulo)
 {
     double seno = sin( angulo * (PI/ 180.0));
     if (fabs(seno) < VALOR_ZERO)
         return (0.0);
     else if (fabs(seno) > VALOR_UM)
         return (1.0);
 		else
     	return (seno);
 }

inline double cosDouble (double angulo)
{
    double coseno = cos( angulo *  (PI/ 180.0) );
    if (fabs(coseno) < VALOR_ZERO)
        return (0.0);
    else if (fabs(coseno) > VALOR_UM)
        return (1.0);
    else
		return (coseno);
}


//funcao usada para criar a estampa AC do transformador
//funcionando
int procuraIndutorTransformador(char nomeElemento[])
{
	int i ;
	char nomeCompara[MAX_NOME+1];

  if (nomeElemento[0] != 'L') // verifica se elemento passado eh um indutor
  {
     printf("Um dos componentes do transformador nao eh um indutor\n");
     getch();
     exit(1);
  }
	else{
		 i = 1;
		 nomeCompara[0] = 'j';
     while(nomeElemento[i-1]!='\0')
     	{
          nomeCompara[i] = nomeElemento[i-1];
          i++;
     	}
		 nomeCompara[i] = '\0';

		 i = 0;
		 do{
			 if( !strcmp(nomeCompara,lista[i]) )
			 		return( i );
			 i++;
		 }while( i <= nv);

  		if (i==nv+1) {
      	printf("Nao encontrado o indutor especificado\n");
      	exit(1);
  		}
	}
}


int resolversistema_DC(void)
{
  int i,j,l, a;
  double t, p;

  for (i=1; i<=nv; i++) {
    t=0.0;
    a=i;
    for (l=i; l<=nv; l++) {
      if (fabs(Yn[l][i])>fabs(t)) {
	a=l;
	t=Yn[l][i];
      }
    }
    if (i!=a) {
      for (l=1; l<=nv+1; l++) {
	p=Yn[i][l];
	Yn[i][l]=Yn[a][l];
	Yn[a][l]=p;
      }
    }
    if (fabs(t)<TOLG_DC) {
      printf("Sistema singular\n");
      return 1;
    }
    for (j=nv+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
      Yn[i][j]/= t;
      p=Yn[i][j];
      if (p!=0)  /* Evita operacoes com zero */
        for (l=1; l<=nv; l++) {
	  if (l!=i)
	    Yn[l][j]-=Yn[l][i]*p;
        }
    }
  }
  return 0;
}

int resolversistema_AC(void)
{
  int i,j,l,a;
  complex <double> t,p;

  for (i = 1; i <= nv; i++) {
    t = 0.0 +  0.0 * J;
	p = 0.0 +  0.0 * J;
    a = i;
    for (l = i; l <= nv; l++)
      if (abs(YnComplex[l][i]) > abs(t)) {
		a = l;
		t = YnComplex[l][i];
      }
    if (i != a)
      for (l = 1; l <= nv+1; l++) {
		p = YnComplex[i][l];
		YnComplex[i][l] = YnComplex[a][l];
		YnComplex[a][l] = p;
      }
    if (abs(t) < TOLG_AC) {
      printf("Sistema singular\n");
      return 1;
    }
    for (j = nv+1; j > 0; j--) {  /* Basta j>i em vez de j>0 */
      YnComplex[i][j] /= t;
      p = YnComplex[i][j];
      if (abs(p) != 0.0)  /* Evita operacoes com zero */
        for (l = 1; l <= nv; l++)
			if (l != i)
				YnComplex[l][j] -= YnComplex[l][i] * p;
    }
  }
  return 0;
}

/* Rotina que conta os nos e atribui numeros a eles */
int numero(char *nome)
{
  int i,achou;

  i=0; achou=0;
  while (!achou && i<=nv)
    if (!(achou=!strcmp(nome,lista[i]))) i++;
  if (!achou) {
    if (nv==MAX_NOS) {
      printf("O programa so aceita ate %d nos\n",nv);
      exit(1);
    }
    nv++;
    strcpy(lista[nv],nome);
    return nv; /* novo no */
  }
  else {
    return i; /* no ja conhecido */
  }
    printf("Dentro da Convergencia: %d", convergiu );

	for (int asd = 0; asd <= nv; asd++)
		printf("%s\n", lista[i]);
}

void controleConvergencia ( double vAtual[], double vProximo[], int iteracoes )
{
  double maxVal = 0;

  printf ("\nCORRECT MODEL: %s\n", mantemModelo?"true":"false");
  printf ("CONVERGIU: %s\n", convergiu?"true":"false");
  getch();
  if (iteracoes < MAX_IT)
  {
     for (int cont = 0; cont <= nv; cont++)
     {
           //define o modo como o erro sera tratado
           if (fabs(vProximo[cont]) > REFVAL)
             tempVar[cont] = fabs( ( vProximo[cont] - vAtual[cont] ) / vProximo[cont]);

           if (vProximo[cont] < REFVAL)
             tempVar[cont] = fabs(vProximo[cont] - vAtual[cont]);

           //pega sempre o mesmo valor
           if (tempVar[cont] > maxVal)
             maxVal = tempVar[cont];
          }

          //faz com o que as tensoes futuras sejam as atuais
          for (int cont = 0; cont <= nv; cont++)
          {
             vAtual[cont] = vProximo[cont];
          }

          //printf("Erro: %f", MAX_ERRO);
         // printf("\nErro atual: %.10f\n", maxVal);
          if (maxVal <= MAX_ERRO)
          {
              mantemModelo = true;
              contadorConv++;
          }
          else if (maxVal >= MAX_ERRO)
          {
              vezNConvergiu++;
          }
          else if ((maxVal >= MAX_ERRO) && (vezNConvergiu >= NCONV))
          {
               mantemModelo = false;
               contadorConv = 0;
               //troca o modelo de transistor
               //e reinicia a contagem
               //se nao for pra manter o modelo
               //iteracoes = 0;
          }

          if (contadorConv >= MINCONV)
          {
               convergiu = true;
			   flagDC = false;
               printf("Convergiu \n");
               printf ("\nCORRECT MODEL: %s\n", mantemModelo?"true":"false");
  			   printf ("CONVERGIU: %s\n", convergiu?"true":"false");
  			  // getch();
          }
               //Mudanca 21/06
              //erroAtual = maxVal;

                //parte onde voce coloca o transistor
          if (mantemModelo)
          {
             iteracoes++;
//             cout << "Mantem Modelo:" << iteracoes <<endl;
//             getch();
          }

           else if (iteracoes >= MAX_IT)
           {
             //troca o modelo do transistor
             //e reinicia a contagem
             //se estourar o limite de iteracoes
             //iteracoes = 0;
             mantemModelo = false;
             contadorConv = 0;
           }
          numMaxIteracoes++;
          numRandIteracoes++;
          if (numMaxIteracoes >= 1000000) // VALOR QUE FUNCIONA 100000
          {
        	  //cout << "Como o Lukita diria: nem fodendo que esta merda converge =(" << endl;
        	  getch();
        	  exit (0);
          }
          if (numRandIteracoes >= 5000)
          {
			  for (int cont=0; cont<=nv;cont++)
			  {
				 if (tempVar[cont]>=MAX_ERRO)
				 {
					 vAtual [cont] = (rand () % (int) exc) - exc/2; //rand between -10 and 10
					 if (vAtual[cont] == 0) vAtual[cont] = 0.1;
					 vProximo [cont] = 0;
//					 cout << "trocando o valor cont:" << cont <<endl;
//					 getch();
				 }
			  }
			  numRandIteracoes = 0;
			  if (exc<10) exc+=0.5;
          }
     }
}

/*void controleConvergencia ( double vAtual[], double vProximo[], int iteracoes)
{

	double max_err = 0;

 	if (iteracoes < MAX_ITER){
		printf ("CORRECT MODEL: %s", correct_model?"true":"false");
		printf ("CONVERGIU: %s", convergiu?"true":"false");
 		for (int counter_var = 0; counter_var <= nv; counter_var++) {
 			//tmp_err = discrepância relativa
 			if (fabs(vProximo[counter_var]) > REF_VAL)
 				tmp_err[counter_var] = fabs( ( vProximo[counter_var] - vAtual[counter_var] ) / vProximo[counter_var]);

 			//tmp_err = discrepância absoluta
 			if (vProximo[counter_var] < REF_VAL)
 				tmp_err[counter_var] = fabs(vProximo[counter_var] - vAtual[counter_var]);

 			//atuliza max_err caso tmp_err seja maior
 			if (tmp_err[counter_var] > max_err)
 				max_err = tmp_err[counter_var];
         }

 		//faz com o que as tensoes futuras sejam as atuais
 		for (int counter_var = 0; counter_var <= nv; counter_var++)
 			vAtual[counter_var] = vProximo[counter_var];

 		//exibe o erro atual entre o nó atual e o nó futuro
 		printf("\nErro atual: %.10f\n", max_err);

 		// se convergir, mantem-se o modelo do transistor e incrementa-se o contador da quantidade de vezes que o algoritmo convergiu
 		if (max_err <= MAX_ERRO) {
 			correct_model = true;
 			count_conv++;
 		}

 		// senao, incrementa-se o contador da quantidade de vezes que o algoritmo NAO convergiu
 		else if (max_err >= MAX_ERRO)
 			count_NOT_conv++;

 		// caso "count_NOT_conv" supere TRY_CONV (quantidade maxima de tentativas para convergencia), troca-se o modelo do transistor e são zerados os contadores para iteracoes
 		if ((max_err >= MAX_ERRO) && (count_NOT_conv >= TRY_CONV)){
 			correct_model = false;
 			count_conv = 0;
 			iteracoes = 0;
 		}

		// caso "count_conv" supere MIN_ITER_CONV (quantidade minima de iteracoes de convergencias corretas), confirma-se que o modelo converge
 		if (count_conv >= MIN_ITER_CONV){
 			convergiu = true;
 			printf("Convergiu \n");
 		}

 		if (correct_model)
 			iteracoes++;
 	}
 	// se "iteracoes" superar MAX_ITER (quantidade maxima de iteracoes que a funcao controleConvergencia pode realizar), assume-se modelo incorreto e zera-se os contadores de iteracao
 	else if (iteracoes >= MAX_ITER) {
 		correct_model = false;
 		count_conv = 0;
 		iteracoes = 0;
 	}
*/

 void montarEstampa(){
	  /* Monta o sistema nodal modificado */
	  //printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
	  //getch();
	  /* Zera sistema */
	  for (i=0; i<=nv; i++) {
		for (j=0; j<=nv+1; j++){
		  Yn[i][j]=0;                              //Inicializacao dos coeficientes
		  YnComplex[i][j] = 0.0 + 0.0 * J;

		}
	  }
	  /* Monta estampas */
	  int numMOS = 0;
	  for (i=1; i<=ne; i++) {
		tipo=netlist[i].nome[0];
		if (tipo=='R') {						//Monta estampa para o RESISTOR, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
			  g=1/netlist[i].valor;
			  Yn[netlist[i].a][netlist[i].a]+=g;
			  Yn[netlist[i].b][netlist[i].b]+=g;
			  Yn[netlist[i].a][netlist[i].b]-=g;
			  Yn[netlist[i].b][netlist[i].a]-=g;
			}
			else{
				 g = 1.0/(netlist[i].valor);
				 YnComplex[netlist[i].a][netlist[i].a]+=g;
				 YnComplex[netlist[i].b][netlist[i].b]+=g;
				 YnComplex[netlist[i].a][netlist[i].b]-=g;
				 YnComplex[netlist[i].b][netlist[i].a]-=g;
			}
		}
		else if (tipo=='K') {
			if(flagDC){
				continue;
			}
			else{
				int indutorLa = procuraIndutorTransformador(netlist[i].nomeLa);
				int indutorLb = procuraIndutorTransformador(netlist[i].nomeLb);

				double valLa = netlist[ indutorLa ].valor;
                double valLb = netlist[ indutorLb ].valor;
                printf("valLa: %.6f valLb: %.6f\n", valLa, valLb);

                double M = netlist[i].valor * sqrt(valLa * valLb);
                YnComplex[indutorLa][indutorLb] += 0.0 + J * 2.0 * PI * frequencia_variavel * M;
                YnComplex[indutorLb][indutorLa] += 0.0 + J * 2.0 * PI * frequencia_variavel * M;
			}
		}
		else if (tipo=='C'){							//Monta estampa para o CAPACITOR, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
				g=netlist[i].valor / FATORDC;
				Yn[netlist[i].a][netlist[i].a]+=g;
				Yn[netlist[i].b][netlist[i].b]+=g;
				Yn[netlist[i].a][netlist[i].b]-=g;
				Yn[netlist[i].b][netlist[i].a]-=g;
			}
			else{
				if(frequenciaHz)
					gComplex = 0.0 + J * 2.0 * PI * netlist[i].valor * frequencia_variavel;  //Para frequencia_variavel em Hz
				else
					gComplex = 0.0 + J * netlist[i].valor * frequencia_variavel;  		//Para frequencia_variavel em Rad/s
				YnComplex[netlist[i].a][netlist[i].a]+=gComplex;
				YnComplex[netlist[i].b][netlist[i].b]+=gComplex;
				YnComplex[netlist[i].a][netlist[i].b]-=gComplex;
				YnComplex[netlist[i].b][netlist[i].a]-=gComplex;
			}
		}
		else if (tipo=='L') {
			if(flagDC){
				g=1/FATORDC;											//FATORDC = 10e8 portando para DC, condutancia muito pequeno (curto)
				Yn[netlist[i].a][netlist[i].x]+=1;
				Yn[netlist[i].b][netlist[i].x]-=1;
				Yn[netlist[i].x][netlist[i].a]-=1;
				Yn[netlist[i].x][netlist[i].b]+=1;
				Yn[netlist[i].x][netlist[i].x]+= g;
			}else{
				if(frequenciaHz){
					gComplex = 0.0  +
					 J * 2.0 * PI * netlist[i].valor * frequencia_variavel;			//Para frequencia_variavel em Hz	    (REVISAR)
					
				}
				else
					gComplex = 0.0 + J * frequencia_variavel * netlist[i].valor;					//Para frequencia_variavel em Rad/s
				YnComplex[netlist[i].a][netlist[i].x]+=1.0 + 0.0*J;
				YnComplex[netlist[i].b][netlist[i].x]-=1.0 + 0.0*J;
				YnComplex[netlist[i].x][netlist[i].a]-=1.0 + 0.0*J;
				YnComplex[netlist[i].x][netlist[i].b]+=1.0 + 0.0*J;
				YnComplex[netlist[i].x][netlist[i].x]+= gComplex;
			}
		}

		else if (tipo=='G') {										//Monta estampa para a FONTE DE CORRENTE CONTROLADA A TENSAO, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
				g=netlist[i].valor;
				Yn[netlist[i].a][netlist[i].c]+=g;
				Yn[netlist[i].b][netlist[i].d]+=g;
				Yn[netlist[i].a][netlist[i].d]-=g;
				Yn[netlist[i].b][netlist[i].c]-=g;
			}else{
				g=netlist[i].valor;
				YnComplex[netlist[i].a][netlist[i].c]+=g;
				YnComplex[netlist[i].b][netlist[i].d]+=g;
				YnComplex[netlist[i].a][netlist[i].d]-=g;
				YnComplex[netlist[i].b][netlist[i].c]-=g;
			}
		}
		else if (tipo=='I') {										//Monta estampa para a FONTE DE CORRENTE, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
				g=netlist[i].valor;
				Yn[netlist[i].a][nv+1]-=g;
				Yn[netlist[i].b][nv+1]+=g;
			}else{
				gComplex = netlist[i].modulo * cosDouble(netlist[i].fase) + J*netlist[i].modulo * sinDouble(netlist[i].fase);
				YnComplex[netlist[i].a][nv+1]-=gComplex;
				YnComplex[netlist[i].b][nv+1]+=gComplex;
			}
		}
		else if (tipo=='V') {										//Monta estampa para a FONTE DE TENSAO, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
				Yn[netlist[i].a][netlist[i].x]+=1;
				Yn[netlist[i].b][netlist[i].x]-=1;
				Yn[netlist[i].x][netlist[i].a]-=1;
				Yn[netlist[i].x][netlist[i].b]+=1;
				Yn[netlist[i].x][nv+1]-=netlist[i].valor;
			}else {
				gComplex = netlist[i].modulo * cosDouble(netlist[i].fase) + J*netlist[i].modulo * sinDouble(netlist[i].fase);
				YnComplex[netlist[i].a][netlist[i].x]+=1;
				YnComplex[netlist[i].b][netlist[i].x]-=1;
				YnComplex[netlist[i].x][netlist[i].a]-=1;
				YnComplex[netlist[i].x][netlist[i].b]+=1;
				YnComplex[netlist[i].x][nv+1]-=gComplex;
			}
		}
		else if (tipo=='E') {										//Monta estampa para a FONTE DE TENSAO CONTROLADA A TENSAO, dependendo do flagDC monta estampa DC ou AC.

		  if(flagDC){
				g=netlist[i].valor;
				Yn[netlist[i].a][netlist[i].x]+=1;
				Yn[netlist[i].b][netlist[i].x]-=1;
				Yn[netlist[i].x][netlist[i].a]-=1;
				Yn[netlist[i].x][netlist[i].b]+=1;
				Yn[netlist[i].x][netlist[i].c]+=g;
				Yn[netlist[i].x][netlist[i].d]-=g;
			} else {
				gComplex=netlist[i].valor;
				YnComplex[netlist[i].a][netlist[i].x]+=1;
				YnComplex[netlist[i].b][netlist[i].x]-=1;
				YnComplex[netlist[i].x][netlist[i].a]-=1;
				YnComplex[netlist[i].x][netlist[i].b]+=1;
				YnComplex[netlist[i].x][netlist[i].c]+=gComplex;
				YnComplex[netlist[i].x][netlist[i].d]-=gComplex;
			}
		}
		else if (tipo=='F') {									//Monta estampa para a FONTE DE CORRENTE CONTROLADA A CORRENTE, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
				g=netlist[i].valor;
				Yn[netlist[i].a][netlist[i].x]+=g;
				Yn[netlist[i].b][netlist[i].x]-=g;
				Yn[netlist[i].c][netlist[i].x]+=1;
				Yn[netlist[i].d][netlist[i].x]-=1;
				Yn[netlist[i].x][netlist[i].c]-=1;
				Yn[netlist[i].x][netlist[i].d]+=1;
			}else {
				gComplex=netlist[i].valor;
				YnComplex[netlist[i].a][netlist[i].x]+=gComplex;
				YnComplex[netlist[i].b][netlist[i].x]-=gComplex;
				YnComplex[netlist[i].c][netlist[i].x]+=1;
				YnComplex[netlist[i].d][netlist[i].x]-=1;
				YnComplex[netlist[i].x][netlist[i].c]-=1;
				YnComplex[netlist[i].x][netlist[i].d]+=1;
			}
		}
		else if (tipo=='H') {									//Monta estampa para a FONTE DE TENSAO CONTROLADA A CORRENTE, dependendo do flagDC monta estampa DC ou AC.
		  if(flagDC){
				g=netlist[i].valor;
				Yn[netlist[i].a][netlist[i].y]+=1;
				Yn[netlist[i].b][netlist[i].y]-=1;
				Yn[netlist[i].c][netlist[i].x]+=1;
				Yn[netlist[i].d][netlist[i].x]-=1;
				Yn[netlist[i].y][netlist[i].a]-=1;
				Yn[netlist[i].y][netlist[i].b]+=1;
				Yn[netlist[i].x][netlist[i].c]-=1;
				Yn[netlist[i].x][netlist[i].d]+=1;
				Yn[netlist[i].y][netlist[i].x]+=g;
			}else{
				gComplex=netlist[i].valor;
				YnComplex[netlist[i].a][netlist[i].y]+=1;
				YnComplex[netlist[i].b][netlist[i].y]-=1;
				YnComplex[netlist[i].c][netlist[i].x]+=1;
				YnComplex[netlist[i].d][netlist[i].x]-=1;
				YnComplex[netlist[i].y][netlist[i].a]-=1;
				YnComplex[netlist[i].y][netlist[i].b]+=1;
				YnComplex[netlist[i].x][netlist[i].c]-=1;
				YnComplex[netlist[i].x][netlist[i].d]+=1;
				YnComplex[netlist[i].y][netlist[i].x]+=gComplex;
			}
		}
		 else if (tipo=='M'){                                    //Monta estampa para o TRANSISTOR PMOS OU NMOS, Por enquanto apenas DC.
            if (flagDC) {
                int numMOS = 0;
                g = ((double) 1.0) / FATORDC;
                netlist[i].gm = 0;
                netlist[i].gDS = 0;
                iO = 0;
                printf("vd %f vs %f\n", vAtual[netlist[i].tD],vAtual[netlist[i].tS]);

                if((vAtual[netlist[i].tD] < vAtual[netlist[i].tS] && netlist[i].pnMOS == nmos) || (vAtual[netlist[i].tD] < vAtual[netlist[i].tS] && netlist[i].pnMOS == nmos)) {
                    printf("Tensao no Drain > Tensao no Source\n");
                    int aux = netlist[i].tD;
                    netlist[i].tD = netlist[i].tS;
                    netlist[i].tS = aux;
                }

                double vDS = vAtual[netlist[i].tD] - vAtual[netlist[i].tS];
                double vGS = vAtual[netlist[i].tG] - vAtual[netlist[i].tS];
                double vBS = vAtual[netlist[i].tB] - vAtual[netlist[i].tS];
                if(netlist[i].pnMOS == pmos){
                    vDS *= -1;
                    vBS *= -1;
                    vGS *= -1;
                }


                double vt = netlist[i].VT + netlist[i].GAMMA * (sqrt(fabs(netlist[i].THETA - (vBS))) - sqrt(fabs(netlist[i].THETA)));
                printf ("vt: %f  vGS: %f vDS: %f vS: %f tS %d\n", vt,vGS,vDS, vAtual[netlist[i].tS], netlist[i].tS);
                if (vGS < vt){
                    netlist[i].operacaoTransistor = corte;
                    printf ("Modo de operacao: corte\n");
                }
                else if (vDS < vGS - vt){
                    netlist[i].operacaoTransistor = triodo;
                    printf("Modo de operacao: triodo\n");
                    netlist[i].gm = netlist[i].K * (netlist[i].W/netlist[i].L)*(2* (vDS))*(1+ netlist[i].LAMBDA* vDS);
                    netlist[i].gDS = netlist[i].K * (netlist[i].W/netlist[i].L) * (2*(vGS - vt) - 2 * vDS + 4* netlist[i].LAMBDA * ( vGS - vt) * (vDS) - 3*netlist[i].LAMBDA * pow ( vDS,2));
                    iO = netlist[i].K * (netlist[i].W/netlist[i].L) * (2* (vGS - vt)*vDS - pow (vDS,2)) - (netlist[i].gm * vGS) - (gDS * vDS);
                }
                else {
                    netlist[i].operacaoTransistor = saturacao;
                    printf("Modo de operacao: saturacao\n");
                    netlist[i].gm = netlist[i].K * (netlist[i].W/netlist[i].L)* 2 *(vGS - vt) * (1 + netlist[i].LAMBDA* vDS);
                    netlist[i].gDS = netlist[i].K * (netlist[i].W/netlist[i].L)* pow ((vGS - vt),2) * netlist[i].LAMBDA;
                    iO = netlist[i].K * (netlist[i].W/netlist[i].L) * pow((vGS - vt),2) * (1 + netlist[i].LAMBDA * vDS)- (netlist[i].gm * vGS) - (netlist[i].gDS * vDS);
                }

                netlist[i].gmB = (netlist[i].gm*netlist[i].GAMMA)/(2*sqrt(fabs(netlist[i].THETA - (vBS))));
                iO-= (netlist[i].gmB * (vBS));
                iO *= (netlist[i].pnMOS == pmos?-1:1);
                printf ("gm %f gmb %f  gDS %f io %f \n", netlist[i].gm, netlist[i].gmB, netlist[i].gDS, iO);

                Yn[netlist[i].tD][netlist[i].tB]+=netlist[i].gmB;
                Yn[netlist[i].tS][netlist[i].tS]+=netlist[i].gmB;
                Yn[netlist[i].tD][netlist[i].tS]-=netlist[i].gmB;
                Yn[netlist[i].tS][netlist[i].tB]-=netlist[i].gmB;

                Yn[netlist[i].tD][netlist[i].tG]+=netlist[i].gm;
                Yn[netlist[i].tS][netlist[i].tS]+=netlist[i].gm;
                Yn[netlist[i].tD][netlist[i].tS]-=netlist[i].gm;
                Yn[netlist[i].tS][netlist[i].tG]-=netlist[i].gm;

                Yn[netlist[i].tD][netlist[i].tD]+=netlist[i].gDS;
                Yn[netlist[i].tS][netlist[i].tS]+=netlist[i].gDS;
                Yn[netlist[i].tD][netlist[i].tS]-=netlist[i].gDS;
                Yn[netlist[i].tS][netlist[i].tD]-=netlist[i].gDS;


                Yn[netlist[i].tD][netlist[i].tD]+=g;
                Yn[netlist[i].tG][netlist[i].tG]+=g;
                Yn[netlist[i].tD][netlist[i].tG]-=g;
                Yn[netlist[i].tG][netlist[i].tD]-=g;

                Yn[netlist[i].tS][netlist[i].tS]+=g;
                Yn[netlist[i].tG][netlist[i].tG]+=g;
                Yn[netlist[i].tS][netlist[i].tG]-=g;
                Yn[netlist[i].tG][netlist[i].tS]-=g;

                Yn[netlist[i].tB][netlist[i].tB]+=g;
                Yn[netlist[i].tG][netlist[i].tG]+=g;
                Yn[netlist[i].tB][netlist[i].tG]-=g;
                Yn[netlist[i].tG][netlist[i].tB]-=g;

                Yn[netlist[i].tD][nv+1]-=iO;
                Yn[netlist[i].tS][nv+1]+=iO;
            }
			else {
				complex <double> gCGS = 0.0 + J * 0.0;
				complex <double> gCGD = 0.0 + J * 0.0;
				complex <double> gCGB = 0.0 + J * 0.0;

				if (frequenciaHz) {
					gCGB = 0.0 + J * 2.0 * PI * netlist[i].cgB * frequencia_variavel;
					gCGS = 0.0 + J * 2.0 * PI * netlist[i].cgS * frequencia_variavel;
					gCGD = 0.0 + J * 2.0 * PI * netlist[i].cgD * frequencia_variavel;

				}
				else {
					gCGB = 0.0 + J * netlist[i].cgB * frequencia_variavel;
					gCGS = 0.0 + J * netlist[i].cgS * frequencia_variavel;
					gCGD = 0.0 + J * netlist[i].cgD * frequencia_variavel;
				}

				YnComplex[netlist[i].tD][netlist[i].tB]+=netlist[i].gmB;
				YnComplex[netlist[i].tS][netlist[i].tS]+=netlist[i].gmB;
				YnComplex[netlist[i].tD][netlist[i].tS]-=netlist[i].gmB;
				YnComplex[netlist[i].tS][netlist[i].tB]-=netlist[i].gmB;

				YnComplex[netlist[i].tD][netlist[i].tG]+=netlist[i].gm;
				YnComplex[netlist[i].tS][netlist[i].tS]+=netlist[i].gm;
				YnComplex[netlist[i].tD][netlist[i].tS]-=netlist[i].gm;
				YnComplex[netlist[i].tS][netlist[i].tG]-=netlist[i].gm;

				YnComplex[netlist[i].tD][netlist[i].tD]+=netlist[i].gDS;
				YnComplex[netlist[i].tS][netlist[i].tS]+=netlist[i].gDS;
				YnComplex[netlist[i].tD][netlist[i].tS]-=netlist[i].gDS;
				YnComplex[netlist[i].tS][netlist[i].tD]-=netlist[i].gDS;

				YnComplex[netlist[i].tD][netlist[i].tD]+=gCGD;
				YnComplex[netlist[i].tG][netlist[i].tG]+=gCGD;
				YnComplex[netlist[i].tD][netlist[i].tG]-=gCGD;
				YnComplex[netlist[i].tG][netlist[i].tD]-=gCGD;

				YnComplex[netlist[i].tS][netlist[i].tS]+=gCGS;
				YnComplex[netlist[i].tG][netlist[i].tG]+=gCGS;
				YnComplex[netlist[i].tS][netlist[i].tG]-=gCGS;
				YnComplex[netlist[i].tG][netlist[i].tS]-=gCGS;

				YnComplex[netlist[i].tB][netlist[i].tB]+=gCGB;
				YnComplex[netlist[i].tG][netlist[i].tG]+=gCGB;
				YnComplex[netlist[i].tB][netlist[i].tG]-=gCGB;
				YnComplex[netlist[i].tG][netlist[i].tB]-=gCGB;
			}
		}
		else if (tipo=='O') {
			if (flagDC) {
				Yn[netlist[i].a][netlist[i].x]+=1;
				Yn[netlist[i].b][netlist[i].x]-=1;
				Yn[netlist[i].x][netlist[i].c]+=1;
				Yn[netlist[i].x][netlist[i].d]-=1;
			}
			else {
				YnComplex[netlist[i].a][netlist[i].x]+=1;
				YnComplex[netlist[i].b][netlist[i].x]-=1;
				YnComplex[netlist[i].x][netlist[i].c]+=1;
				YnComplex[netlist[i].x][netlist[i].d]-=1;
			}
		}
	}
//	getch();
}

void calculoCapacitanciasParasitas(elemento *netlist)
{
	if(netlist->operacaoTransistor == corte){
		netlist->cgB = netlist->ALPHA * netlist->W * netlist->L;
		netlist->cgS = netlist->ALPHA * netlist->W * netlist->LD;
		netlist->cgD = netlist->ALPHA * netlist->W * netlist->LD;
	}
	else if (netlist->operacaoTransistor == triodo){
		netlist->cgB = 0;
		netlist->cgS = (1.0/2 * netlist->ALPHA * netlist->W * netlist->L) + (netlist->ALPHA * netlist->W * netlist->LD);
		netlist->cgD = (1.0/2 * netlist->ALPHA * netlist->W * netlist->L) + (netlist->ALPHA * netlist->W * netlist->LD);
	}
	else{
		netlist->cgB = 0;
		netlist->cgS = (2.0/3 * netlist->ALPHA * netlist->W * netlist->L) + (netlist->ALPHA * netlist->W * netlist->LD);
		netlist->cgD = netlist->ALPHA * netlist->W * netlist->LD;
	}
 }

void escreverTabelaFreq ()
{
   

    // substitui em uma string auxiliar o nome do arquivo de ".net" para ".tab"
    if (first_run) {
        strcpy (temp_char,nomearquivo);
        *strstr(temp_char,".net") = '\0';
        strcat (temp_char,".tab");
    }

    // ao iniciar a analise de pequenos sinais, cria-se o arquivo de tabela e coloca-se na 1a linha os nomes das colunas da tabela
    if (first_run) {
        tabelaFreq = fopen (temp_char,"w");
        fprintf(tabelaFreq,"f ");
        for(int i = 1; i < nv + 1; i++)
            fprintf(tabelaFreq, "%sm %sf ", lista[i], lista[i]);
        fprintf(tabelaFreq, "\n");
        fprintf (tabelaFreq, "%g ", frequencia_variavel);
        for(int i = 1; i < nv + 1; i++) {
            fprintf(tabelaFreq, "%g ", abs(YnComplex[i][nv+1]));
            fprintf(tabelaFreq, "%g ", ((180.0 / PI) * arg(YnComplex[i][nv+1])));
        }
        fprintf(tabelaFreq, "\n");
        fclose(tabelaFreq);
        first_run = false;
    }
    // escreve a frequencia atual e os valores de modulo/fase para as variaveis do sistema
    else {
    	tabelaFreq = fopen(temp_char,"a");
        fprintf (tabelaFreq, "%g ", frequencia_variavel);
        for(int i = 1; i < nv + 1; i++) {
            fprintf(tabelaFreq, "%g ", abs(YnComplex[i][nv+1]));
            fprintf(tabelaFreq, "%g ", ((180.0 / PI) * arg(YnComplex[i][nv+1])));
        }
        fprintf(tabelaFreq, "\n");
        fclose(tabelaFreq);
    }
}

void lerNetlist(){

  printf("Lendo netlist:\n");
  fgets(txt,MAX_LINHA,arquivo);
  printf("Titulo: %s",txt);
  while (fgets(txt,MAX_LINHA,arquivo)) {
    ne++; /* Nao usa o netlist[0] */
    if (ne>MAX_ELEM) {
      printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
      exit(1);
    }
    txt[0]=toupper(txt[0]);
    tipo=txt[0];
    sscanf(txt,"%10s",netlist[ne].nome);
    p=txt+strlen(netlist[ne].nome); /* Inicio dos parametros */
    /* O que e lido depende do tipo */
		if(tipo=='.'){
				ne--;	
				sscanf(p,"%10s%d%d%d",&escala_frequencia, &ptsFreq, &freqInicialHz, &freqFinalHz);
				
				printf("Analise AC\n");
				printf("Modo : %s\n", escala_frequencia);
				printf("Pontos : %d\n", ptsFreq);
				printf("Frequencia Inicial: %d\n", freqInicialHz);
				printf("Frequencia Final: %d\n", freqFinalHz);
				printf("\n");
				
		}
    else if (tipo=='R' || tipo=='L' || tipo=='C') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo == 'V' || tipo=='I'){
    sscanf(p,"%10s%10s%lg%lg%lg",na,nb,&netlist[ne].modulo,&netlist[ne].fase,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);

	}
    else if (tipo=='K') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
    }
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[ne].valor);
      printf("%s %s %s %s %s %g\n",netlist[ne].nome,na,nb,nc,nd,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
    else if (tipo=='O') {
      sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
      printf("%s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
    }
	else if (tipo=='M') {
	  sscanf (p, "%10s%10s%10s%10s%10s L=%lf W=%lf %lf%lf%lf%lf%lf%lf", ntD, ntG, ntS, ntB, ntTipo, &netlist[ne].L, &netlist[ne].W, &netlist[ne].K,
	  																	&netlist[ne].VT, &netlist[ne].LAMBDA, &netlist[ne].GAMMA, &netlist[ne].THETA, &netlist[ne].LD);
      printf("%s %s %s %s %s %s %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e\n",netlist[ne].nome,ntD, ntG, ntS, ntB, ntTipo, netlist[ne].L, netlist[ne].W, netlist[ne].K,
	  																	netlist[ne].VT, netlist[ne].LAMBDA, netlist[ne].GAMMA, netlist[ne].THETA, netlist[ne].LD);
      netlist[ne].pnMOS = ((ntTipo[0] =='N')?nmos:pmos);
      netlist[ne].tD=numero (ntD);
      netlist[ne].tG=numero (ntG);
      netlist[ne].tS=numero (ntS);
      netlist[ne].tB=numero (ntB);
      netlist[ne].operacaoTransistor = corte;
      double t = (netlist[ne].pnMOS == nmos?0.05:0.02);
      netlist[ne].ALPHA = 2* netlist[ne].K/t;
      linear = false;
	}
    else if (tipo=='*') { /* Comentario comeca com "*" */
      printf("Comentario: %s",txt);
      ne--;
    }
    else {
      printf("Elemento desconhecido: %s\n",txt);
      getch();
      exit(1);
    }
  }
  fclose(arquivo);



  /* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */
  nn=nv;
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo =='L') {
      nv++;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv],"j"); /* Tem espaco para mais dois caracteres */
      strcat(lista[nv],netlist[i].nome);
      netlist[i].x=nv;
    }
    else if (tipo=='H') {
      nv=nv+2;
      if (nv>MAX_NOS) {
        printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
        exit(1);
      }
      strcpy(lista[nv-1],"jx"); strcat(lista[nv-1],netlist[i].nome);
      netlist[i].x=nv-1;
      strcpy(lista[nv],"jy"); strcat(lista[nv],netlist[i].nome);
      netlist[i].y=nv;
    }
  }
 // getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  //getch();
  printf("Netlist interno final\n");
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='R' || tipo=='I' || tipo=='V') {
      printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
    }
	else if(tipo == 'C'){
		printf("%s %d %d %g\n","RC",netlist[i].a,netlist[i].b,netlist[i].valor * FATORDC);
	}
	else if(tipo == 'L'){
		printf("%s %d %d %g\n","RL",netlist[i].a,netlist[i].b,netlist[i].valor / FATORDC);
	}
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
    }
    else if (tipo=='O') {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
      printf("Corrente jx: %d\n",netlist[i].x);
    else if (tipo=='H')
      printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
	}
}

int main(void)
{
	
  srand(time(NULL));
  system ("cls");
  printf("Analise nodal modificada\n");
  denovo:
  /* ############################## Leitura do NETLIST ############################## */
  ne=0; nv=0; strcpy(lista[0],"0");
  printf("Nome do arquivo com o netlist (ex: mna.net): ");
  scanf("%50s",nomearquivo);
  arquivo=fopen(nomearquivo,"r");
  if (arquivo==0) {
    printf("Arquivo %s inexistente\n",nomearquivo);
    goto denovo;
  }

  for(int indice = 0; indice < MAX_NOS; indice ++)
	{
		vAtual[indice] = 0.1;
		vProximo[indice] = 0;
	}

  lerNetlist();

   int numMOS = 0;

	while(!convergiu){ // montar sistema nodal modificado -- MODIFIQUEI!!! (ESTAVA ANTES DO FOR QUE INICIALIZA Yn)

		for (int indice=0; indice<=nv; indice++)   //inicializa os vetores utilizados na analise de convergencia
			for (int j=0; j<=nv+1; j++)
				Yn[indice][j]=0;

		montarEstampa();
							#ifdef DEBUG
								/* Opcional: Mostra o sistema apos a montagem da estampa */
								printf("Sistema apos a estampa de %s\n",netlist[i].nome);
								for (k=1; k<=nv; k++) {
								  for (j=1; j<=nv+1; j++)
									if (Yn[k][j]!=0) {
										printf("%+3.1e ",Yn[k][j]);
									}
									else printf(" ........ ");
								  printf("\n");
								}
								getch();
							#endif
		/* Resolve o sistema */
		if (resolversistema_DC())
			return 1;
							#ifdef DEBUG
							  /* Opcional: Mostra o sistema resolvido */
							  printf("Sistema resolvido:\n");
							  for (i=1; i<=nv; i++) {
								  for (j=1; j<=nv+1; j++){
									if ((Yn[i][j] >= -1/FATORDC*1000) && (Yn[i][j] <= 1/FATORDC*1000)) {
									  Yn[i][j] *= 0;
									}
									if (Yn[i][j]!=0) {
										printf("%+3.1g ",Yn[i][j]);
									}
									else{
										printf(" ... ");
									}
								  }
								  printf("\n");
							  }
							  getch();
							#endif
		/* Mostra solucao */
		printf("Solucao:\n");
		strcpy(txt,"Tensao");
		for (i=1; i<=nv; i++) {
			if (i==nn+1) strcpy(txt,"Corrente");
				printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
			vProximo[i] = Yn[i] [nv+1];
		}
		vezes++;
		//Caso nao linear, utilizar Newton-Raphson
		if(!linear)
			controleConvergencia(vAtual,vProximo,iteracoes);
		else
		  convergiu = true;
	}

	for (int i = 1; i < ne; i++){
		tipo = netlist[i].nome[0];
		if(tipo == 'M'){
			calculoCapacitanciasParasitas(&netlist[i]);
			printf("Gm = %e Gds= %e Gmb = %e\n",netlist[i].gm,netlist[i].gDS,netlist[i].gmB);
			printf("Cgs=%e Cgd= %e Cgb = %e\n",netlist[i].cgS,netlist[i].cgD,netlist[i].cgB);
		}
	}

	// inicio da analise de pequenos sinais
	if (convergiu)
	{
		printf ("\n\n\nCONVERGIU ESSA PORRA!!!\n\n\n");

		int limitador = 0;                            //limitar o numero de linhas na tabela que aparecerá na tela.
		printf("f ");
		for(int i = 1; i < nv + 1; i++ )
			printf("%sm %sf ",lista[i],lista[i]);
		printf("\n");



		// escala linear
		if ( !strcmp(escala_frequencia, "LIN") ) {
			double passo = (freqFinalHz - freqInicialHz) / (ptsFreq - 1);
			for ( frequencia_variavel = freqInicialHz;  frequencia_variavel <= freqFinalHz; frequencia_variavel += passo) {
				// resolucao do sistema AC
				montarEstampa ();
				resolversistema_AC();
				if(limitador <= LIMITE_LINHAS){
					printf("%g " ,frequencia_variavel);
					for(int coluna = 1; coluna < nv + 1 ; coluna++)
						printf("%g %g ",abs(YnComplex[coluna][nv+1]), (180.0/ PI) * arg(YnComplex[coluna][nv+1]));
					printf("\n");
				limitador++;
				} 
				escreverTabelaFreq();				
			}
		}
		// escala logaritmica
		else if ( !strcmp(escala_frequencia, "LOG") ) {
			double passo = 1.0 / (ptsFreq - 1);
			for ( frequencia_variavel = freqInicialHz;  frequencia_variavel <= freqFinalHz; frequencia_variavel *= pow(10, passo)) {
				// resolucao do sistema AC
				montarEstampa ();
				resolversistema_AC();
				if(limitador <= LIMITE_LINHAS){
					printf("%.2e " ,frequencia_variavel);
					for(int coluna = 0; coluna < nv + 1 ; coluna++)
						printf("%.2e %.2e ",abs(YnComplex[i][nv+1]), (180.0/ PI) * arg(YnComplex[i][nv+1]));
					printf("/n");
				}
				escreverTabelaFreq();
			}
			fclose(tabelaFreq);
		}
		// escala octal
		else if ( !strcmp(escala_frequencia, "OCT") ) {
			double passo = 1.0 / (ptsFreq - 1);
			for ( frequencia_variavel = freqInicialHz;  frequencia_variavel <= freqFinalHz; frequencia_variavel *= pow(2, passo)) {
				// resolucao do sistema AC
				montarEstampa ();
				resolversistema_AC();
				if(limitador <= LIMITE_LINHAS){
					printf("%.2e " ,frequencia_variavel);
					for(int coluna = 0; coluna < nv + 1 ; coluna++)
						printf("%.2e %.2e ",abs(YnComplex[i][nv+1]), (180.0/ PI) * arg(YnComplex[i][nv+1]));
					printf("/n");
				}
				escreverTabelaFreq();
			}
			fclose(tabelaFreq);
		}
	}
	//getch();
	
	return 0;
}
