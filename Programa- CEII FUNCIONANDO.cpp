/* Análise nodal modificada de:
-  Ponto de Operação
-  Análise no Estado Permanente

Elementos aceitos:
Resistor: R<nome> <nó +> <nó -> <Resistência>
Indutor: L<nome> <nó +> <nó -> <Indutância>
Acoplamento entre indutores: K<nome> <La> <Lb> <k> (La e Lb nomes de indutores já declarados.)
Capacitor: C<nome> <nó +> <nó -> <Capacitância>
Fonte de tensão controlada a tensão: E<nome> <nó V+> <nó V-> <nó v+> <nó v-> <Av>
Fonte de corrente controlada a corrente: F<nome> <nó I+> <nó I-> <nó i+> <nó i-> <Ai>
Fonte de corrente controlada a tensão: G<nome> <nó I+> <nó I-> <nó v+> <nó v-> <Gm>
Fonte de tensão controlada a corrente: H<nome> <nó V+> <nó V-> <nó i+> <nó i-> <Rm>
Fonte de corrente: I<nome> <nó +> <nó -> <módulo> <fase (graus)> <valor contínuo>
Fonte de tensão: V<nome> <nó +> <nó -> <módulo> <fase (graus)> <valor contínuo>
Amplificador operacional ideal: O<nome> <nó saída +> <nó saída -> <nó entrada +> <nó entrada ->
Transistor MOS: M<nome> <nó drain> <nó gate> <nó source> <nó base> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt 0> <lambda> <gamma> <theta> <Ld>
*/

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <complex>
#include <iostream>
#include <time.h>

#define DEBUG

#define LIMITE_LINHAS			5
#define MAX_LINHA 				80
#define MAX_NOME 				11
#define MAX_ELEM 				50
#define MAX_NOS 				50
#define TOLG_DC 				1e-9
#define TOLG_AC 				1e-30
#define FATORDC 				1e9
#define MAX_ITER_FAIL 			1000000 			// numero maximo de iteracoes que a funcao controleConvergencia pode realizar ate abortar o codigo
#define MIN_ITER 				2 					// numero minimo de iteracoes para considerar como estavel a solucao
#define TRY_ITER 				1000 			 	// numero de tentativas minimas para substiuicao dos valores de vAtual que nao convergiram
#define MAX_ERROR 				1e-9				// valor maximo de erro aceitavel para considerar solucao convergente
#define REF_ERR					1 					// valor de referencia utilizado nos calculos de convergencias
#define J 						doubleComplex(0.0,1.0)		// definindo J para facilitar exibicao do codigo
#define PI 						3.141592653589793
#define VALOR_UM 	 			0.9999999999999999999999999999999999999999
#define VALOR_ZERO 	 			0.0000000000000000000000000000000000000001

using namespace std;

enum pontoOperacao{corte, triodo, saturacao};

enum tipoMOS{nmos,pmos};

enum escolha_exibicao{basico,estampas_completo,resultados_completo,tudo};

typedef std::complex<double> doubleComplex;

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
escolha_exibicao varChoice = estampas_completo;

int
  ne, 												/* N?mero de Elementos */
  nv, 												/* N?mero de Variaveis */
  nn, 												/* N?mero de N?s */
  i,j,k;

char

/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo,
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  ntD[MAX_NOME], ntG[MAX_NOME], ntS[MAX_NOME], ntB[MAX_NOME],
  ntTipo[MAX_NOME],nL[MAX_NOME],nW[MAX_NOME],nK[MAX_NOME],nVT[MAX_NOME],nLAMBDA[MAX_NOME],nGAMMA[MAX_NOME],nTHETA[MAX_NOME],nLD[MAX_NOME], //Vari?veis extras para o transistor
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;

char temp_char[MAX_LINHA+1];					// variavel que guarda o nome do arquivo .tab  
  
FILE *arquivo;
FILE *tabelaFreq;

doubleComplex gComplex;
doubleComplex YnComplex[MAX_NOS+1][MAX_NOS+2];

double
  g = 0,
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


bool first_run = true;   			// Primeira vez que roda a analise de pequenos sinais (escrever a primeira linha de informacoes da tabela)
bool frequenciaHz = true; 			// Utiliza a frequencia_variavel em Hz por default
bool linear = true;
bool flagDC = true; 				// Para comecar pela analise do ponto de operacao
bool convergiu = false; 			// Comecar com false para entrar na primeira checagem
bool switch_vAtual = false;			// determina se é necessário trocar o valor de vAtual 

double rand_adjust = 0;				// variavel que ajusta o valor aleatorio dado a vAtual
int rand_iter = 0;					// numero de iteracoes que a funcao controleConvergencia realiza (sem troca aleatoria de vAtual)
int max_rand_iter = 0;				// numero de iteracoes que a funcao controleConvergencia realiza
int count_conv = 0; 				// Contador para convergencia
int restart_counter = 0;			// Contador para quantidade de vezes que o controle de convergencia reiniciou com valores aleatorios
		
char escala_frequencia[3] = {'L','I','N'} ;				// default para escolha de escala para frequencias
int freqInicialHz = 1 ;		 					// default de frequencia inicial para análise de pequenos sinais
int freqFinalHz = 10000 ;						// default de frequencia final para análise de pequenos sinais
int ptsFreq = 100 ;								// default para quantidade de pontos obtidos entre limites de frequencia_variavel

 
// Protótipos das funções
 
double sinDouble (double angulo);
double cosDouble (double angulo);
int resolversistema_DC (void);
int resolversistema_AC (void);
int numero (char *nome);
void controleConvergencia (double vAtual[], double vProximo[]);
void montarEstampa();
void lerNetlist();
void lerArquivo();
void operacoes_DC();
void operacoes_AC();
void calculoCapacitanciasParasitas(elemento netlist[]);
int procuraIndutorTransformador(char *nomeElemento);



double sinDouble (double angulo)
 {
     double seno = sin( angulo * (PI/ 180.0));
    if (fabs(seno) < VALOR_ZERO)
         return (0.0);
    else if (fabs(seno) > VALOR_UM)
         return (1.0);
 	else
     	return (seno);
 }

double cosDouble (double angulo)
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
     printf("Um dos componentes do transformador nao e um indutor\n");
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

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
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
		switch_vAtual = true;
		return 0;
      //printf("Sistema singular\n");
      //return 1;
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

void controleConvergencia ( double vAtual[], double vProximo[] )
{
	double maxErr = 0;
                                                              
    for (int counter_var = 0; counter_var <= nv; counter_var++){
		//define o modo como o erro sera tratado
		if (fabs(vProximo[counter_var]) > REF_ERR)
			tmp_err[counter_var] = fabs( ( vProximo[counter_var] - vAtual[counter_var] ) / vProximo[counter_var]);
		if (vProximo[counter_var] < REF_ERR)
            tmp_err[counter_var] = fabs(vProximo[counter_var] - vAtual[counter_var]);
        if (tempVar[counter_var] > maxErr)
            maxErr = tempVar[counter_var];
    }      

    //faz com o que as tensoes futuras sejam as atuais
	for (int counter_var = 0; counter_var <= nv; counter_var++)
		vAtual[counter_var] = vProximo[counter_var];
        
    //printf("Erro: %f", MAX_ERRO);
    // printf("\nErro atual: %.10f\n", maxErr);
         
	if (maxErr <= MAX_ERROR)
		count_conv++;
    if (count_conv >= MIN_ITER){
		convergiu = true;
		flagDC = false;
		printf("Número de iteracoes: %d || Reinicios: %d\n",max_rand_iter,restart_counter);
		printf("Convergiu! Ponto de Operacao encontrado!\n");
    }
    
	rand_iter++;
	max_rand_iter++;
    if (max_rand_iter >= MAX_ITER_FAIL) {
		printf ("Passou do numero maximo de iteracoes e nao convergiu");
		getch();
		exit(1);
	}
	if (rand_iter >= TRY_ITER || switch_vAtual) {
		for (int counter_var = 0; counter_var <= nv; counter_var++)
			if (tmp_err[counter_var] >= MAX_ERROR){
				vAtual [counter_var] = (rand () % (int)rand_adjust) - rand_adjust/2;   //Utiliza um valor aleatorio entre 0 e 20
				if (vAtual[counter_var] == 0)
					vAtual[counter_var] = 0.1;
				vProximo [counter_var] = 0;
				restart_counter++;
			}
		rand_iter = 0;
		if (rand_adjust < 20) 
			rand_adjust += 0.5;
		if (switch_vAtual)
			max_rand_iter--;
		switch_vAtual = false;
    } 
}

 void montarEstampa()
{
	  /* Monta o sistema nodal modificado */
	  /* Zera sistema */
	  for (i=0; i<=nv; i++) {
		for (j=0; j<=nv+1; j++){
		  Yn[i][j]=0;                              //Inicializacao dos coeficientes
		  YnComplex[i][j] = 0.0 + 0.0 * J;

		}
	  }
	  /* Monta estampas */
	for (i=1; i<=ne; i++) {
		tipo=netlist[i].nome[0];
		
		if (tipo=='R') {						//Monta estampa para o RESISTOR, dependendo do flagDC monta estampa DC ou AC.
			if(flagDC){
			  g= 1.0 / netlist[i].valor;
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
				g = 1.0 / FATORDC;
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
					gComplex = 0.0  + J * (2.0 * PI * netlist[i].valor * (frequencia_variavel));			//Para frequencia_variavel em Hz	    (REVISAR)	
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
      			Yn[netlist[i].x][netlist[i].b]+=1;
      			Yn[netlist[i].x][netlist[i].a]-=1;
      			Yn[netlist[i].x][nv+1]-=netlist[i].valor;
			}else {
				gComplex = netlist[i].modulo * cosDouble(netlist[i].fase) + J * netlist[i].modulo * sinDouble(netlist[i].fase);
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
		 else if (tipo=='M'){                                    //Monta estampa para o TRANSISTOR PMOS OU NMOS,
            if (flagDC) {
                g = ((double) 1.0) / FATORDC;
                netlist[i].gm = 0.0;
                netlist[i].gDS = 0.0;
				netlist[i].gmB =0.0;
                iO = 0.0;
                printf("vd %f vs %f\n", vAtual[netlist[i].tD],vAtual[netlist[i].tS]);
				
				 double vDS = vAtual[netlist[i].tD] - vAtual[netlist[i].tS];

                if((vDS < 0 && netlist[i].pnMOS == nmos) || (vDS > 0 && netlist[i].pnMOS == pmos)) {
                    printf("Tensao no Drain > Tensao no Source\n");
                    int aux = netlist[i].tD;
                    netlist[i].tD = netlist[i].tS;
                    netlist[i].tS = aux;
                }

               
                double vGS = vAtual[netlist[i].tG] - vAtual[netlist[i].tS];
				vDS = vAtual[netlist[i].tD]-vAtual[netlist[i].tS];
                double vBS = vAtual[netlist[i].tB] - vAtual[netlist[i].tS];
				
                if(netlist[i].pnMOS == pmos){
                    vDS *= -1;
                    vBS *= -1;
                    vGS *= -1;
                }
				
				vBS = (vBS>netlist[i].THETA/2.0?netlist[i].THETA/2.0:vBS);

                double vt = netlist[i].VT + netlist[i].GAMMA * (sqrt((netlist[i].THETA - vBS)) - sqrt((netlist[i].THETA)));
                printf ("vt: %f  vGS: %f vDS: %f vS: %f tS %d\n", vt,vGS,vDS, vAtual[netlist[i].tS], netlist[i].tS);
                if (vGS < vt){
                    netlist[i].operacaoTransistor = corte;
                    printf ("Modo de operacao: corte\n");
                }
                else if (vDS < vGS - vt){
                    netlist[i].operacaoTransistor = triodo;
                    printf("Modo de operacao: triodo\n");
                    netlist[i].gm = netlist[i].K * (netlist[i].W/netlist[i].L)*(2.0* vDS)*(1.0+ netlist[i].LAMBDA* vDS);
                    netlist[i].gDS = netlist[i].K * (netlist[i].W/netlist[i].L) * (2.0*(vGS - vt) - 2.0 * vDS + 4.0* netlist[i].LAMBDA * ( vGS - vt) * (vDS) - 3.0*netlist[i].LAMBDA * pow ( vDS,2.0));
                    iO = netlist[i].K * (netlist[i].W/netlist[i].L) * (2.0* (vGS - vt)*vDS - pow (vDS,2.0)) - (netlist[i].gm * vGS) - (gDS * vDS);
                }
                else {
                    netlist[i].operacaoTransistor = saturacao;
                    printf("Modo de operacao: saturacao\n");
                    netlist[i].gm = netlist[i].K * (netlist[i].W/netlist[i].L)* 2.0 *(vGS - vt) * (1.0 + netlist[i].LAMBDA* vDS);
                    netlist[i].gDS = netlist[i].K * (netlist[i].W/netlist[i].L)* pow ((vGS - vt),2.0) * netlist[i].LAMBDA;
                    iO = netlist[i].K * (netlist[i].W/netlist[i].L) * pow((vGS - vt),2.0) * (1.0 + netlist[i].LAMBDA * vDS)- (netlist[i].gm * vGS) - (netlist[i].gDS * vDS);
                }

                netlist[i].gmB = (netlist[i].gm*netlist[i].GAMMA)/(2.0*sqrt(fabs(netlist[i].THETA - (vBS))));
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
		
		
		if (varChoice == estampas_completo || varChoice == tudo)
			if (netlist[i].nome[0] != '.'){
				printf("Sistema apos estampa de: %s\n",netlist[i].nome);
				if (flagDC) {
					for (k=1; k<=nv; k++){
						for (j=1; j<=nv+1; j++)
							if (Yn[k][j]!=0) 
								printf("%+3.1e ",Yn[k][j]);
							else 
								printf(" ........ ");
						printf("\n");
					}
				}
				/*else { // consertar e escolher se vai ser modulo/fase!!!
					for (k=1; k<=nv; k++){
						for (j=1; j<=nv+1; j++)
							if (Yn[k][j]!=0) 
								printf("%+3.1e ",Yn[k][j]);
							else 
								printf(" ........ ");
						printf("\n");
					}
				}*/
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
      printf("%s %s %s %g %g %g\n",netlist[ne].nome,na,nb,netlist[ne].modulo,netlist[ne].fase,netlist[ne].valor);
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
    if (tipo=='R') {
      	printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
    }
    else if (tipo=='I' || tipo=='V') {
      	printf("%s %d %d %g %g %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].modulo,netlist[i].fase,netlist[i].valor);
	}
	else if(tipo == 'C' || tipo == 'L'){
		printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
	}
    else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H') {
      printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
    }
	else if (tipo =='M'){
		printf("%s %s %s %s %s L=%lf W=%lf %lf %lf %lf %lf %lf %lf",netlist[ne].nome,ntD, ntG, ntS, ntB, ntTipo, netlist[ne].L, netlist[ne].W, netlist[ne].K,
	  																	netlist[ne].VT, netlist[ne].LAMBDA, netlist[ne].GAMMA, netlist[ne].THETA, netlist[ne].LD);
	}    else if (tipo=='O') {
      printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
    }
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
      printf("Corrente jx: %d\n",netlist[i].x);
    else if (tipo=='H')
      printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
	}

	printf("\nO circuito tem %d nos, %d variaveis e %d elementos.\nCalculando...\n",nn,nv,ne);
}

void lerArquivo ()
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

}

void operacoes_DC ()
{
	while(!convergiu){ // montar sistema nodal modificado -- MODIFIQUEI!!! (ESTAVA ANTES DO FOR QUE INICIALIZA Yn)

		//for (int indice=0; indice<=nv; indice++)   //inicializa os vetores utilizados na analise de convergencia
			//for (int j=0; j<=nv+1; j++)
				//Yn[indice][j]=0;

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
			exit(1);
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
							
		for (i=1; i<=nv; i++) 
				vProximo[i] = Yn[i] [nv+1];
		
		// Caso nao linear, utilizar Newton-Raphson
		if(!linear)
			controleConvergencia(vAtual,vProximo);
		else {
			convergiu = true;
			flagDC = false;
		}
		  
	}
}

void operacoes_AC ()
{
	
	/* Mostra solucao */
	printf("--------------------\n\n\n\nSolucao:\n");
	strcpy(txt,"Tensao");
	for (i=1; i<=nv; i++) {
		if (i==nn+1) strcpy(txt,"Corrente");
			printf("%s %s: %g\n",txt,lista[i],Yn[i][nv+1]);
		vProximo[i] = Yn[i] [nv+1];
	}	
		
	for (int i = 1; i < ne; i++){
		tipo = netlist[i].nome[0];
		if(tipo == 'M'){
			calculoCapacitanciasParasitas(&netlist[i]);
			if (varChoice != basico){
				printf("Gm = %e Gds= %e Gmb = %e\n",netlist[i].gm,netlist[i].gDS,netlist[i].gmB);
				printf("Cgs=%e Cgd= %e Cgb = %e\n",netlist[i].cgS,netlist[i].cgD,netlist[i].cgB);
			}
		}
	}

	// inicio da analise AC
	if (convergiu)
	{
		printf ("\n\n\n\n");

		int limitador = 0;                            //limitar o numero de linhas na tabela que aparecer� na tela.
		if(first_run)
			printf("Gerando %s , aguarde.\n--------------------\n", temp_char);
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
				if((varChoice == resultados_completo) || (varChoice == tudo)){
					printf("%g " ,frequencia_variavel);
					for(int coluna = 1; coluna < nv + 1 ; coluna++)
						printf("%g %g ",abs(YnComplex[coluna][nv+1]), (180.0/ PI) * arg(YnComplex[coluna][nv+1]));
					printf("\n");
				}
				/*
				if(limitador <= LIMITE_LINHAS){
					printf("%g " ,frequencia_variavel);
					for(int coluna = 1; coluna < nv + 1 ; coluna++)
						printf("%g %g ",abs(YnComplex[coluna][nv+1]), (180.0/ PI) * arg(YnComplex[coluna][nv+1]));
					printf("\n");
				limitador++; 
				}
				*/
				escreverTabelaFreq();
			}
			printf("--------------------\nTabela gerada com sucesso!.\n", temp_char);
		}
		// escala logaritmica
		else if ( !strcmp(escala_frequencia, "DEC") ) {
			double passo = 1.0 / (ptsFreq - 1);
			for ( frequencia_variavel = freqInicialHz;  frequencia_variavel <= freqFinalHz; frequencia_variavel *= pow(10, passo)) {
				// resolucao do sistema AC
				montarEstampa ();
				resolversistema_AC();
				if((varChoice == resultados_completo) || (varChoice == tudo)){
					printf("%g " ,frequencia_variavel);
					for(int coluna = 1; coluna < nv + 1 ; coluna++)
						printf("%g %g ",abs(YnComplex[coluna][nv+1]), (180.0/ PI) * arg(YnComplex[coluna][nv+1]));
					printf("\n");
				}
				escreverTabelaFreq();
			}
			printf("--------------------\nTabela gerada com sucesso!.\n", temp_char);
		}
		// escala octal
		else if ( !strcmp(escala_frequencia, "OCT") ) {
			double passo = 1.0 / (ptsFreq - 1);
			for ( frequencia_variavel = freqInicialHz;  frequencia_variavel <= freqFinalHz; frequencia_variavel *= pow(2, passo)) {
				// resolucao do sistema AC
				montarEstampa ();
				resolversistema_AC();
				if((varChoice == resultados_completo) || (varChoice == tudo)){
					printf("%g " ,frequencia_variavel);
					for(int coluna = 1; coluna < nv + 1 ; coluna++)
						printf("%g %g ",abs(YnComplex[coluna][nv+1]), (180.0/ PI) * arg(YnComplex[coluna][nv+1]));
					printf("\n");
				}
				escreverTabelaFreq();
			}
			printf("--------------------\nTabela gerada com sucesso!.\n", temp_char);
		}
	}
}


int main(void)
{

	lerArquivo();
	lerNetlist();
	operacoes_DC();
	operacoes_AC();
	
	getch();
	
	return 0;
}

