/* An�lise nodal modificada de:
-  Ponto de Opera��o
-  An�lise no Estado Permanente

Elementos aceitos:
Resistor: R<nome> <n� +> <n� -> <Resist�ncia>
Indutor: L<nome> <n� +> <n� -> <Indut�ncia>
Acoplamento entre indutores: K<nome> <La> <Lb> <k> (La e Lb nomes de indutores j� declarados.)
Capacitor: C<nome> <n� +> <n� -> <Capacit�ncia>
Fonte de tens�o controlada a tens�o: E<nome> <n� V+> <n� V-> <n� v+> <n� v-> <Av>
Fonte de corrente controlada a corrente: F<nome> <n� I+> <n� I-> <n� i+> <n� i-> <Ai>
Fonte de corrente controlada a tens�o: G<nome> <n� I+> <n� I-> <n� v+> <n� v-> <Gm>
Fonte de tens�o controlada a corrente: H<nome> <n� V+> <n� V-> <n� i+> <n� i-> <Rm>
Fonte de corrente: I<nome> <n� +> <n� -> <m�dulo> <fase (graus)> <valor cont�nuo>
Fonte de tens�o: V<nome> <n� +> <n� -> <m�dulo> <fase (graus)> <valor cont�nuo>
Amplificador operacional ideal: O<nome> <n� sa�da +> <n� sa�da -> <n� entrada +> <n� entrada ->
Transistor MOS: M<nome> <n� drain> <n� gate> <n� source> <n� base> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt 0> <lambda> <gamma> <theta> <Ld>
*/

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
#define DEBUG
#define FATORDC 10e12
#define MAX_ERRO 10e-9			//valor maximo de erro para iteracoes
#define MAX_ITER 50 			//maximo de iteracoes
#define REF_VAL 1 				//valor de referencia utilizado nos calculos de convergencias
#define MIN_ITER_CONV 2 		//numero minimo de iteracoes para considerar como estavel a solucao
#define TRY_CONV 5 				//o numero de tentativas que o algoritmo faz para calcular uma solucao, mesmo nao convergindo

#define FATORAC bla

enum pontoOperacao{
	corte, triodo, saturacao
};

enum tipoMOS{ 
	nmos,pmos};

typedef struct elemento { /* Defini��o de Elemento */
  char nome[MAX_NOME];
  double valor;
  double cgb,cgs,cgd;
  int a,b,c,d,x,y,tD,tG,tS,tB;  // Nos dos elementos 
  double L,W,K,VT,LAMBDA,GAMMA,THETA,LD, ALPHA;
  char nomeA[MAX_NOME], nomeB[MAX_NOME], NPMOS[MAX_NOME];
  pontoOperacao operacaoTransistor;
  tipoMOS pnMOS;
} elemento;

elemento netlist[MAX_ELEM]; /* Lista de Elementos -> Netlist */



int
  ne, /* N�mero de Elementos */
  nv, /* N�mero de Variaveis */
  nn, /* N�mero de N�s */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo, 
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  ntD[MAX_NOME], ntG[MAX_NOME], ntS[MAX_NOME], ntB[MAX_NOME],
  ntTipo[MAX_NOME],nL[MAX_NOME],nW[MAX_NOME],nK[MAX_NOME],nVT[MAX_NOME],nLAMBDA[MAX_NOME],nGAMMA[MAX_NOME],nTHETA[MAX_NOME],nLD[MAX_NOME], //Vari�veis extras para o transistor
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

double
  g,
  gm,     //
  gDS,   // // Variaveis para o transistor.
  gmB,  //
  iO,  //
  
  Yn[MAX_NOS+1][MAX_NOS+2],
  vAtual[MAX_NOS+1],
  vProximo[MAX_NOS+1];
  
  pontoOperacao operacaoTransistorAtual [MAX_NOS +1];
  pontoOperacao operacaoTransistorProximo [MAX_NOS +1];
  
  bool linear = true;
  bool correct_model;
  bool pOperacao = true; // Para comecar com analise de ponto de operacao
  bool convergiu = false; // Comecar com false para entrar na primeira checagem
  int count_conv = 0; // Contador para convergencia
  int iteracoes = 0; // Numero de iteracoes
  int count_NOT_conv = 0; // Quantas vezes o algoritmo nao converge
  int vezes = 0;

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
   
   void montarEstampaDC();
   
   
   
   
int resolversistema(void)
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
    if (fabs(t)<TOLG) {
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
}

void controleConvergencia ( double vAtual[], double vProximo[], int iteracoes)
{
	double max_err = 0;
	double tmp_err;
	if (iteracoes < MAX_ITER)
	{
		for (int counter_var = 0; counter_var <= nv; counter_var++)
		{
			//tmp_err = discrep�ncia relativa
			if (fabs(vProximo[counter_var]) > REF_VAL)
				tmp_err = fabs( ( vProximo[counter_var] - vAtual[counter_var] ) / vProximo[counter_var]);

			//tmp_err = discrep�ncia absoluta
			if (vProximo[counter_var] < REF_VAL)
				tmp_err = fabs(vProximo[counter_var] - vAtual[counter_var]);

			//atuliza max_err caso tmp_err seja maior
			if (tmp_err > max_err)
				max_err = tmp_err;
        }

		//faz com o que as tensoes futuras sejam as atuais
		for (int counter_var = 0; counter_var <= nv; counter_var++)
			vAtual[counter_var] = vProximo[counter_var];

		//exibe o erro atual entre o n� atual e o n� futuro
		printf("\nErro atual: %.10f\n", max_err);
		
		// se convergir, mantem-se o modelo do transistor e incrementa-se o contador da quantidade de vezes que o algoritmo convergiu
		if (max_err <= MAX_ERRO)
		{
			correct_model = true;
			count_conv++;
		}
		
		// senao, incrementa-se o contador da quantidade de vezes que o algoritmo NAO convergiu
		else if (max_err >= MAX_ERRO)
			count_NOT_conv++;
		
		// caso "count_NOT_conv" supere TRY_CONV (quantidade maxima de tentativas para convergencia), troca-se o modelo do transistor e s�o zerados os contadores para iteracoes
		if ((max_err >= MAX_ERRO) && (count_NOT_conv >= TRY_CONV))
		{
			correct_model = false;
			count_conv = 0;
			iteracoes = 0;
		}

		// caso "count_conv" supere MIN_ITER_CONV (quantidade minima de iteracoes de convergencias corretas), confirma-se que o modelo converge 
		if (count_conv >= MIN_ITER_CONV)
		{
			convergiu = true;
			printf("Convergiu \n");
		}
		
		if (correct_model)
			iteracoes++;
		
		return;
	}
	
	// se "iteracoes" superar MAX_ITER (quantidade maxima de iteracoes que a funcao controleConvergencia pode realizar), assume-se modelo incorreto e zera-se os contadores de iteracao
	else if (iteracoes >= MAX_ITER)
	{
		correct_model = false;
		count_conv = 0;
		iteracoes = 0;
		
		return;
	}
	 
}

int main(void)
{
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
    if (tipo=='R' || tipo=='I' || tipo=='V' || tipo=='L' || tipo=='C') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
    }
    else if (tipo=='K') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      
      for(int count = 0; count < MAX_NOME; count++){
    	netlist[ne].nomeA[count]=na[count];
      	netlist[ne].nomeB[count]=nb[count];
	  }
	  
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
	else if (tipo=='M'){
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
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='K') {
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
  getch();
  /* Lista tudo */
  printf("Variaveis internas: \n");
  for (i=0; i<=nv; i++)
    printf("%d -> %s\n",i,lista[i]);
  getch();
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
  getch();
  	
	while(!convergiu){
		montarEstampaDC();
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
		  }
		  /* Resolve o sistema */
		  if (resolversistema()) {
			return 1;
		  }
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
	  if(!linear){
		  controleConvergencia(vAtual,vProximo,iteracoes);
	  } else{
		  convergiu = true;
	  }
	  for (int i = 1;i<ne;i++){
		  tipo=netlist[i].nome[0];
		  if(tipo == 'M'){
			  printf("Gm = %e Gds= %e Gmb = %e\n",gm,gDS,gmB);
			  printf("Cgs=%e Cgd= %e Cgb = %e\n",netlist[i].cgs,netlist[i].cgd,netlist[i].cgb);
		  }
	  }
	  
}

 void montarEstampaDC(){
	  /* Monta o sistema nodal modificado */
	  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
	  getch();
	  /* Zera sistema */
	  for (i=0; i<=nv; i++) {
		for (j=0; j<=nv+1; j++)
		  Yn[i][j]=0;
	  }
	  /* Monta estampas */
	  int numMOS = 0;
	  for (i=1; i<=ne; i++) {
		tipo=netlist[i].nome[0];
		if (tipo=='R') {
		  g=1/netlist[i].valor;
		  Yn[netlist[i].a][netlist[i].a]+=g;
		  Yn[netlist[i].b][netlist[i].b]+=g;
		  Yn[netlist[i].a][netlist[i].b]-=g;
		  Yn[netlist[i].b][netlist[i].a]-=g;
		}
		else if (tipo=='K') {
		  g=netlist[i].valor;
		  
		  for (int count1 = 1; count1 <= ne; count1++){
			if (netlist[i].nomeA == netlist[count1].nome){
				Yn[netlist[i].x][netlist[count1].a]+=g;
				Yn[netlist[i].x][netlist[count1].b]-=g;
				Yn[netlist[count1].x][netlist[i].x]-=g;
				Yn[netlist[count1].x][netlist[i].x]+=g;
				break;
			}
		  }
		  
		for (int count2 = 1; count2 <= ne; count2++){
			if (netlist[i].nomeB == netlist[count2].nome){
				
				Yn[netlist[i].x][netlist[count2].a]-=1;
				Yn[netlist[i].x][netlist[count2].b]+=1;
				Yn[netlist[count2].x][netlist[i].x]+=1;
				Yn[netlist[count2].x][netlist[i].x]-=1;
				break;
			}
		  }

		}
		else if (tipo=='C') {
		  g=netlist[i].valor / FATORDC;
		  Yn[netlist[i].a][netlist[i].a]+=g;
		  Yn[netlist[i].b][netlist[i].b]+=g;
		  Yn[netlist[i].a][netlist[i].b]-=g;
		  Yn[netlist[i].b][netlist[i].a]-=g;
		}
		else if (tipo=='L') {
		  g=netlist[i].valor * FATORDC;
		  Yn[netlist[i].a][netlist[i].a]+=g;
		  Yn[netlist[i].b][netlist[i].b]+=g;
		  Yn[netlist[i].a][netlist[i].b]-=g;
		  Yn[netlist[i].b][netlist[i].a]-=g;
		}
		else if (tipo=='G') {
		  g=netlist[i].valor;
		  Yn[netlist[i].a][netlist[i].c]+=g;
		  Yn[netlist[i].b][netlist[i].d]+=g;
		  Yn[netlist[i].a][netlist[i].d]-=g;
		  Yn[netlist[i].b][netlist[i].c]-=g;
		}
		else if (tipo=='I') {
		  g=netlist[i].valor;
		  Yn[netlist[i].a][nv+1]-=g;
		  Yn[netlist[i].b][nv+1]+=g;
		}
		else if (tipo=='V') {
		  Yn[netlist[i].a][netlist[i].x]+=1;
		  Yn[netlist[i].b][netlist[i].x]-=1;
		  Yn[netlist[i].x][netlist[i].a]-=1;
		  Yn[netlist[i].x][netlist[i].b]+=1;
		  Yn[netlist[i].x][nv+1]-=netlist[i].valor;
		}
		else if (tipo=='E') {
		  g=netlist[i].valor;
		  Yn[netlist[i].a][netlist[i].x]+=1;
		  Yn[netlist[i].b][netlist[i].x]-=1;
		  Yn[netlist[i].x][netlist[i].a]-=1;
		  Yn[netlist[i].x][netlist[i].b]+=1;
		  Yn[netlist[i].x][netlist[i].c]+=g;
		  Yn[netlist[i].x][netlist[i].d]-=g;
		}
		else if (tipo=='F') {
		  g=netlist[i].valor;
		  Yn[netlist[i].a][netlist[i].x]+=g;
		  Yn[netlist[i].b][netlist[i].x]-=g;
		  Yn[netlist[i].c][netlist[i].x]+=1;
		  Yn[netlist[i].d][netlist[i].x]-=1;
		  Yn[netlist[i].x][netlist[i].c]-=1;
		  Yn[netlist[i].x][netlist[i].d]+=1;
		}
		else if (tipo=='H') {
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
		}
		else if (tipo=='M'){
			g = ((double) 5.0) / FATORDC;
			gm = 0;
			gDS = 0;
			iO = 0;
			printf("vd %f vs %f\n", vAtual[netlist[i].tD],vAtual[netlist[i].tS]);
			if(vAtual[netlist[i].tD] < vAtual[netlist[i].tS])
			{
				printf("Tensao no Drain > Tensao no Source\n");
				int aux = netlist[i].tD;
				netlist[i].tD = netlist[i].tS;
				netlist[i].tS = aux;
			}
			double vGS = vAtual[netlist[i].tG] - vAtual[netlist[i].tS];
			double vDS = vAtual[netlist[i].tD] - vAtual[netlist[i].tS];
			double vt = netlist[i].VT + netlist[i].GAMMA * (sqrt(fabs(netlist[i].THETA - (vAtual[netlist[i].tB] -vAtual[netlist[i].tS]))) - sqrt(fabs(netlist[i].THETA)));
			printf ("vt: %f  vGS: %f vDS: %f vS: %f tS %d\n", vt,vGS,vDS, vAtual[netlist[i].tS], netlist[i].tS);
			if (vGS < vt)
			  {		
				   netlist[i].operacaoTransistor = corte;
				   printf ("Modo de operacao: corte\n");
			  }
			else if (vDS < vGS - vt)
			{
				netlist[i].operacaoTransistor = triodo;
				printf("Modo de operacao: triodo\n");
				gm = netlist[i].K * (netlist[i].W/netlist[i].L)*(2* (vDS))*(1+ netlist[i].LAMBDA* vDS);
				   gDS = netlist[i].K * (netlist[i].W/netlist[i].L) * (2*(vGS - vt) - 2 * vDS + 4* netlist[i].LAMBDA * ( vGS - vt) * (vDS) - 3*netlist[i].LAMBDA * pow ( vDS,2));
				   iO = netlist[i].K * (netlist[i].W/netlist[i].L) * (2* (vGS - vt)*vDS - pow (vDS,2)) - (gm * vGS) - (gDS * vDS);
			  
			}
			else
			{
				netlist[i].operacaoTransistor = saturacao;
				printf("Modo de operacao: saturacao\n");
				gm = netlist[i].K * (netlist[i].W/netlist[i].L)* 2 *(vGS - vt) * (1 + netlist[i].LAMBDA* vDS);
				
				gDS = netlist[i].K * (netlist[i].W/netlist[i].L)* pow ((vGS - vt),2) * netlist[i].LAMBDA;
				   
				iO = netlist[i].K * (netlist[i].W/netlist[i].L) * pow((vGS - vt),2) * (1 + netlist[i].LAMBDA * vDS) 
					- (gm * vGS) - (gDS * vDS);
			}
			
			double gmB = (gm*netlist[i].GAMMA)/(2*sqrt(fabs(netlist[i].THETA - (vAtual[netlist[i].tB] -vAtual[netlist[i].tS]))));
			
			iO-= (gmB * (vAtual [netlist[i].tB] -vAtual[netlist[i].tS]));
			
			printf ("gm %f gmb %f  gds %f io %f \n", gm, gmB, gDS, iO);
			
			
			  Yn[netlist[i].tD][netlist[i].tB]+=gmB;
			  Yn[netlist[i].tS][netlist[i].tS]+=gmB;
			  Yn[netlist[i].tD][netlist[i].tS]-=gmB;
			  Yn[netlist[i].tS][netlist[i].tB]-=gmB;

			  Yn[netlist[i].tD][netlist[i].tG]+=gm;
			  Yn[netlist[i].tS][netlist[i].tS]+=gm;
			  Yn[netlist[i].tD][netlist[i].tS]-=gm;
			  Yn[netlist[i].tS][netlist[i].tG]-=gm;

			  Yn[netlist[i].tD][netlist[i].tD]+=gDS;
			  Yn[netlist[i].tS][netlist[i].tS]+=gDS;
			  Yn[netlist[i].tD][netlist[i].tS]-=gDS;
			  Yn[netlist[i].tS][netlist[i].tD]-=gDS;


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
		else if (tipo=='O') {
		  Yn[netlist[i].a][netlist[i].x]+=1;
		  Yn[netlist[i].b][netlist[i].x]-=1;
		  Yn[netlist[i].x][netlist[i].c]+=1;
		  Yn[netlist[i].x][netlist[i].d]-=1;
		}
	  }
	}


