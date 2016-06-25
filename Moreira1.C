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
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
#define DEBUG
#define FATORDC 10e9
#define MAX_ERRO 10e-9			//valor maximo de erro para iteracoes
#define MAX_ITER 50 			//maximo de iteracoes
#define REF_VAL 1 				//valor de referencia utilizado nos calculos de convergencias
#define MIN_ITER_CONV 2 		//numero minimo de iteracoes para considerar como estavel a solucao
#define TRY_CONV 5 				//o numero de tentativas que o algoritmo faz para calcular uma solucao, mesmo nao convergindo

typedef struct elemento { /* Definição de Elemento */
  char nome[MAX_NOME];
  double valor;
  int a,b,c,d, e,f,g,h,i,j,k,l,m,x,y;
} elemento;

elemento netlist[MAX_ELEM]; /* Lista de Elementos -> Netlist */

int
  ne, /* Número de Elementos */
  nv, /* Número de Variaveis */
  nn, /* Número de Nós */
  i,j,k;

char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
  nomearquivo[MAX_LINHA+1],
  tipo, 
  na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
  nee[MAX_NOME],nef[MAX_NOME],neg[MAX_NOME],neh[MAX_NOME],nei[MAX_NOME],nej[MAX_NOME],nek[MAX_NOME],nel[MAX_NOME],nem[MAX_NOME], //Variáveis extras para o transistor
  lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome */
  txt[MAX_LINHA+1],
  *p;
FILE *arquivo;

double
  g,
  Yn[MAX_NOS+1][MAX_NOS+2];

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */
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
		for (int counter_var = 0; counter_var <= nv; cont++)
		{
			//tmp_err = discrepância relativa
			if (fabs(vProximo[counter_var]) > REF_VAL)
				tmp_err = fabs( ( vProximo[counter_var] - vAtual[counter_var] ) / vProximo[counter_var]);

			//tmp_err = discrepância absoluta
			if (vProximo[counter_var] < REF_VAL)
				tmp_err = fabs(vProximo[counter_var] - vAtual[counter_var]);

			//atuliza max_err caso tmp_err seja maior
			if (tmp_err > max_err)
				max_err = tmp_err;
        }

		//faz com o que as tensoes futuras sejam as atuais
		for (int counter_var = 0; counter_var <= nv; counter_var++)
			vAtual[counter_var] = vProximo[counter_var];

		//exibe o erro atual entre o nó atual e o nó futuro
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
		
		// caso "count_NOT_conv" supere TRY_CONV (quantidade maxima de tentativas para convergencia), troca-se o modelo do transistor e são zerados os contadores para iteracoes
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
    if (tipo=='R' || tipo=='I' || tipo=='V' || tipo=='L' || tipo=='C' || tipo=='K') {
      sscanf(p,"%10s%10s%lg",na,nb,&netlist[ne].valor);
      printf("%s %s %s %g\n",netlist[ne].nome,na,nb,netlist[ne].valor);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
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
	  sscanf(p,"%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s",na,nb,nc,nd,nee,nef,neg,neh,nei,nej,nek,nel,nem);
      printf("%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",netlist[ne].nome,na,nb,nc,nd,nee,nef,neg,neh,nei,nej,nek,nel,nem);
      netlist[ne].a=numero(na);
      netlist[ne].b=numero(nb);
      netlist[ne].c=numero(nc);
      netlist[ne].d=numero(nd);
      netlist[ne].e=numero(nee);
      netlist[ne].f=numero(nef);
      netlist[ne].g=numero(neg);
      netlist[ne].h=numero(neh);
      netlist[ne].i=numero(nei);
      netlist[ne].j=numero(nej);
      netlist[ne].k=numero(nek);
      netlist[ne].l=numero(nel);
      netlist[ne].m=numero(nem);
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
    if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O') {
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
  /* Monta o sistema nodal modificado */
  printf("O circuito tem %d nos, %d variaveis e %d elementos\n",nn,nv,ne);
  getch();
  /* Zera sistema */
  for (i=0; i<=nv; i++) {
    for (j=0; j<=nv+1; j++)
      Yn[i][j]=0;
  }
  /* Monta estampas */
  for (i=1; i<=ne; i++) {
    tipo=netlist[i].nome[0];
    if (tipo=='R') {
      g=1/netlist[i].valor;
      Yn[netlist[i].a][netlist[i].a]+=g;
      Yn[netlist[i].b][netlist[i].b]+=g;
      Yn[netlist[i].a][netlist[i].b]-=g;
      Yn[netlist[i].b][netlist[i].a]-=g;
    }
    if (tipo=='C') {
      g=netlist[i].valor / FATORDC;
      Yn[netlist[i].a][netlist[i].a]+=g;
      Yn[netlist[i].b][netlist[i].b]+=g;
      Yn[netlist[i].a][netlist[i].b]-=g;
      Yn[netlist[i].b][netlist[i].a]-=g;
    }
    if (tipo=='L') {
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
    else if (tipo=='O') {
      Yn[netlist[i].a][netlist[i].x]+=1;
      Yn[netlist[i].b][netlist[i].x]-=1;
      Yn[netlist[i].x][netlist[i].c]+=1;
      Yn[netlist[i].x][netlist[i].d]-=1;
    }
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
    getch();
    exit;
  }
#ifdef DEBUG
  /* Opcional: Mostra o sistema resolvido */
  printf("Sistema resolvido:\n");
  for (i=1; i<=nv; i++) {
      for (j=1; j<=nv+1; j++){
      	if ((Yn[i][j] >= -0.000000001) && (Yn[i][j] <= 0.000000001)) {
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
  }
  getch();
  return 0;
}
