/*
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "params.h"

double t,tf[nfoto+2],v[3+1][npart+1],temp[nfoto+3],momento[4], ruido;
double tempR[nfoto+3];

int it,i,xi,j,xj,npcol,col,foto,cfoto,c,nc;
double idcol; // contador de colisiones
double amp,fase,dtm,wr,wrmax,dcol,res,norma,prob,v2,tempi;
double sigmaij[4];

double w[3+1][npart+1], sw,ss,sv,vij,wij,w2,u[3+1];
double vc[3+1],wc[3+1];


int ai;	// nai es el numero de diferentes configuraciones (ej., distinto alfa) a simular
double alfa, beta, fa, fb, fab; //	Este es el parametro bucle de la simulacion (aqui, es la inelasticidad)
double a20[nfoto+3], a11[nfoto+3], a02[nfoto+3],KB[nfoto+3],KA[nfoto+3];
double a2v[nfoto+3],a2w[nfoto+3];
double a30[nfoto+3],a03[nfoto+3],a21[nfoto+3],a12[nfoto+3];
double a40[nfoto+3],a04[nfoto+3],a31[nfoto+3],a13[nfoto+3],a22[nfoto+3];

// magnitudes de la funcion de distribucion
double fv[npart],fw [npart],fvwdot[npart],fvw[npart],fcos[npart];
int itfv;

int c68; // indica si se midieron los cumulantes de orden superior o no

int nsim;

int ipv[6+1]={0,3,2,1,3,2,1}; // indices producto vectorial

// variables de la funcion de distribucion

FILE *armensajes;
char nommensajes[72];
char nomfv[72]; // nombre del archivo de fv/fvr
FILE *farchivo;

char fname[72];
int filed;
long nb;

int LR;


double varruido[] = {0., 1.}; // ruido
double varai[]={0, 0.8}; //alfa
double vartempi[]= {0, 9.807618218181856}; // T ROTACIONAL


// Variable para contar el tiempo de ejecucion
time_t hora0,horaF;


int main(void)
{
  // ************		COMIENZO DE  SIMULACION	************

  time(&hora0);

  sprintf(nommensajes,"numsim.info");
  armensajes=fopen(nommensajes,"r");
  fscanf(armensajes,"%d",&nsim);
  fclose(armensajes);
  ++nsim;
  armensajes=fopen(nommensajes,"w");
  fprintf(armensajes,"%.4d \n",nsim); // actualiza el numero de simulacion
  fclose(armensajes);


  iff=0;		// inicializa el generador de numero aleatorio con semilla predeterminada LR (solo una vez en cada simulacion)
  LR = LRi + nth_prime(3+nsim); // change rndm no. seed
  //printf ("nth prime: %d= %d\n", 3+nsim,  nth_prime(3+nsim));
  LR = LRi;

  for (ai=1;ai<=nai;ai++)
    {
      //	AQUI DEFINIMOS CUAL ES EL PARAMETRO EN EL QUE HACEMOS EL BUCLES,
      for(xi=1;xi<=3;xi++) momento[xi]=0.; // se inicializa el momento
      beta = 0.;
      alfa = varai[ai];
      fa = 0.5*(1+alfa);fb=0.5*(1+beta)/(1+ka);fab=fa-ka*fb;
      ruido = gi2dt * varruido[ai];
      tempi = vartempi[ai];
 
      cout << "************************" << endl;
      cout << "alfa: " << alfa << endl;
      cout << "beta: " << beta << endl;
      cout << "k: " << ka << endl;
      cout <<"npart:" << npart <<endl;
      cout << "xi2: " << ruido/dt << endl;
      cout << "************************" << endl << endl;
	
      if (restart){
	
	/*
	//sprintf(fname, "/home/malopez/Datos/RoughWN/sim0087/%.2i/restart.file", ai);
	sprintf(fname, "/home/malopez/Codigos/RoughWN/Datos/ceta/sim0075/02/restart.file");
	filed = open(fname, O_RDONLY);
	
	nb = read(filed, (char *) &v[1], npart*sizeof(double));
	nb = read(filed, (char *) &v[2], npart*sizeof(double));
	nb = read(filed, (char *) &v[3], npart*sizeof( double));

	nb = read(filed, (char *) &w[1], npart*sizeof(double));
	nb = read(filed, (char *) &w[2], npart*sizeof(double));
	nb = read(filed, (char *) &w[3], npart*sizeof( double));

	
	close(filed);
	*/

	if (pulso) {
	  inicgauss(2.0-1.e-6);
	  inicgaussR(1.e-6);
	  midemomento();midetemp();termo();
	  escala(2.0-1.e-6, 1.e-6);
	}

      }   
      else {

	inicgauss(0.9807618218181844);  // distribucion inicial gaussiana, a temperatura 1.0
	inicgaussR(tempi);  // distribucion inicial gaussiana, a temperatura tempi
	midemomento();midetemp();termo();
	escala(0.9807618218181844, tempi); // asegura temps iniciales exactas
      }
      

      // Inicia paso de tiempo reducido, y prob max de colision
      // Reinicializa los indices de recuento de replica, en las variables que se guardan.
      foto=0; cfoto=0;t=0;
      idcol=0.; // pon contador de colisiones a cero
      it=0;  //iterador de tiempo a cero

      itfv=0; // inicializa la funcion de distribucion. esta accion es necesaria para medir bien f(v)

      midemomento();
      midetemp();termo();


      cout << "T inicial: " << (2./3.)*v2 << endl;
      cout << "TR inicial: " << (2.*ka/3.)*w2 << endl;
      cout << "itfv: "<< itfv << endl;
      mensajesR();
      direcs();salidatR();
	
      // *********************************************************************************************************
      //		FUNCION evolucion: evoluciona el estado inicial
      // *********************************************************************************************************


      //		INICIA BUCLE DE PASOS DE TIEMPO
      for (it=1;it<=nt;++it){
	// tiempo, en unidades de t0 
	t+=dt;

	colision();
	ruidogauss(ruido);

	midemomento(); midetemp(); 
	if (it/(utermo*1.) == it/utermo){
	  termo(); mensajesR();salidatR();
	  //grafica();
	}
	
	//escala();

	/*
	// condicion de enfriamiento kovacks
	if (idcol >= 1.0) {
	  termo(); mensajesR();salidatR();
	  //grafica();
	  break; // fin de protocolo
	}
	*/

      } // FIN BUCLE it

      
      // //normaliza todas las iteraciones de la funcion de distribucion y sacalas a archivo,
      // // (solo si se han llegado a hacer medidas de fv)
      // if(itfv) { 
      // 	for(i=0;i<nbin;i++) {
      // 	  fv[i]=(fv[i]*ndeltabin)/(1.0*itfv*npart);
      // 	  fw[i]=(fw[i]*ndeltabin)/(1.0*itfv*npart);
      // 	  fvwdot[i]=(fvwdot[i]*ndeltabin)/(1.0*itfv*npart);
      // 	  fvw[i]=(fvw[i]*ndeltabin)/(1.0*itfv*npart);
      // 	  fcos[i]=(fcos[i]*ndeltabin/(1.0*itfv*npart));
      // 	}
      // 	salidafv(); 
      // 	printf("itfv: %d\n", itfv);
      // 	cout << "itfv: "<< itfv << endl;
      // }


      // if(microsave){
      // 		sprintf(fname, "/home/malopez/Datos/RoughWN/sim%.4i/%.2i/restart.file", nsim,ai);
      // 	//sprintf(fname, "/home/malopez/Datos/RoughWN/sim%.4i/restart.file", nsim);
      // 	filed = creat(fname, 0755);

      // 	nb = write(filed, (char *) &v[1], npart*sizeof(double));
      // 	nb = write(filed, (char *) &v[2], npart*sizeof(double));
      // 	nb = write(filed, (char *) &v[3], npart*sizeof(double));

      // 	close(filed);

      // }


    } // FIN BUCLE en parametro alpha (ai)

  // Mensajes finales
  cout << endl << "\a\a\a";
  time(&horaF);
  cout << "Tiempo de ejecucion: "<< (float) (horaF-hora0)/60. << " minutos" << endl;
  cout << endl << "Pulsa Intro para terminar";
  while((c=getchar())!=EOF) return 0;

} // FIN DE PROGRAMA PRINCIPAL

