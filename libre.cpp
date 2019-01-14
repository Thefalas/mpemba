/*
 *  libre.cpp
 *  RHCS
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "params.h"

// ******************************************************************************************************************
//		FUNCION calcogauss: excita un sistema de npart particulas una temperatura norma, con distribucion gaussiana
// ******************************************************************************************************************


int inicgauss(double norma)
{
		double x1,x2,wi;

		cout << "momento inicial residual (componentes): " << endl;

		for (xi=1;xi<=3;xi++){

		momento[xi]=0.;	// inicializa momento de primer orden


			for (i=2;i<=npart;i+=2)
		 	{


			do{
				x1=2.0*aleat(LR)-1.0;
				x2=2.0*aleat(LR)-1.0;
				wi=x1*x1+x2*x2;
			}while(wi>=1.0);
			//w=sqrt((-2.0*norma*log(w))/w);
			wi=sqrt((-norma*log(wi))/wi);
			v[xi][i-1]=x1*wi;
			v[xi][i]=x2*wi;



			momento[xi]+=(v[xi][i-1]+v[xi][i]);	// $ comienza calculo del momento de 1er orden


			}

		momento[xi]/=npart;	// $ acaba calculo del momento de 1er orden

		cout << momento[xi] << endl;	// presenta el momento de 1er orden residual (velocidad
											// media) de la distribucion gaussiana inicial

		for (i=1;i<=npart;++i) v[xi][i]-=momento[xi]; // $ corrige la distribucion de velocidad para evitar el flujo
													 // comentar lineas marcadas $ para suprimir correccion de momento
		}

     
return 0;

}


int inicgaussR(double normaR)
{
		double x1,x2,wi;

		// se hace la parte rotacional de la temperatura 
		cout << "momento inicial residual ROTACIONAL (componentes): " << endl;

		for (xi=1;xi<=3;xi++){

		momento[xi]=0.;	// inicializa momento de primer orden


			for (i=2;i<=npart;i+=2)
		 	{


			do{
				x1=2.0*aleat(LR)-1.0;
				x2=2.0*aleat(LR)-1.0;
				wi=x1*x1+x2*x2;
			}while(wi>=1.0);
			//w=sqrt((-2.0*norma*log(w))/w);
			wi=sqrt((-normaR*log(wi))/(ka*wi)); //borrar el 0.05
			w[xi][i-1]=x1*wi;
			w[xi][i]=x2*wi;



			momento[xi]+=(w[xi][i-1]+w[xi][i]);	// $ comienza calculo del momento de 1er orden


			}

		momento[xi]/=npart;	// $ acaba calculo del momento de 1er orden

		cout << momento[xi] << endl;	// presenta el momento de 1er orden residual (velocidad
											// media) de la distribucion gaussiana inicial

		for (i=1;i<=npart;++i) w[xi][i]-=momento[xi]; // $ corrige la distribucion de velocidad para evitar el flujo
													 // comentar lineas marcadas $ para suprimir correccion de momento
		}
	
return 0;

}



int ruidogauss(double normaR)
{
		double x1,x2,w;
 
for (xi=1;xi<=3;++xi)
		{


                for (i=2;i<=npart;i+=2)
                 {

				do{
					x1=2.0*aleat(LR)-1.0;
					x2=2.0*aleat(LR)-1.0;
					w=x1*x1+x2*x2;
				}while(w>=1.0);
				w=sqrt((-normaR*log(w))/w);// era 0.5*normaR...

				v[xi][i-1]+=x1*w;
				v[xi][i]+=x2*w;
		            

                  } 


		}

	return 0;
	
}


int prepara_pulso(double reducc){

  for (i=1;i<=npart;i++) for(xi=1;xi<=3;xi++) w[xi][i]*=reducc;
  inicgauss(1.0);
  return 0;

}

//       **** **** **** **** **** **** **** **** ****

// esta funcion escala a valores determinados de T_t y T_r
int escala(double fact_t, double fact_r){

  for (xi=1;xi<=3;xi++) {
    for (i=1;i<=npart;i++) {
      v[xi][i] = sqrt(fact_t)*v[xi][i]/sqrt((2./3.)*v2);
      w[xi][i] = sqrt(fact_r)*w[xi][i]/sqrt((2.*ka/3.)*w2);
    }
  }

  return 0;

}

//       **** **** **** **** **** **** **** **** ****


// esta funcion produce los numeros primos hasta 200

int nth_prime( int i_prime) 
{

  int n_prime,ii;

  n_prime = i_prime;
  if (n_prime < 3){
    if (n_prime == 1 ) i = 2;
    if (n_prime == 2 ) i = 3;
  }
  else{
    ii = 2;
    i = 1;
    while (ii < n_prime){
      i++;
      for (int j=2; j*j<=i; j++)
	{
	  if (i % j == 0) 
	    break;
	  else if (j+1 > sqrt(i)) {
	    //printf(" %d\n", i);
	    ii ++;
	  }

	}
    }
  }

  //  printf("prime number is: %d\n", i);
  return i;

}
