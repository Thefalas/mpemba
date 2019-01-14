/*
 *  colision.cpp
 *  RHCS
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "params.h"


//       **** FUNCION colision ****
int colision(void )
{

prob=2.;
//do{

wrmax=dr2*fwr*sqrt((2./3.)*v2);  // actualiza la prob. máxima de colisión // aqui añadi factor 2 ante fwr. REVISAR

//		Numero m�ximo de colisiones
//    	dcol=res+npart*dr2*fwr*dt/2.; // dcol=N omegamax h/2 =
    	dcol=res+npart*sqrt(2.*v2/3.)*dr2*fwr*dt/2.; // dcol=N omegamax h/2 = // inicialmente era 
	npcol=(int)floor(dcol);
        res=dcol-floor(dcol);


for(col=1;col<=npcol;col++)
	{

//		Vector de colision. Genera una direccion de colision aleatoria: las 3 coordenadas de sigmaij
        amp=2.*aleat(LR)-1.;
        fase=2.*PI*aleat(LR);
	sigmaij[1]=sqrt(1.-pow(amp,2))*cos(fase);
        sigmaij[2]=sqrt(1.-pow(amp,2))*sin(fase);
        sigmaij[3]=amp;

   //   Consigue aleatoriamente un par de particulas (i,j)
        i=(int)ceil(aleat(LR)*npart);	// devuelve un numero aleatorio de particula entre 1 y npart
et2:	j=(int)ceil(aleat(LR)*npart);// devuelve un 2o idem
		if (i==j) goto et2; // repite el 2o paso si las particulas coinciden

   //   Calcula wr, la probabilidad de colision del par (i,j), que es g�sigma
		prob=0.;

        wr=(v[1][i]-v[1][j])*sigmaij[1]+(v[2][i]-v[2][j])*sigmaij[2]+(v[3][i]-v[3][j])*sigmaij[3];
	

		if (wr > 0.)	// solo si ambas particulas se aproximan, considera posibilidad de colision
		{
         prob=dr2*wr/wrmax;
         if(prob > 1.)	// corrige wrmax si no suficientemente grande
		 {
          cout << "CUIDADO!, wr= " << wr << "\tmayor que wrmax= " << wrmax << endl ;
          wrmax=wr;
//		  break;
		 }
		}

        if (prob > aleat(LR))	// si la condicion de colision se cumple
		{						// actualiza las velocidades del par de particulas

		  idcol+=deltacol;  //actualiza el numero de colisiones por particula

ss=(w[1][i]+w[1][j])*sigmaij[1]+(w[2][i]+w[2][j])*sigmaij[2]+(w[3][i]+w[3][j])*sigmaij[3]; //REVISAR
		  
		for (xi=1;xi<=3;++xi)
        {


       sw=(w[ipv[2*xi-1]][i]+w[ipv[2*xi-1]][j])*sigmaij[ipv[2*xi]]-(w[ipv[2*xi]][i]+w[ipv[2*xi]][j])*sigmaij[ipv[2*xi-1]];
	  sv=(v[ipv[2*xi-1]][i]-v[ipv[2*xi-1]][j])*sigmaij[ipv[2*xi]]-(v[ipv[2*xi]][i]-v[ipv[2*xi]][j])*sigmaij[ipv[2*xi-1]];
	  vij=v[xi][i]-v[xi][j];wij=w[xi][i]+w[xi][j];
	  // terminos de correccion de las velocidades, se guardan aparte para evitar solapamientos de memoria
	  vc[xi] = fab*wr*sigmaij[xi]+fb*ka*(vij-sw);
	  wc[xi] = fb*(sv-sigmaij[xi]*ss+wij);//+fa*wr*sigmaij[xi]; //REVISAR

	}

		for (xi=1;xi<=3;++xi){
		  v[xi][i] -= vc[xi]; v[xi][j] += vc[xi];
		  w[xi][i] -= wc[xi]; w[xi][j] -= wc[xi];
		}

		} // fin de la condicion de colision


	} //  fin del bucle col (de colisiones)


//}while(prob>1.); // condici�n de asignaci�n condicional de numero de colisiones


		return 0;

}


