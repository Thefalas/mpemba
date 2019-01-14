/*
 *  salida.cpp
 *  RHCS
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "params.h"


// *********************************************************************************************************
//  SALIDA DE DATOS: FUNCION salida
// *********************************************************************************************************

int salidatR(void){

    FILE *archivo1;
    char nombre[72];

    // escribe las temperaturas traslacional, rotacional y total y el ratio en RgradsHCS.dat
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/RgradsHCS.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/RgradsHCS.dat",nsim);
    archivo1=fopen(nombre,"a");
    fprintf(archivo1,"%d \t %g \t %g\t  %g \t %g \t %g \t %g \n", foto, idcol, tf[foto], temp[foto],tempR[foto],tempR[foto]/temp[foto], 0.5*(temp[foto]+tempR[foto]));
    fclose(archivo1);	

    // salida de los cumulantes a20, a11, a02
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/aijHCS.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/aijHCS.dat",nsim);
    archivo1=fopen(nombre,"a");
    fprintf(archivo1,"%d \t %g \t %g\t  %g \t %g \t %g \t %g \t %g \n", foto, idcol, tf[foto], a20[foto], a11[foto], a02[foto], KB[foto], KA[foto]);
    fclose(archivo1);	

    // salida de los cumulantes de segundo orden v2, s2sim%.4i/%.2i/c2HCS.dat",nsim,ai);
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/c2HCS.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/c2HCS.dat",nsim);
    archivo1=fopen(nombre,"a");
    fprintf(archivo1,"%d \t %g \t %g \t %g \t %g \n", foto, idcol, tf[foto], a2v[foto], a2w[foto]);
    fclose(archivo1);	

    if(c68){ // imprime cumulantes de orden 6 y 8 si antes se han calculado

    // salida de los cumulantes de orden superior, orden 6
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/n6HCS.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/n6HCS.dat",nsim);
    archivo1=fopen(nombre,"a");
    fprintf(archivo1,"%d \t %g \t %g \t %g \t %g \t %g \t %g \n", foto, idcol, tf[foto], a30[foto], a03[foto], a21[foto],a12[foto]);
    fclose(archivo1);	

    // salida de los cumulantes de orden superior, orden 8
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/n8HCS.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/n8HCS.dat",nsim);
    archivo1=fopen(nombre,"a");
    fprintf(archivo1,"%d \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", foto, idcol, tf[foto], a40[foto], a04[foto], a22[foto], a31[foto],a13[foto]);
    fclose(archivo1);	
    }

return 0;

}


int salidafv(void){

    FILE *archivo1;
    char nombre[72];

  int ibin;
   // salida de los cumulantes de orden superior, orden 8
if (nai>1) sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/%.2i/fmarginal.dat",nsim,ai);
    else sprintf(nombre, "/home/malopez/mpemba/Datos/sim%.4i/fmarginal.dat",nsim);
    archivo1=fopen(nombre,"w");
    for(ibin=0;ibin<nbin;ibin++) fprintf(archivo1,"%g \t %g \t %g\t %g\t %g\t %g \n", (ibin+0.5)/ndeltabin, fv[ibin], fw[ibin], fvwdot[ibin], fvw[ibin], fcos[ibin]);
    fclose(archivo1);	
  
  return 0;

}

int direcs (void){
  
  char nommensajes[72];

		//	Creamos el directorio donde se guardaran los datos
  //if(ai==1){ // para que no lo esté escribiendo en cada punto ai
	sprintf(nommensajes, "/home/malopez/mpemba/Datos/sim%.4i",nsim);
	mkdir(nommensajes,S_IRWXU);
	// }
		
	if(nai>1) {  // si es una serie en parametro, crea subdirectorios para datos
	sprintf(nommensajes, "/home/malopez/mpemba/Datos/sim%.4i/%.2i",nsim,ai);
	mkdir(nommensajes,S_IRWXU);	
	}
	
		
	return 0;
}

