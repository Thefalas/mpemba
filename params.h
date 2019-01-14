/*
 *  params.h
 *  RHCS. Estado de enfriamiento homogeneo esferas duras rugosas
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

//**********************************************************************************
//*********   SIMULACION MONTE CARLO de  GRANULODINAMICA ESTOCASTICA	   *********
//***************** *****************************************************************


// paquetes C++
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <iostream>
#include<stdlib.h>
#include<fcntl.h>
#include<sys/stat.h>
#include<time.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;



//************************************************************
//*       PARAMETROS DEFINIDOS EN LA SIMULACION
//************************************************************


//************************************************************
//*      Constantes Matematicas
//************************************************************

#define PI acos(-1.)
#define rPI sqrt(PI)
#define r2 sqrt(2.)
#define dr2 2.*r2


//	semilla del generador de numero aleatorio number (aleat.cpp)
#define LRi 102133*10


//************************************************************
//*       Parametros generales
//************************************************************

//	numero de particulas
#define npart 200000//was 2000000 // hasta dos millones de particulas!!
#define deltacol 1./npart
//	numero de pasos de tiempo
#define nt 1000000 // was 800000
// unidad de paso de tiempo
#define dt 0.00025 // was 0.00025 
//	pasos de tiempo entre medidas de temperatura (utermo debe ser menor o igual que nt)
#define utermo 200
// pasos de tiempo entre reescalado de velocidades
#define uescala 1000
//	numero de medidas guardadas de la temperatura
#define nfoto nt/utermo
//	pasos de tiempo entre medidas de las magnitudes colisiones, relativas a utermo
//#define ucol 20*utermo
#define ucol 200
//	numero de medidas guardadas de la temperatura
//#define ncfoto nt/(ucol)
#define ncfoto nt/ucol
// numero de valores del parametro variable
#define nai 1

#define ka 2./5 // k entre 0 y 2/3, 2/5: masa uniforme, 2/3: masa en la superficie; k=0, masa en el centro


//************************************************************
//*      Parametros teoricos y de simulacion
//************************************************************

// RESTART. si 1, lee empieza desde estado microscopico grabado en restart.file.
#define restart 0
// microsave. si 1, guarda estado miscroscopico en restart.file
#define microsave 1

#define pulso 0// 0 si no se quiere pulso calenton (solo funciona si restart = 1)

// fwr
#define fwr 20.5


//*******
#define gi2dt 1.0 * dt // es el nivel de ruido al cuadrado por dt 
//*******

// PARAMETROS PARA LAS FUNCIONES DE DISTRIBUCION MARGINALES
#define tfumbral nt*(1+0.1)//umbral de tiempo para medir f(v)'s. 
//si (tfumbral > nt) no se mide f(v) 
#define vfumbral 10. // umbral de velocidad termica para medir funciones de distribucion marginales
#define ndeltabin  50. // numero de bines por unidad de velocidad termica 
#define nbin vfumbral*ndeltabin // numero de bines en total

//************************************************************************************
//      VARIABLES GENERALES DE LA SIMULACION: v: velocidad de las particulas; momento: momentos generalizados;
//		t: tiempo, temp: temperatura.
//		tsal: tiempo medido en tempsal: temperatura. Estas dos variables son para almacenar.

extern double t,tf[nfoto+2],v[3+1][npart+1],v2, temp[nfoto+3], momento[4], ruido;
extern double tempR[nfoto+3];


extern double w[3+1][npart+1], sw, ss, sv, vij, wij, w2;
extern double vc[3+1],wc[3+1];

//magnitudes de cumulantes
extern double a20[nfoto+3], a11[nfoto+3], a02[nfoto+3], KB[nfoto+3], KA[nfoto+3];
extern double a2w[nfoto+3], a2v[nfoto+3];
extern double a30[nfoto+3], a03[nfoto+3], a21[nfoto+3],a12[nfoto+3];
extern double a40[nfoto+3],a04[nfoto+3],a31[nfoto+3],a13[nfoto+3],a22[nfoto+3];


//************************************************************************************

//		VARIABLES DE ITERACION

extern int it,i,xi,j,xj,npcol,col,foto,cfoto,c,nrcol;
extern double idcol; // contador de colisiones
extern int nsim;
	/*  indices (enteros):

		replica: numero de replica estadistica
		it: paso temporal
		i: numero de particula
		xi: numero de coordenada (1,2 o 3 en 3D)
		j: par de colision de la particula
		npcol: numero total de colisiones en cada paso de tiempo
		col: numero de colision
	*/

extern double amp,fase,dtm,wr,wrmax,dcol,res,prob,sigmaij[4],tempi; // *** REVISAR dim de wc, wgc ***


	/*	magnitudes de doble precision:

		amp: amplitud del vector de velocidad aleatorio
		fase: fase del vector de velocidad aleatorio
		t: tiempo (unidades reducidas)
		dtm: incremento de tiempo (unidades reducidas)
		wr: probabilidad de colision del par, calculada de la direccion de velocidad aleatoria (unidades reducidas)
		wrmax: probabilidad maxima de colision (unidades reducidas)
		dcol: estimacion del numero de colisiones (en doble precision)
		res: diferencia entre dcol y npcol
		temp: temperatura
		sigmaij[4]: direccion de colision normalizada generada

	*/

//		MAGNITUDES double SIN ESTABILIZAR
extern double alfa, beta, varai[],fa, fb, fab;  // estas se usan si se hace bucle en parametro, por ej., alfa
extern int ai; // variable de iteracion, bien de alfa, bien de beta, dependiendo del valor de
extern int ipv[6+1]; // indices producto vectorial

// magnitudes de la funcion de distribucion
extern double fv[npart],fw[npart],fvwdot[npart],fvw[npart],fcos[npart];
extern int itfv;

extern int c68; // variable logica. si =1: se han medido los cumulantes de orden superior, si=0: no.

extern int LR;

//************************************************************************************
// CONJUNTO DE FUNCIONES LLAMADAS POR EL PROGRAMA PRINCIPAL
extern int inicgauss(double norma );	// inicia una distribucion gaussiana a temperatura 'norma'
extern int inicgaussR(double normaR);
extern int colision(void );	// selecciona las coisiones. actualiza las velocidades para las colisiones aceptadas
extern int midetemp(void);	//	mide la temperatura
extern int termo(void );	// mide la temperatura y la guarda en una variable, cada cierto tiempo
extern int termo_long(void );	// idem y ademas mide cumulantes de orden superior
extern int salidatR(void);
extern int direcs(void);
extern int mensajesR(void);
extern int midemomento(void);
//************************************************************************************

//
extern int ruidogauss(double normaR);
//

//
//extern int midetempfv(void); //mide la temperatura Y LA FUNCION DE DISTRIBUCION
extern int salidafv(void); // salida de datos correspondientes a la funcion de distr.

// funcion graficadora
extern int grafica(void);

// funcion escala velocidades, en caso de HCS
extern int escala(double fact_t, double fact_r);

// funcion que inicializa el pulso (estado inicial gaussiano descompensado)
extern int prepara_pulso (double reducc);


//************************************************************************************
//		FUNCIONES ANEXAS (y parametros) AL CODIGO PRINCIPAL y variables externas asociadas
// funciones:
extern double aleat(int idum);		// aleat: generador del numero aleatorio; idum: semilla
											// a: coeficiente de restitucion normal, fracsol: fraccion solida
// parametros:
extern double al;			// el numero aleatorio entre 0 y 1 devuelto por aleat es almacenado en esta variable
extern int iff;				// iff inicializa, si nulo, el generador de numero aleatorio con la semilla 'idum'

extern int nth_prime (int i_prime);
//************************************************************************************


