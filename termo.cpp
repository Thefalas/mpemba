/*
 *  termo.cpp
 *  RHCS
 *
 *  Created by Francisco Vega Reyes on 05/06/14.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */


#include "params.h"


int midemomento(void){

  	for(xi=1;xi<=3;xi++) momento[xi]=0.;
	for (i=1;i<=npart;i++) for(xi=1;xi<=3;xi++) momento[xi]+=v[xi][i];
	for(xi=1;xi<=3;xi++) momento[xi]/=npart;
	return 0;

}


int midetemp(void){

	v2=0.0; //v2 es (d)(T/m)
	for (i=1;i<=npart;i++) for(xi=1;xi<=3;xi++) v2+=pow((v[xi][i]-momento[xi]),2);
	v2/=npart;

	w2=0.0;
	for (i=1;i<=npart;i++) for(xi=1;xi<=3;xi++) w2+=pow(w[xi][i],2);
	w2/=npart;
	
	return 0;
}

/*
int midetempfv(void){
  
  int ibin;
  double ecv, ecw; // energías cineticas de las particulas
	v2=0.0; //v2 es (d)(T/m)
	for (i=1;i<=npart;i++) {
	  ecv=0.; 
	  for(xi=1;xi<=3;xi++) ecv+=pow((v[xi][i]-momento[xi]),2);
	  v2+=ecv;
	  fv[sqrt(ecv)/nbin]++;
	}
	v2/=npart;
	for(ibin=1;ibin<=nbin;ibin++) fv[ibin]/=npart;

	w2=0.0;
	for (i=1;i<=npart;i++) for(xi=1;xi<=3;xi++) w2+=pow(w[xi][i],2);
	w2/=npart;
	
	return 0;
}
*/

int termo(void ){

  double vv, vw, ww;
  int ibin;

  c68 = 0;

  foto+=1;
  tf[foto]=t;

  temp[foto]= (2/3.)*v2;
  tempR[foto]=(2.*ka/3.)*w2;

  // inicializacion de cumulantes
  a2v[foto]=0.0;a2w[foto]=0.0;
  a20[foto]=0.0;a02[foto]=0.0;a11[foto]=0.0;
  KA[foto]=0.0;KB[foto]=0.0; //KA=a_00^1; KB=<cos^2\theta>;
  //inicializa solo si es la primera medida
  if(it==tfumbral) {for(ibin=0;ibin<nbin;ibin++) {fv[ibin]=0.; fw[ibin]=0.;fvwdot[ibin]=0.;fvw[ibin]=0.;
      fcos[ibin]=0.; }}
  
  // calculo de cumulantes
  for (i=1;i<=npart;i++){ 
    vv=v[1][i]*v[1][i]+v[2][i]*v[2][i]+v[3][i]*v[3][i];
    ww=w[1][i]*w[1][i]+w[2][i]*w[2][i]+w[3][i]*w[3][i];
    vw=v[1][i]*w[1][i]+v[2][i]*w[2][i]+v[3][i]*w[3][i];

    if(it>=tfumbral){
      ibin=(int)(sqrt(vv/temp[foto])*ndeltabin); fv[ibin]++;  // actualiza f(c)
      ibin=(int)(sqrt(ww*ka/tempR[foto])*ndeltabin); fw[ibin]++;  // actualiza f(w)
      ibin=(int)(vw*vw*ka*ndeltabin/(temp[foto]*tempR[foto])); fvwdot[ibin]++; // actualiza f(s)
      ibin=(int)(ww*ka*vv*ndeltabin/(temp[foto]*tempR[foto])); fvw[ibin]++;  // actualiza f(r)
      ibin=(int)(vw*vw*ndeltabin/(vv*ww));fcos[ibin]++;
      itfv++; // actualiza el numero de iteraciones en medida de las f's
    }

    a2v[foto]+=vv;a2w[foto]+=ww;
    a20[foto]+=(vv*vv);
    a02[foto]+=(ww*ww);
    a11[foto]+=(vv*ww);

    KA[foto]+=(vw*vw);
    KB[foto]+=((vw*vw)/(ww*vv));
  }

  a2v[foto]/=(npart*temp[foto]);a2w[foto]/=(npart*tempR[foto]);a2w[foto]*=ka;

  // preparacion de a20
  a20[foto]/=(npart*temp[foto]*temp[foto]);// <v^4>

   // preparacion de a02
  a02[foto]/=(npart*tempR[foto]*tempR[foto]); 
  a02[foto]*=(ka*ka); //este escalamiento pasa a w desde S (q es la variable de DSMC)// <w^4>

  // preparacion de a11
  a11[foto]/=(npart*tempR[foto]*temp[foto]); a11[foto]*=ka; // <v^2w^2>
 
  //preparacion y calculo de KA y KB. KA es ahora directamente el a_00^(1)

  KA[foto]/=(npart*tempR[foto]*temp[foto]);KA[foto]*=ka; //<(v·w)^2>
  KA[foto]-=(a11[foto]/3.);
  KA[foto]*=(8./15.); // a_00^(1)
  KB[foto]/=npart;KB[foto]-=(1./3.);KB[foto]*=5.; // este es el b prima definido en las notas

 // calculo de a20 y a02
  a20[foto]*=(4./15);a20[foto]-=1.;
  a02[foto]*=(4./15.);a02[foto]-=1.;

  // calculo de a11, b y b_BP
  //  KA[foto]/=a11[foto]; //linea antigua, a mantener comentada
  a11[foto]*=(4./9.);a11[foto]-=1.;
  //KB[foto]/=npart; //linea antigua, tb a mantener comentada


  return 0;

}



int termo_long(void ){

  double vv, vw, ww;

  c68=1;

  foto+=1;
  tf[foto]=t;

  temp[foto]= (2/3.)*v2;
  tempR[foto]=(2.*ka/3.)*w2;

  // inicializacion de cumulantes
  a2v[foto]=0.0;a2w[foto]=0.0;
  a20[foto]=0.0;a02[foto]=0.0;a11[foto]=0.0;
  KA[foto]=0.0;KB[foto]=0.0;

  a30[foto]=0.;a03[foto]=0.; a21[foto]=0.;a12[foto]=0.;
  a40[foto]=0.;a04[foto]=0.; a31[foto]=0.;a13[foto]=0.;a22[foto]=0.;

  // calculo de cumulantes
  for (i=1;i<=npart;i++){ 
    vv=v[1][i]*v[1][i]+v[2][i]*v[2][i]+v[3][i]*v[3][i];
    ww=w[1][i]*w[1][i]+w[2][i]*w[2][i]+w[3][i]*w[3][i];
    vw=v[1][i]*w[1][i]+v[2][i]*w[2][i]+v[3][i]*w[3][i];
    a2v[foto]+=vv;a2w[foto]+=ww;
    a20[foto]+=(vv*vv);
    a30[foto]+=(vv*vv*vv);
    a02[foto]+=(ww*ww);
    a03[foto]+=(ww*ww*ww);
    a11[foto]+=(vv*ww);
    a21[foto]+=(vv*vv*ww);
    a12[foto]+=(vv*ww*ww);

    a40[foto]+=(vv*vv*vv*vv);
    a04[foto]+=(ww*ww*ww*ww);
    a31[foto]+=(vv*vv*vv*ww);
    a13[foto]+=(vv*ww*ww*ww);
    a22[foto]+=(vv*vv*ww*ww);

    KA[foto]+=(vw*vw);
    KB[foto]+=((vw*vw)/(ww*vv));
  }

  a2v[foto]/=(npart*temp[foto]);a2w[foto]/=(npart*tempR[foto]);a2w[foto]*=ka;

  // preparacion de a40, a30 y a20
  a20[foto]/=(npart*temp[foto]*temp[foto]);// <v^4>
  a30[foto]/=(npart*temp[foto]*temp[foto]*temp[foto]);   // <v^6>
  a40[foto]/=(npart*temp[foto]*temp[foto]*temp[foto]*temp[foto]);   // <v^8>

   // preparacion de a04, a03 y a02
  a02[foto]/=(npart*tempR[foto]*tempR[foto]); 
  a02[foto]*=(ka*ka); //este escalamiento pasa a w desde S (q es la variable de DSMC)// <w^4>
  a03[foto]/=(npart*tempR[foto]*tempR[foto]*tempR[foto]);
  a04[foto]/=(npart*tempR[foto]*tempR[foto]*tempR[foto]*tempR[foto]);
  a03[foto]*=(ka*ka*ka); //idem (escalamiento de w a S) // <w^6>
  a04[foto]*=(ka*ka*ka*ka); //idem (escalamiento de w a S) // <w^8>

  // preparacion de a11 y a21 y a12
  a11[foto]/=(npart*tempR[foto]*temp[foto]); a11[foto]*=ka; // <v^2w^2>
  a21[foto]/=(npart*temp[foto]*temp[foto]*tempR[foto]);a21[foto]*=ka; // <v^4w^2>
  a12[foto]/=(npart*temp[foto]*tempR[foto]*tempR[foto]);a12[foto]*=(ka*ka); //<v^2w^4>
 
  // preparacion de a31, a13 y a22
  // <v^6w^2>
  a31[foto]/=(npart*temp[foto]*temp[foto]*temp[foto]*tempR[foto]);a31[foto]*=ka; 
  // <v^2w^6>
  a13[foto]/=(npart*temp[foto]*tempR[foto]*tempR[foto]*tempR[foto]);
  a13[foto]*=(ka*ka*ka); 
  // <v^4w^4>
  a22[foto]/=(npart*temp[foto]*temp[foto]*tempR[foto]*tempR[foto]);
  a22[foto]*=(ka*ka);

  //preparacion y calculo de KA y KB. KA es ahora directamente el a_00^(1)

  KA[foto]/=(npart*tempR[foto]*temp[foto]);KA[foto]*=ka; //<(v·w)^2>
  KA[foto]-=(a11[foto]/3.);
  KA[foto]*=(8./15.); // a_00^(1)
  KB[foto]/=npart;KB[foto]-=(1./3.);KB[foto]*=5.; // este es el b prima definido en las notas

  // calculo de a40
  a40[foto]*=(16./945);
  a40[foto]+=((8./5)*a20[foto]-(32./105)*a30[foto]-3.);

  // calculo de a04
  a04[foto]*=(16./945);
  a04[foto]+=((8./5)*a02[foto]-(32./105)*a03[foto]-3.);

  // calculo de a31
  a31[foto]*=(16./315);
  a31[foto]+=((4./5)*a20[foto]-(8./105)*a30[foto]+(4./3)*a11[foto]-(8./15)*a21[foto]-3.);

  // calculo de a13
  a13[foto]*=(16./315);
  a13[foto]+=((4./5)*a02[foto]-(8./105)*a03[foto]+(4./3)*a11[foto]-(8./15)*a12[foto]-3.);

  // calculo de a22
  a22[foto]*=((16./225));
  a22[foto]+=((4./15)*a20[foto]+(4./15)*a02[foto]+(16./9)*a11[foto]-(16./45)*a21[foto]-(16./45)*a12[foto]-3.);

  // calculo de a30
  a30[foto]*=(-8./105);
  a30[foto]+=((4./5)*a20[foto]-2.); 

  //calculo de a03
  a03[foto]*=(-8./105); 
  a03[foto]+=((4./5)*a02[foto]-2.); 

  // calculo de a21
  a21[foto]*=(-8./45);
  a21[foto]+=((4./15)*a20[foto]+(8./9)*a11[foto]-2);  

  // calculo de a12
  a12[foto]*=(-8./45);
  a12[foto]+=((4./15)*a02[foto]+(8./9)*a11[foto]-2);  

  // calculo de a20 y a02
  a20[foto]*=(4./15);a20[foto]-=1.;
  a02[foto]*=(4./15.);a02[foto]-=1.;

  // calculo de a11, b y b_BP
  //  KA[foto]/=a11[foto]; //linea antigua, a mantener comentada
  a11[foto]*=(4./9.);a11[foto]-=1.;
  //KB[foto]/=npart; //linea antigua, tb a mantener comentada

  return 0;

}



