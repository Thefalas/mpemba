#include <stdlib.h>
#include <stdio.h>
#include "params.h"
#include <string.h>


#define NUM_COMMANDS 2

int grafica(void)
{

  char comando[72];

  if (nai>1) sprintf(comando, "gnuplot -p '/home/fvega/Datos/RoughWN/sim%.4i/%.2i/RgradsHCS.dat' u 2:4 w l",nsim,ai);
  else sprintf(comando, "gnuplot -p '/home/fvega/Datos/RoughWN/sim%.4i/RgradsHCS.dat' u 2:5 w l", nsim);
  
  system("killall gnuplot_qt > /dev/null");
  system("gnuplot -e 'nsim=64; ai=1' -p instant_graph.gnuplot");
  return 0;

}
