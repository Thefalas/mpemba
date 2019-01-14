### Nombre del compilador
CC=icpc	

# librerias para el enlazador
#mkl version 2011
LPMKL=/opt/intel/composerxe/mkl/lib/intel64

## includes para el compilador
IMKL=/opt/intel/composerxe/include # version 2011

### Opciones del compilador
CFLAGS=-c -g -O3 #-vec-report

#para mkl 11.1, estatico, enteros de 32 bits
#LDFLAGS=-static -L$(LPMKL) -I$(IMKL) -Wl,--start-group $(LPMKL)/libmkl_intel_lp64.a $(LPMKL)/libmkl_intel_thread.a $(LPMKL)/libmkl_core.a -Wl,--end-group -openmp -lpthread

#linkado dinamico, version en paralelo (thread). enteros de 64 bits
LDFLAGS=-L$(LPMKL) -I$(IMKL)-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm 

### Objetos a compilar
OBJS= main.o aleat.o libre.o termo.o mensajes.o colision.o salida.o grafica.o


RWN: $(OBJS)
	$(CC) $(OBJS) -o RWN
#	$(CC) -o RWN -static $(OBJS)

main.o: ../main.cpp ../params.h
	$(CC) $(CFLAGS) ../main.cpp

aleat.o: ../aleat.cpp ../params.h 
	$(CC) $(CFLAGS) ../aleat.cpp 

mensajes.o: ../mensajes.cpp ../params.h 
	$(CC) $(CFLAGS) ../mensajes.cpp 

auxiliar.o: ../auxiliar.cpp ../params.h
	$(CC) $(CFLAGS) ../auxiliar.cpp  

termo.o: ../termo.cpp ../params.h
	$(CC) $(CFLAGS) ../termo.cpp 

libre.o: ../libre.cpp ../params.h
	$(CC) $(CFLAGS) ../libre.cpp

colision.o: ../colision.cpp ../params.h
	$(CC) $(CFLAGS) ../colision.cpp  

salida.o: ../salida.cpp ../params.h 
	$(CC) $(CFLAGS) ../salida.cpp

grafica.o: ../grafica.cpp ../params.h
	$(CC) $(CFLAGS) ../grafica.cpp

	 
all: RWN $(OBJS)

.PHONY: clean
clean: 
	rm RWN $(OBJS)

