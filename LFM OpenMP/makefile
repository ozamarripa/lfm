#CC=mpic++
CC=g++
#CC=tau_cxx.sh
LOP= -o
#CFLAGS = -g 
#CFLAGS+=-gstabs+
CFLAGS+=-O3
#CFLAGS+=-O2
CFLAGS+=-fopenmp
CFLAGS+=-funroll-loops
MAIN=./Library/lfm
TAG=lfm


$(MAIN).o :
	$(CC) $(LOP) $(TAG) $(CFLAGS) $(MAIN).cpp 


