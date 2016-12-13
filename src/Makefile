USE_MPI = 0

ifeq ($(USE_MPI),0)
  CC = g++
  OMPFLAGS = -DUSE_MPI=0
else
  CC = mpic++
  #OMPFLAGS = -fopenmp -DUSE_MPI=1
  OMPFLAGS = -DUSE_MPI=1
endif

CFLAGS= -g -lm -Wall 
PROGRAM = EMOOPIC

include Makefile.obj

SRCDIR = $(PWD)
LCONFIG_DIR = ../libconfig-1.5
LCONFIG = $(LCONFIG_DIR)/lib/libconfig++.a

LIBS = $(LCONFIG)

all: $(PROGRAM) test

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

IO/IO.o: IO/IO.cpp $(LCONFIG)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I $(LCONFIG_DIR)/include

$(LCONFIG):
	cd ..; \
	tar xf libconfig-1.5.tar.gz; \
	cd libconfig-1.5; \
	./configure --prefix=$(PWD)/$(LCONFIG_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)
depend:
	$(CC) -MM $(CFLAGS) *.cpp > .depend

clean:
	rm -f $(PROGRAM) $(OBJ) .depend

test:
	cd ../test; \
	make; \
	./IO_unittests; \
	./oGrid_unittests; \
	cd ../src;

-include .depend

