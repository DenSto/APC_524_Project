# for serial version:
USE_MPI = 0
# for parallel version:
#USE_MPI = 1
#Note: be sure to make clean after changing USE_MPI!

# run test:
#RUN_TEST = 1
# only compile:
RUN_TEST = 0

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

LCONFIG_DIR = $(SRCDIR)/../libconfig-1.5
LCONFIG_LIB = $(LCONFIG_DIR)/lib/libconfig++.a

HDF5_DIR = $(SRCDIR)/../hdf5-1.8.12
HDF5_LIB = $(HDF5_DIR)/hdf5/lib/libhdf5.a

LIBS = $(LCONFIG_LIB) $(HDF5_LIB)

ifeq ($(RUN_TEST),1)
  all: $(PROGRAM) run_tests
else
  all: $(PROGRAM)
endif

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I ../googletest/include

IO/input.o: IO/input.cpp $(LCONFIG_LIB)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS) -I $(LCONFIG_DIR)/include

$(LCONFIG_LIB):
	cd ..; \
	tar xf libconfig-1.5.tar.gz; \
	cd libconfig-1.5; \
	./configure --prefix=$(LCONFIG_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

$(HDF5_LIB):
	cd ..; \
	[ -d $(HDF5_DIR) ] || tar xfv hdf5-1.8.12.tar.gz; \
	cd $(HDF5_DIR); \
	CC=mpic++ ./configure --enable-parallel --prefix=$(HDF5_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

depend:
	$(CC) -MM $(CFLAGS) *.cpp > .depend

clean:
	rm -f $(PROGRAM) $(OBJ) .depend


run_tests:
ifeq ($(USE_MPI),0)
	cd ../test; \
	make; \
	make run_tests; \
	cd ../src;
endif

-include .depend

