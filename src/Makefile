CC = g++
MPICC = mpic++
MPIFLAGS = -DUSE_MPI=1

CFLAGS= -g  -Wall -O3

include Makefile.obj

SRCDIR = $(PWD)

GTEST_DIR = $(SRCDIR)/../googletest
GTEST_INCL = $(GTEST_DIR)/include

LCONFIG_DIR = $(SRCDIR)/../libconfig-1.5
LCONFIG_INCL = $(LCONFIG_DIR)/include
LCONFIG_LIB = $(LCONFIG_DIR)/lib/libconfig++.a

HDF5_DIR = $(SRCDIR)/../hdf5-1.8.12
HDF5PAR_DIR = $(SRCDIR)/../hdf5-1.8.12-parallel

HDF5_INCL = $(HDF5_DIR)/include
HDF5PAR_INCL = $(HDF5PAR_DIR)/include

HDF5_LIB = $(HDF5_DIR)/lib/libhdf5.a
HDF5PAR_LIB = $(HDF5PAR_DIR)/lib/libhdf5.a

LIBS = $(LCONFIG_LIB) $(HDF5_LIB)
LIBS_PAR = $(LCONFIG_LIB) $(HDF5PAR_LIB)

PROGRAM = EMOOPIC
PROGRAM_PAR = $(PROGRAM)_par

all: $(PROGRAM) $(PROGRAM_PAR)

$(PROGRAM): $(OBJ) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) -lm -lz -ldl 

$(PROGRAM_PAR): $(OBJ_PAR) $(LIBS_PAR)
	$(MPICC) -o $@ $^ $(CFLAGS) -lm -lz -ldl $(MPIFLAGS)

%.o: %.cpp $(LIBS)
	$(CC) -c -o $@ $< $(CFLAGS) -I $(GTEST_INCL) -I $(HDF5_INCL)

%_par.o : %.cpp $(LIBS_PAR)
	$(MPICC) -c -o $@ $< $(CFLAGS) $(MPIFLAGS) -I $(GTEST_INCL) -I $(HDF5PAR_INCL)

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
	[ -d $(HDF5_DIR) ] || tar xf hdf5-1.8.12.tar.gz; \
	cd $(HDF5_DIR); \
	CC=cc ./configure --prefix=$(HDF5_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

$(HDF5PAR_LIB):
	cd ..; \
	[ -d $(HDF5_DIR) ] || tar xf hdf5-1.8.12.tar.gz; \
	cd $(HDF5_DIR); \
	CC=mpicc ./configure --enable-parallel --prefix=$(HDF5PAR_DIR); \
	make; \
	make install; \
	cd $(SRCDIR)

clean:
	rm -f $(PROGRAM) $(PROGRAM_PAR) $(OBJ) $(OBJ_PAR) .depend

run_tests:
	cd ../test; \
	make; \
	make run_tests; \
	cd ../src;

