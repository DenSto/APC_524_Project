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

IOOBJ = IO.o
DOMOBJ = domain.o ghosts.o pass_fields.o pass_parts.o
PUSHOBJ = boris.o
B_PARTICLEOBJ = bc_p_MPI.o bc_p_periodic.o bc_p_reflecting.o
BOUNDARYOBJ = boundary_particles.o \
			  $(addprefix b_particles/, $(B_PARTICLEOBJ))
GRIDOBJ = grid.o oGrid.o spookyGrid.o
PARTOBJ = particle.o particle_handler.o interpolate.o deposit.o

OBJ = $(addprefix ./IO/, $(IOOBJ)) \
      $(addprefix ./domain/, $(DOMOBJ)) \
      $(addprefix ./pusher/, $(PUSHOBJ)) \
      $(addprefix ./boundaries/, $(BOUNDARYOBJ)) \
      $(addprefix ./grid/, $(GRIDOBJ)) \
      $(addprefix ./particles/, $(PARTOBJ)) \
      main.o 

SRCDIR = $(PWD)
LCONFIG_DIR = ../libconfig-1.5
LCONFIG = $(LCONFIG_DIR)/lib/libconfig++.a

LIBS = $(LCONFIG)

all: $(PROGRAM)

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

-include .depend

