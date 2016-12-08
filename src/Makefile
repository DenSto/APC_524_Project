USE_MPI = yes

ifeq ($(USE_MPI),no)
  CC = g++
  OMPFLAGS = -DUSE_MPI=0
else
  CC = mpic++
  #OMPFLAGS = -fopenmp -DUSE_MPI=1
  OMPFLAGS = -DUSE_MPI=1
endif

CFLAGS= -g -lm -Wall -stdlib=libstdc++  
PROGRAM = EMOOPIC

IOOBJ = IO.o
DOMOBJ = domain.o
PUSHOBJ = boris.o
GRIDOBJ = grid.o oGrid.o
PARTOBJ = particle.o particle_handler.o interpolate.o deposit.o

OBJ = $(addprefix ./IO/, $(IOOBJ)) \
      $(addprefix ./domain/, $(DOMOBJ)) \
      $(addprefix ./pusher/, $(PUSHOBJ)) \
      $(addprefix ./grid/, $(GRIDOBJ)) \
      $(addprefix ./particles/, $(PARTOBJ)) \
      main.o 

all: $(PROGRAM)

$(PROGRAM): $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAGS)

%.o: %.cpp 
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAGS)

clean:
	rm -f $(PROGRAM) $(OBJ) 
