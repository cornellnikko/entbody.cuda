DOPLOT = 0
DOFILE = 1

GCC=nvcc
EXE=entbody
SRC=main.c util.c
FLAGS=-O3
LIBFLAGS=-lm -lrt

# we want the compile line to be essentially
# nvcc main.cu -arch sm_12 -DPLOT -lGL -lGLU -lglut
FLAGS += -x cu -DCUDA  #-arch=sm_13

ifeq ($(DOPLOT), 1)
    SRC += plot.c
    FLAGS += -DPLOT
    LIBFLAGS += -lGL -lGLU -lglut
endif

ifeq ($(DOFILE), 1)
	SRC += io.c
endif

# default super-target
all: $(EXE)

# the standard executable
$(EXE): $(SRC)
	$(GCC) $(FLAGS) $^ -o $@ $(LIBFLAGS)

clean:
	rm -rf $(EXE)
	rm -rf particles.csv
.PHONY: clean
