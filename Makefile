DOPLOT = 1
FPS    = 1
POINTS = 0

GCC = nvcc 
EXE = entbody
SRC = main.c util.c 
FLAGS = -O3 
LIBFLAGS = -lm

# we want the compile line to be essentially
# nvcc main.cu -arch sm_12 -DPLOT -lGL -lGLU -lglut
FLAGS += -x cu -DCUDA -arch sm_11

ifeq ($(DOPLOT), 1)
    SRC += plot.c
    FLAGS += -DPLOT
    LIBFLAGS += -lGL -lGLU -lglut
endif

ifeq ($(FPS),1)
    LIBFLAGS += -lrt
    FLAGS += -DFPS
endif

ifeq ($(POINTS), 1)
    FLAGS += -DPOINTS
endif

# default super-target
all: $(EXE)

# the standard executable
$(EXE): $(SRC) 
	$(GCC) $(FLAGS) $^ -o $@ $(LIBFLAGS)

clean: 
	rm -rf $(EXE)

.PHONY: clean
