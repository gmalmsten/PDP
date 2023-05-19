CC = mpicc
CCFLAGS = -g -std=c99 -Ofast -fopt-info-vec
LIBS = -lm 

CFLAGS = $(CCFLAGS)

BINS = mc

.PHONY: clean

all: $(BINS) 

$(BINS): %: %.c 
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	$(RM) $(BINS) prop.o
