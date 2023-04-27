CC = mpicc
CCFLAGS = -g -std=c99 -Ofast 
LIBS = -lm 

CFLAGS = $(CCFLAGS)

BINS = main gillespie

.PHONY: clean

all: $(BINS) 

$(BINS): %: %.c prop.o
	$(CC) $(CFLAGS) $< prop.o -o $@ $(LIBS)

prop.o: prop.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(BINS) prop.o
