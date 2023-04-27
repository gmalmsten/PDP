CC = mpicc
CCFLAGS = -g -std=c99 -Ofast 
LIBS = -lm 

CFLAGS = $(CCFLAGS)

BINS = find0

.PHONY: clean

all: $(BINS) 

$(BINS): %: %.c 
	$(CC) $(CFLAGS) $<  -o $@ $(LIBS)

clean:
	$(RM) $(BINS) prop.o
