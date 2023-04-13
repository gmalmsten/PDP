############################
# Makefile for lab 1
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -g -std=c99 -O3
LIBS = -lm

BINS = conti array_sect

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

