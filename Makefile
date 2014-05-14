CC = gcc
DEPS = ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c mysolver.c
OUT_EXE = mysolver
CFLAGS = -lm
build: $(DEPS)
	$(CC) -fPIC $(DEPS) $(CFLAGS) -o $(OUT_EXE) -Ofast -ftree-vectorizer-verbose=1 -pg
debug: $(DEPS)
	$(CC) -Wall -fPIC $(DEPS) $(CFLAGS) -o -g $(OUT_EXE) 
clean:
	-rm -f *.o core

rebuild: clean build

.PHONY: clean 
