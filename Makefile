CC = gcc
DEPS = ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c mysolver.c
OUT_EXE = mysolver
CFLAGS = -lm -g
build: $(DEPS)
	$(CC) -fPIC $(DEPS) $(CFLAGS) -o $(OUT_EXE) -O3 -ftree-vectorizer-verbose=1
##	$(CC) -fPIC $(DEPS) $(CFLAGS) -o $(OUT_EXE)
clean:
	-rm -f *.o core

rebuild: clean build

.PHONY: clean 
