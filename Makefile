CC = gcc
DEPS = ran1.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb_proper.c mysolver.c
OUT_EXE = mysolver.out
CFLAGS = -lm -lgsl -lgslcblas
build: $(DEPS)
	$(CC) -fPIC $(DEPS) $(CFLAGS) -o $(OUT_EXE) -Ofast -ftree-vectorizer-verbose=1 
debug: $(DEPS)
	$(CC) $(DEPS) $(CFLAGS) -o $(OUT_EXE) -g
buildProf: $(DEPS)
	$(CC) -fPIC $(DEPS) $(CFLAGS) -o $(OUT_EXE) -Ofast -ftree-vectorizer-verbose=1  -pg
clean:
	-rm -f *.o core

rebuild: clean build

.PHONY: clean 
