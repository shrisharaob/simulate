CC = gcc
DEPS = ran2.c gasdev.c nrutil.c rk4.c auxFunctions.c rkdumb.c mysolver.c
OUT_EXE = mysolver
CFLAGS = -lm
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
