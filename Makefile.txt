CC=gcc

LDFLAGS=-lpthread 

CFLAGS=-Wall -g -std=gnu99 -o

JacobiC.o: JacobiC.c
	$(CC) JacobiC.c $(LDFLAGS) $(CFLAGS) -c

JacobiC: JacobiC.c
	$(CC) $(CFLAGS) JacobiC JacobiC.c $(LDFLAGS)

Jacobi.java: Jacobi.java
	javac Jacobi.java

JacobiJava:
	javac Jacobi.java

clean: 
	rm -rf JacobiC.o JacobiC JacobiC.log *.class

all:
	$(CC) $(CFLAGS) JacobiC JacobiC.c $(LDFLAGS)
	javac Jacobi.java
