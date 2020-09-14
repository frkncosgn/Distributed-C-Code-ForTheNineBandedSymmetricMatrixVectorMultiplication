CC=mpicc

CFLAGS=-c -Wall

NAME=-o mainP.exe

NoC=-np 4

CRUN=mpirun

all: sparseProgram

sparseProgram: mainP.exe
	$(CRUN) $(NoC) mainP.exe

mainP.exe: mainP.o
	$(CC) $(NAME) mainP.o -lm 

mainP.o:	mainP.c
	$(CC) $(CFLAGS) mainP.c

clean:
	rm -rf *o mainP.exe

