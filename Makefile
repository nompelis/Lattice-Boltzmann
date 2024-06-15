 CC = gcc
 COPTS = -g -O0 -Wall -fPIC

 CXX = g++
 CXXOPTS = -g -O0 -Wall -fPIC

 DEBUG = -D  _DEBUG_  -D  _DEBUG2_  -D  _DEBUG3_

 LIBS = -lm

##### targets #####
all:
	$(CC) -c $(COPTS) $(DEBUG) mesh.c
	$(CC) -c $(COPTS) $(DEBUG) lb.c
	$(CC)    $(COPTS) $(DEBUG) main.c \
                 mesh.o lb.o $(LIBS)

clean:
	rm -f a.out *.o

