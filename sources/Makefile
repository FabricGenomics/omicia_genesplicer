CC=g++

all: genesplicer_omicia genesplicer

genesplicer_omicia: genesplicer_omicia.o sites.o
	test -e ../bin || mkdir ../bin
	$(CC) -o ../bin/genesplicer_omicia genesplicer_omicia.o sites.o -lm

genesplicer: genesplicer.o sites.o
	$(CC) -o genesplicer genesplicer.o sites.o -lm

genesplicer_omicia.o: genesplicer_omicia.cpp
	$(CC) -c genesplicer_omicia.cpp

genesplicer.o: genesplicer.cpp
	$(CC) -c genesplicer.cpp

sites.o: sites.c
	$(CC) -c sites.c

clean:
	 $(RM) *.o *~ $(ALL)

