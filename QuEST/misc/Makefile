CC = g++
#CCFLAGS = -Wall -ansi -pedantic -D DEBUG -O4 -funroll-loops

include ../config.mk

TARGETS = seq_contig sequence_utility_functions ticker utils array_2d params
SOURCES = $(TARGETS:=.cpp)
HEADERS = $(TARGETS:=.h)
OBJECTS = $(TARGETS:=.o)

all: $(OBJECTS)

%.o : %.cpp %.h 
	$(CC) $(INCLUDE) $(CCFLAGS) -c $<

clean: 
	rm -rf *.o 

#all: seq_contig.o contig_set.o sequence_utility_functions.o ticker.o utils.o params.o array_2d.o

#%.o: 
#%.o : %.cpp %.h 
#	$(CC) $(INCLUDE) $(CCFLAGS) -c $<
#-o $@.o
#$(OBJ_DIR)/$@.o

#$(TARGETS): $(SOURCES) $(HEADERS)
#	$(CC) $(INCLUDE) $(CCFLAGS) -c $(SOURCES)

#params.o: params.cpp params.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c params.cpp
#
#seq_contig.o: seq_contig.cpp seq_contig.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c seq_contig.cpp
#
#contig_set.o: contig_set.cpp contig_set.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c contig_set.cpp 
#
#sequence_utility_functions.o: sequence_utility_functions.cpp sequence_utility_functions.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c sequence_utility_functions.cpp
#
#ticker.o: ticker.cpp ticker.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c ticker.cpp
#
#utils.o: utils.cpp utils.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c utils.cpp
#
#array_2d.o: array_2d.cpp array_2d.h
#	$(CC) $(INCLUDE) $(CCFLAGS) -c array_2d.cpp
#
#clean: 
#	rm -rf *.o 
