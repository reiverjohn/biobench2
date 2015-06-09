CC = g++
#CCFLAGS = -Wall -ansi -pedantic -D DEBUG -O4 -funroll-loops
OBJECTS := $(patsubst %.cpp, %.o,$(wildcard *.cpp))
MISC_DIR = ./../misc
INCLUDE := $(MISC_DIR)

include ../config.mk

all: $(OBJECTS)

%.o : %.cpp
	$(CC) -I$(INCLUDE) $(CCFLAGS) -c $<

clean: 
	rm -rf *.o 
