#! gnumake
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Makefile for qft++ tutorials directory.
#
# Author: Mike Williams
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# declare some variables
#
LD      = g++ 
CC      = gcc
FLAGS   = -O2 -Wall
INCLUDE = -I ./ -I ../include/
LIBS    = ../lib/libqft++.a
SOURCES = $(wildcard *.cxx)
EXES    = $(addprefix bin/,$(basename $(SOURCES)))
DEPENDS = $(addprefix depends/,$(SOURCES:.cxx=.d))
#
# commands
#
.PHONY: all clean
all: exes
#
# general pattern to write depency files
#
depends/%.d: %.cxx
	@echo "Generating dependencies for $*.C ..."
	$(SHELL) -ec '$(CC) -MM $(FLAGS) $(INCLUDE) $< | ../src/depends.pl $@'
	@echo "done."
#
# general pattern to build objects from source files 
#
objects/%.o: %.cxx
	@echo "Building $*.o ..."
	$(LD) $(FLAGS) $(INCLUDE) -c -o objects/$*.o $*.cxx
	@echo "done."
#
# general pattern to build executables from object files
#
bin/%: objects/%.o 
	@echo "Linking $* ..."
	$(LD) $(FLAGS) objects/$*.o $(LIBS) -o bin/$*
	@echo "done."
#
# build the executables
# 
include $(DEPENDS)
exes: $(EXES)
#
# rule for 'cleaning'
#
clean:
	@rm -f bin/*
	@rm -f depends/*.d
	@rm -f objects/*.o
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
