#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Rules for compiling qft++ targets.
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# define some globals
#
LD    = g++ 
CC    = gcc
FLAGS = -O2 -fPIC -Wall 
#
# commands
#
all: objs
#
# general pattern to write depency files
#
depends/%.d: %.C
	@echo "Generating dependencies for $*.C ..."
	$(SHELL) -ec '$(CC) -MM $(FLAGS) $(INCLUDE) $< | ../depends.pl $@'
	@echo "done."
#
# general pattern to build objects from source files 
#
objects/%.o: %.C
	@echo "Building $*.o ..."
	$(LD) $(FLAGS) $(INCLUDE) -c -o objects/$*.o $*.C 
	@echo "done."
#
# local variables
#
INCLUDE = -I ./
SOURCES = $(wildcard *.C)
OBJS    = $(addprefix objects/,$(SOURCES:.C=.o))
DEPENDS = $(addprefix depends/,$(SOURCES:.C=.d))
include $(DEPENDS)
objs: $(OBJS)
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
