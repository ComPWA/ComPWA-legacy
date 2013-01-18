#! gnumake
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# Global Makefile for qft++ src directory.
#
# Author: Mike Williams
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
# declare some variables
#
LD           = g++
SUBDIRS      = tensor relativistic-quantum-mechanics
MAKE_OPTIONS = --no-print-directory
#
# targets
#
.PHONY: all ${SUBDIRS} clean
all: ${SUBDIRS} lib
#
# rule for building a subdirectory
#
$(SUBDIRS):
	@echo "Making $@ directory ..."
	@$(MAKE) -f ../Module.mk -C $@ $(MAKE_OPTIONS)
	@echo "done"
#
# rule for linking the library
#
REL_QM_SRCS = $(wildcard relativistic-quantum-mechanics/*.C)
OBJS        = $(subst mechanics,mechanics/objects,$(REL_QM_SRCS:.C=.o))
TENSOR_SRCS = $(wildcard tensor/*.C)
OBJS       += $(subst tensor,tensor/objects,$(TENSOR_SRCS:.C=.o))
lib: ../lib/libqft++.so ../lib/libqft++.a
../lib/libqft++.so: ${OBJS}
	@echo "Linking libqft++.so ..."
	${LD} -shared ${OBJS} -o ../lib/libqft++.so
	@chmod 555 ../lib/libqft++.so
	@echo "done."
../lib/libqft++.a: $(OBJS)
	@echo "Archiving libqft++.a ..."
	@ar rv ../lib/libqft++.a $(OBJS)
	@echo "done."
#
# rule for 'cleaning'
#
clean:
	@rm -f tensor/objects/*.o
	@rm -f relativistic-quantum-mechanics/objects/*.o
	@rm -f ../lib/*.a ../lib/*.so
	@rm -f tensor/depends/*.d
	@rm -f relativistic-quantum-mechanics/depends/*.d
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
