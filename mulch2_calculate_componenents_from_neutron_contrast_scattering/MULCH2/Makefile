#+++++++++++++++++++++++++++++++
# MULCh program suite Makefile
#++++++++++++++++++++++++++++++++

SHELL = /bin/sh
TARGETContrast = Contrast2
TARGETRg = Rg2
TARGETCompost = Compost2
HOSTNAMESHORT := $(shell hostname)

#++++++++++++++++++++++++++++++++
# IMPORTANT: uncomment one of the 
#            following lines only
#++++++++++++++++++++++++++++++++

#MYTYPE = other
MYTYPE = gcc

#++++++++++++++++++++++++++++++++
# Settings for gcc (gnu compilers)
#++++++++++++++++++++++++++++++++
ifeq ($(MYTYPE),gcc)
TYPE = "GCC (GNU compilers)"
CC = cc
CFLAGS = -O3
#EXECEXT = .exe
# Uncomment EXECEXT if using cygwin
endif
#++++++++++++++++++++++++++++++++

CSRC = mathFunctions2.c contrast2.c atom2.c molecule2.c compost2.c Rg2.c

COBJS := $(CSRC:.c=.o)
LIB = -lm

.c.o:
	$(CC) -o $*.o -c $(CFLAGS) $*.c

$(TARGETContrast): print $(COBJS)
	$(CC) $(CFLAGS) -o $(TARGETContrast) mathFunctions2.o contrast2.o atom2.o molecule2.o $(LIB)
	@echo "" 
	@echo "This version of" $@ "was built for" $(TYPE) "on" $(HOSTNAMESHORT)
$(TARGETRg): print $(COBJS)
	$(CC) $(CFLAGS) -o $(TARGETRg) mathFunctions2.o Rg2.o $(LIB)
	@echo "" 
	@echo "This version of" $@ "was built for" $(TYPE) "on" $(HOSTNAMESHORT)
$(TARGETCompost): print $(COBJS)
	$(CC) $(CFLAGS) -o $(TARGETCompost) mathFunctions2.o compost2.o $(LIB)
	@echo "" 
	@echo "This version of" $@ "was built for" $(TYPE) "on" $(HOSTNAMESHORT)
$(TARGETCompost): print $(COBJS)
	$(CC) $(CFLAGS) -o $(TARGETCompost) mathFunctions2.o compost2.o $(LIB)
	@echo "" 
	@echo "This version of" $@ "was built for" $(TYPE) "on" $(HOSTNAMESHORT)
print:  
	@echo ""
	@echo "Building" $(TARGET) "for the following system:" $(TYPE) "on" $(HOSTNAMESHORT)
	@echo ""
clean :
	-@rm -f $(COBJS)
distclean : clean
	-@rm -f $(TARGETContrast)$(EXECEXT) $(TARGETRg)$(EXECEXT) $(TARGETCompost)$(EXECEXT)

