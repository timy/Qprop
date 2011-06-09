DIRROOT		= ../../
DIROBJ		= $(DIRROOT)obj/
DIRLIB		= $(DIRROOT)lib/

PACKAGE  	= $(DIRLIB)/libqprop.a

COPT		= -c -O3 -I.
FOPT		= -c -O3 -I.
CC 		= g++
FC		= gfortran

all:	$(PACKAGE)
$(PACKAGE):
	mkdir -p $(DIROBJ)
	mkdir -p $(DIRLIB)
	$(CC) $(COPT) bar.cc -o $(DIROBJ)bar.o
	$(CC) $(COPT) ylm.cc -o $(DIROBJ)ylm.o
	$(CC) $(COPT) wavefunction.cc -o $(DIROBJ)wavefunction.o
	$(CC) $(COPT) hamop.cc -o $(DIROBJ)hamop.o
	$(CC) $(COPT) grid.cc -o $(DIROBJ)grid.o
	$(CC) $(COPT) fluid.cc -o $(DIROBJ)fluid.o
	$(CC) $(COPT) ned.cc -o $(DIROBJ)ned.o
	$(CC) $(COPT) cmatrix.cc -o $(DIROBJ)cmatrix.o
	$(CC) $(COPT) xml_parse.cc -o $(DIROBJ)xml_parse.o
	$(CC) $(COPT) factorial.cc -o $(DIROBJ)factorial.o
	$(CC) $(COPT) winop.cc -o $(DIROBJ)winop.o
	-ar ruvs $(DIRLIB)libqprop.a \
		$(DIROBJ)bar.o \
		$(DIROBJ)ylm.o \
		$(DIROBJ)wavefunction.o \
		$(DIROBJ)hamop.o \
		$(DIROBJ)grid.o \
		$(DIROBJ)fluid.o \
		$(DIROBJ)ned.o \
		$(DIROBJ)cmatrix.o \
		$(DIROBJ)xml_parse.o \
		$(DIROBJ)factorial.o \
		$(DIROBJ)winop.o 
clean:
	-rm -f *~
	-rm -f $(DIROBJ)*.o
	-rm -f $(DIROBJ)*~