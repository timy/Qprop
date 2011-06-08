DIROBJ		= ../../obj/
DIRLIB		= ../../lib/

PACKAGE  	= $(DIRLIB)libutil.a

COPT		= -c -O3 -Wall -I. -I./math -I../core
FOPT		= -c -O3 -I.
CC 		= g++
FC		= gfortran

all:	$(PACKAGE)
$(PACKAGE):

	$(FC) $(FOPT) ./math/adgm.f -o $(DIROBJ)adgm.o
	$(FC) $(FOPT) ./math/msphi.for -o $(DIROBJ)msphi.o
	$(CC) $(COPT) ./math/sphi.cc -o $(DIROBJ)sphi.o
	$(CC) $(COPT) vp_type.cc -o $(DIROBJ)vp_type.o
	$(CC) $(COPT) sp_type.cc -o $(DIROBJ)sp_type.o
	$(CC) $(COPT) ip_type.cc -o $(DIROBJ)ip_type.o
	$(CC) $(COPT) fd_type.cc -o $(DIROBJ)fd_type.o
	$(CC) $(COPT) wf_type.cc -o $(DIROBJ)wf_type.o
	-ar ruvs $(DIRLIB)libutil.a \
		$(DIROBJ)adgm.o \
		$(DIROBJ)msphi.o \
		$(DIROBJ)sphi.o \
		$(DIROBJ)vp_type.o \
		$(DIROBJ)sp_type.o \
		$(DIROBJ)ip_type.o \
		$(DIROBJ)fd_type.o \
		$(DIROBJ)wf_type.o

clean:
	-rm -f *~
	-rm -f $(DIROBJ)*.o
	-rm -f $(DIROBJ)*~
	-rm -f $(PACKAGE)