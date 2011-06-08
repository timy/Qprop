# user defined target name
TARGET	=	h1s
# the source file needed
OPTSRC	=	parameters.cc toolbox.cc $(TARGET).cc ../../lib/libqprop.a \
		../../lib/libutil.a

# other options
CC	=	g++
FC	=	gfortran
COPTS	=	-O3 -Wall
OPTINC	=	-I. -I../core -I../util
OPTLIB	=	-lgfortran -llapack

all: libqprop.a $(TARGET)

$(TARGET) : libqprop.a libutil.a $(TARGET).cc
	$(CC) $(COPTS) -o $(TARGET) $(OPTSRC) $(OPTINC) $(OPTLIB)

libqprop.a:
	cd ../core/ && make 

libutil.a:
	cd ../util/ && make
clean:
	rm -f ./$(TARGET)
	rm -f *~
	rm -f core
	rm -f screen.log

cleanall: clean
	cd ../core/ && make clean
	cd ../util/ && make clean
	rm -f ../../lib/libqprop.a
	rm -f ../../lib/libutil.a

cleandat:
	cd ./res && rm -f *.dat

run:
	qsub6 run.sh