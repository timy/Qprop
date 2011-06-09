# user defined
TARGET	=	winop_lin_my
OPTSRC	=	winop_lin_my.cc ../../lib/libqprop.a ../main/util/lib/libutil.a

# other options
CC	=	g++
COPTS	=	
OPTINC	=	-I. -I../core
OPTLIB	=	-llapack


all: libqprop.a $(TARGET)

winop_lin_my : libutil.a libqprop.a
	$(CC) $(COPTS) -o $(TARGET) $(OPTSRC) $(OPTINC) $(OPTLIB)

libqprop.a:
	cd ../core/ && make 

libutil.a:
	cd ../main/util/ && make

clean:
	rm -f ./$(TARGET)
	rm *~

cleanall:
	cd ../core/ && make clean
	rm -f ../../lib/libqprop.a

run:
	qsub run.sh