#/a makefile for modules

include ../make.inc

# location of needed modules
modflags=$(basemod_flags) \
         $(mod_flag)../elpa/src

# list of modules

modules = \
io_files.o \
io_global.o  \
parameters.o \
parser.o \
read_cards.o \
read_input.o \
read_namelists.o \
io_base.o 
#hdf5_qe.o\
#h5_module.o\


# list of subroutines and functions (not modules) 

objs = \
inpfile.o \
plot_io.o \
test_input_file.o 

all : libio08.a

## the following is needed only for lapack compiled from sources

libio08.a: 	$(modules) $(objs)
	$(ar) $(arflags) $@ $?       
	$(ranlib) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.l

include make.depend
