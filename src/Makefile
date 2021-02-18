
FC 			= ifort
FCFLAGS 	= -O2

LIBS 		= -L/home/neelravi/devel/libfdf-original/lib
INCLUDES 	= -I/home/neelravi/devel/libfdf-original/include

LDFLAGS 	= -lfdf

SOURCE 		= iochamp.f90
OBJECTS		= iochamp.o

#default: parser


#parser: $(OBJECTS)
#	$(FC) -o parser $(OBJECTS) 

#$(OBJECTS): $(SOURCE)
$(FC) $(FCFLAGS) $(SOURCE) $(LDFLAGS) $(LIBS) $(INCLUDES) -o $(OBJECTS) 

#clean:
#	$(RM) $(OBJECTS)

