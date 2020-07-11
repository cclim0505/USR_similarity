PROG =	USR.out

SRCS =  

OBJS =  constants.o initialise.o simil.o simulate.o main.o

### gfortran options
LIBS =
FC = gfortran
FFLAGS = -g -ffree-form -Wall -Wextra -O2 -fbounds-check 
LDFLAGS= 
RUNCPP = 
   
### ifort options   
#FC = ifort
#FFLAGS = -g -O2 -free -ftz -ip -ipo -qopenmp -parallel -prec-div -prec-sqrt

all: $(PROG)

$(PROG): $(OBJS) $(EXTRAS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(EXTRAS) $(LIBS) 

clean:
	rm *.o *.mod

cleaner:
	rm $(PROG)

