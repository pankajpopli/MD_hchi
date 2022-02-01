IDIR=../include
CC=gcc
CFLAGS=-I$(IDIR) -O3 -ffast-math

ODIR=../obj
LDIR =../lib

LIBS=-lm -lgsl -lgslcblas
#-I/usr/include/gsl -L/usr/lib 

_DEPS = global.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main_MD.o cell_list.o generate_latt2.o global.o initializer.o n_list_static.o periodic.o force.o cell_list.o get_quant.o potential_polarForm.o \
		X.o x_honeycomb.o x_kagome.o x_sq.o x_triang.o x_rectangle.o init_Y_matrix.o hchi_force.o hchi_space_potential.o
#dhch.o idchi_square.o dchi_honeycomb.o dchi_kagome.o dchi_triang.o dchi_rectangle.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

rolat.out: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

	
.PHONY: clean install

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
	

install:
	rm -rf ../results/
	mkdir -p ../results
	cp ../lib/p_of_chi_square.c ../results/
	cp ../lib/live_plot.gnu ../results/
	mv rolat.out ../results/
