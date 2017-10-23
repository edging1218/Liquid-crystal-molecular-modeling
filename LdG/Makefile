CC=mpicc
CFLAGS=-O3
DEPS = finite.h
OBJ = main.o initial.o scatter.o energy.o relax.o output.o read_param.o conf.o func.o initial_bulk.o initial_drop.o initial_channel.o initial_cylinder.o initial_ellip.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

iso: $(OBJ)
	mpicc -o $@ $^ $(CFLAGS)
