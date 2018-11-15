CC = g++
BIN = peax
OBJ = peax.o complex.o fft.o ims_core.o preprocessing.o candidate-detection.o picking.o modeling.o
OPT = -O3 -fopenmp -static

main: $(OBJ)
	$(CC) $(OPT) -o $(BIN) $(OBJ) $(LIBDIR) $(LIBS)

%.o: %.cpp
	$(CC) $(OPT) -o $@ -c $< $(CFLAGS)

new: clean main

clean:
	rm -f *.o $(BIN) *~