# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/

# librairies de SuiteSparse
L1 = SuiteSparse/UMFPACK/Lib/libumfpack.a
L2 = SuiteSparse/CHOLMOD/Lib/libcholmod.a
L3 = SuiteSparse/AMD/Lib/libamd.a
L4 = SuiteSparse/CAMD/Lib/libcamd.a
L5 = SuiteSparse/COLAMD/Lib/libcolamd.a
L6 = SuiteSparse/CCOLAMD/Lib/libccolamd.a
L7 = SuiteSparse/metis-4.0/libmetis.a
L8 = SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a

# toutes les librairies
LIB = $(LIBP) $(L1) $(L2) $(L3) $(L4) $(L5) $(L6) $(L7) $(L8) -lm -lblas -llapack

COPT = -O3  # Compilation flags

default: main

clean:
	rm *.o
	rm main
	rm src/*.o

main: main.c main.h src/prob.o src/time.o src/funcTools.o src/smoothing.o src/multiGridMethod.o src/geometry_vec.o src/reduction_mat.o src/prolongation_mat.o primme.o umfpack.o
	cc $(COPT) $^ -o $@ $(LIB)

%.o: %.c
	cc $(COPT) -c $< -o $@ $(INCP)

src/%.o: src/%.c
	cc $(COPT) -c $< -o $@ $(INCP)


umfpack.o: umfpack.c
	cc $(COPT) -c $< -o $@ -ISuiteSparse/UMFPACK/Include \
  -ISuiteSparse/SuiteSparse_config  -ISuiteSparse/AMD/Include
