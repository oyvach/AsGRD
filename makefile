# programming environment
COMPILER     := mpicxx 
LFFLD         = /Users/oyvinch/Projects/software/LATfield2

# target and source
EXEC         := AsGRD
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

INCLUDE      :=  -I$(LFFLD) -I/usr/local/include
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas  -L/usr/local/lib

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5
DGEVOLUTION  := 
# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  += -DCYCLE_INFO_INTERVAL=10
DGEVOLUTION  += -DPHINONLINEAR  # non-linear as default
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DHAVE_HEALPIX
DGEVOLUTION  += -DGW
DGEVOLUTION  += -DCIC_PROJECT_TIJ

# further compiler options
OPT          := -O3 -std=c++11 

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $(OPTIM)  $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

lccat: lccat.cpp
	$(COMPILER) $(OPTIM) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE) $(LIB)

lcmap: lcmap.cpp
	$(COMPILER) $(OPTIM) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

clean:
	-rm -f $(EXEC) lccat lcmap force

