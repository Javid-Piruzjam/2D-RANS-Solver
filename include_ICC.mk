CC  = icc
CXX = icpc
FC  = ifort
LINKER = $(FC)

CFLAGS   = -O3 -std=c99 -xHost
CXXFLAGS = $(CFLAGS)
FCFLAGS  = -O3 -xHost
LFLAGS   = 
DEFINES  = -D_GNU_SOURCE
INCLUDES = 
LIBS     = 


