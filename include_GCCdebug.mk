CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

ANSI_CFLAGS  = -ansi
ANSI_CFLAGS += -std=c++0x
ANSI_CFLAGS += -pedantic
ANSI_CFLAGS += -Wextra

CFLAGS   = -g -O3 -Wno-format  -Wall $(ANSI_CFLAGS)
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
CPPFLAGS = -std=c++0x
LFLAGS   =  
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
