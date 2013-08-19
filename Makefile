# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linux

ObjSuf        = o
SrcSuf        = C
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -g
LD            = g++
LDFLAGS       = 
SOFLAGS       = -shared


ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    := $(shell $(ROOTSYS)/bin/root-config --glibs)

CXXFLAGS     += $(ROOTCFLAGS) 
LIBS          = $(ROOTLIBS) -lHtml -lMinuit 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    

all:            mixingfit

mixingfit:	MixingFit.C MixingFit.h Graphics.o 
		$(LD) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
		@echo "$@ done"

Graphics.o:	Graphics.C MixingFit.h 
		$(CXX) $(CXXFLAGS) -c -o $@ $<
		@echo "$@ done"

clean:
		@rm -f mixingfit $(OBJS) core

.SUFFIXES: .$(SrcSuf)

###


