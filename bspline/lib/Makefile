# $Id$
#

CXX = /scr/pigpen/granger/egcs-1.0.3a-5.6/bin/gcc
#CXX = /scr/pigpen/granger/egcs-1.0.1/bin/gcc

CXXDEBUGFLAGS = -g -Wall

CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES)

.SUFFIXES: .cxx

OBJS = BSpline.o

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: $(OBJS)

clean::
	rm -f $(OBJS)

depend::

BSpline.o: BandedMatrix.h BSplineLU.h BSpline.h BSplineSolver.h
