# $Id$
#

CXX = /scr/pigpen/granger/egcs-1.0.3a-5.6/bin/gcc
#CXX = /scr/pigpen/granger/egcs-1.0.1/bin/gcc

CXXDEBUGFLAGS = -g -Wall

INCLUDES = -I./lib
CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES)

.SUFFIXES: .cxx

OBJS = driver.o
LIBS = lib/BSpline.o -lstdc++

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: driver

driver: $(OBJS)
	$(CXX) -o driver $(LDFLAGS) $(OBJS) $(LIBS)

clean::
	rm -f $(OBJS)

depend::

driver.o: lib/BSpline.h 
