# $Id$
#

CXX = /scr/pigpen/granger/egcs-1.0.3a-5.6/bin/gcc
#CXX = /scr/pigpen/granger/egcs-1.0.1/bin/gcc

CXXDEBUGFLAGS = -g -O -Wall -pg

# Include comparison with Ooyama FORTRAN implementation on sun
DEFINES = -DOOYAMA
INCLUDES = -I./lib
CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES) $(DEFINES)
LDFLAGS = -pg
FFLAGS = -O

.SUFFIXES: .cxx

FLIBS = -L/opt/SUNWspro/SC4.0/lib -lF77 -lV77 -lM77 -lsunmath -lm
BSPLINE = lib/BSpline.o
OBJS = $(BSPLINE) driver.o ooyama.o
LIBS =  -lstdc++ $(FLIBS)

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: driver

driver: subdirs $(OBJS) 
	$(CXX) -o driver $(LDFLAGS) $(OBJS) $(LIBS)

subdirs:
	@set -x; cd lib; \
	$(MAKE) $(MFLAGS) LDFLAGS='$(LDFLAGS)' \
		CXXDEBUGFLAGS='$(CXXDEBUGFLAGS)' all

clean::
	rm -f $(OBJS)

clean::
	cd lib; $(MAKE) $(MFLAGS) clean

depend::

driver.o: lib/BSpline.h 
