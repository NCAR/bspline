# $Id$
#

CXX = /scr/pigpen/granger/egcs-1.0.3a-5.6/bin/gcc
#CXX = /scr/pigpen/granger/egcs-1.0.1/bin/gcc

CXXDEBUGFLAGS = -g -Wall

CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES)

.SUFFIXES: .cxx

DOCXX = doc++
DOCXX_HEADERS = BSpline.h
DOCXX_HTMLDIR = /net/global/www/homes/granger/bspline

OBJS = BSpline.o

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: $(OBJS)

html::
	doc++ -d $(DOCXX_HTMLDIR) $(DOCXX_HEADERS)

clean::
	rm -f $(OBJS)

depend::

BSpline.o: BandedMatrix.h BSpline.h
