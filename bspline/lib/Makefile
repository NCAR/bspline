# $Id$
#

INCLUDES = -I../BSpline

CXX = g++

CXXDEBUGFLAGS = -g -Wall

CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES)

.SUFFIXES: .cpp

DOCXX = doc++
DOCXX_HEADERS = ../BSpline/BSpline.h
DOCXX_HTMLDIR = ../doc

OBJS = BSplineF.o BSplineD.o

LIB = libbspline.a

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: $(LIB)

$(LIB): $(OBJS)
	rm -f $(LIB)
	ar crv $(LIB) $(OBJS)
	ranlib $(LIB)

html::
	doc++ -d $(DOCXX_HTMLDIR) $(DOCXX_HEADERS)

clean::
	rm -f $(OBJS) $(LIB) core *~

depend::

$(OBJS): BSpline.cpp BandedMatrix.h ../BSpline/BSpline.h
