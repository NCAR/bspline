# $Id$
#

CXX = g++

CXXDEBUGFLAGS = -g -Wall

CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES)

.SUFFIXES: .cxx

DOCXX = doc++
DOCXX_HEADERS = BSpline.h
DOCXX_HTMLDIR = ../doc

OBJS = BSplineF.o BSplineD.o

LIB = libbspline.a

.cxx.o:
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

$(OBJS): BSpline.cxx BandedMatrix.h BSpline.h
