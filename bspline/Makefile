# $Id$
#

CXX = g++

CXXDEBUGFLAGS = -g -O -Wall

# Include comparison with Ooyama FORTRAN implementation on sun
# DEFINES = -DOOYAMA
# FLIBS = -L/opt/SUNWspro/SC4.0/lib -lF77 -lV77 -lM77 -lsunmath -lm
# OOYAMA_OBJS = ooyama/ooyama.o

INCLUDES = -I./lib
CXXFLAGS = $(CXXDEBUGFLAGS) $(INCLUDES) $(DEFINES)
LDFLAGS = -pg
FFLAGS = -O

.SUFFIXES: .cxx

OBJS = driver.o $(OOYAMA_OBJS)
LIBS =  -lstdc++ $(FLIBS)

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: driver

# SUBDIRS = subdirs

driver: $(SUBDIRS) $(OBJS) 
	$(CXX) -o driver $(LDFLAGS) $(OBJS) $(LIBS)

subdirs:
	@set -x; cd lib; \
	$(MAKE) $(MFLAGS) LDFLAGS='$(LDFLAGS)' \
		CXXDEBUGFLAGS='$(CXXDEBUGFLAGS)' all

TARFILE = bspline
TARDEST = .
ZIPPIPE = gzip -c
ZIPEXT = gz
LN = ln -s
DISTFILES = \
     README COPYRIGHT Makefile driver.cxx bspline.dsw driver.dsp \
     lib/BSpline.cxx lib/BSpline.h lib/BSplineD.cxx lib/BSplineF.cxx \
     lib/BandedMatrix.h lib/Makefile lib/BSpline.dsp lib/BSpline.dsw doc

DISTSUBDIRS = 

tarfile: $(DISTFILES)
	@echo Making file $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	if $(LN) . $(TARFILE) ; then \
	  tar cvhf - `ls $(DISTFILES) | sed -e 's,^,$(TARFILE)/,'` | \
		  $(ZIPPIPE) > $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	  rm -f $(TARFILE) ; \
	  echo Done making $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	else \
	  echo "Cannot make link '$(TARFILE)'!  Aborting."; \
	  false ; \
	fi

clean::
	rm -f $(OBJS)

clean::
	cd lib; $(MAKE) $(MFLAGS) clean

depend::

driver.o: lib/BSpline.cxx lib/BSpline.h 
