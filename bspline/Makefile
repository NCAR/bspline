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
# LDFLAGS = -pg
FFLAGS = -O

.SUFFIXES: .cpp

OBJS = driver.o $(OOYAMA_OBJS)
LIBS =  -lstdc++ $(FLIBS)

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

all:: driver

# SUBDIRS = subdirs

driver: $(SUBDIRS) $(OBJS) 
	$(CXX) -o driver $(LDFLAGS) $(OBJS) $(LIBS)

clean::
	rm -f driver $(OBJS) *~

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
     README COPYRIGHT Makefile driver.cpp bspline.dsw driver.dsp \
     lib/BSpline.cpp lib/BSpline.h lib/BSplineD.cpp lib/BSplineF.cpp \
     lib/BandedMatrix.h lib/Makefile lib/BSpline.dsp lib/BSpline.dsw doc

DISTSUBDIRS = 

tarfile: $(DISTFILES)
	@echo Making file $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	if $(LN) . $(TARFILE) ; then \
	  tar cvhf - `ls -d $(DISTFILES) | sed -e 's,^,$(TARFILE)/,'` | \
		  $(ZIPPIPE) > $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	  echo Done making $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	  echo Making $(TARFILE).zip ; \
	  zip -r $(TARFILE).zip \
	    `ls -d $(DISTFILES) | sed -e 's,^,$(TARFILE)/,'` ; \
	  rm -f $(TARFILE) ; \
	else \
	  echo "Cannot make link '$(TARFILE)'!  Aborting."; \
	  false ; \
	fi

clean::
	cd lib; $(MAKE) $(MFLAGS) clean

depend::

driver.o: lib/BSpline.cpp lib/BSpline.h 
