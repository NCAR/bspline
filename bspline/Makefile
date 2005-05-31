# $Id$
#

PUBDIR=/net/www/docs/homes/granger/bspline

CXX = g++

CXXDEBUGFLAGS = -g -O -Wall

# Include comparison with Ooyama FORTRAN implementation on sun
# DEFINES = -DOOYAMA
# FLIBS = -L/opt/SUNWspro/SC4.0/lib -lF77 -lV77 -lM77 -lsunmath -lm
# OOYAMA_OBJS = ooyama/ooyama.o

INCLUDES = -I./lib -I./BSpline
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

VERSION = dev
TARFILE = bspline-$(VERSION)
TARDEST = .
ZIPPIPE = gzip -c
ZIPEXT = gz
LN = ln -s
DISTFILES = \
     README COPYRIGHT Makefile driver.cpp bspline.dsw driver.dsp \
     lib/BSpline.cpp BSpline/BSpline.h lib/BSplineD.cpp lib/BSplineF.cpp \
     lib/BandedMatrix.h lib/Makefile lib/BSpline.dsp lib/BSpline.dsw \
     sample.temps sample.txt sample.wdir sample.wspd plots.R \
     lib/Doxyfile doc

DISTSUBDIRS = 

tarfile: $(DISTFILES)
	@echo Making file $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	if $(LN) . $(TARFILE) ; then \
	  tar cvhf - `ls -d $(DISTFILES) | sed -e 's,^,$(TARFILE)/,'` | \
		  $(ZIPPIPE) > $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	  echo Done making $(TARDEST)/$(TARFILE).tar.$(ZIPEXT) ; \
	  echo Making $(TARFILE).zip ; \
	  zip -r $(TARDEST)/$(TARFILE).zip \
	    `ls -d $(DISTFILES) | sed -e 's,^,$(TARFILE)/,'` ; \
	  rm -f $(TARFILE) ; \
	else \
	  echo "Cannot make link '$(TARFILE)'!  Aborting."; \
	  false ; \
	fi

IMAGES = plot-1.png plot-2.png plot-3.png plot-4.png plot-5.png plot-6.png

doc::
	cd lib; doxygen
	test -d doc && cp -pr $(IMAGES) doc

clean::
	cd lib; $(MAKE) $(MFLAGS) clean

depend::

driver.o: lib/BSpline.cpp BSpline/BSpline.h 

tag::
	@if [ x$(VERSION) = xdev ]; \
	then echo "Need a VERSION string."; exit 1 ; \
	else cvs tag -l $(CVSOPTS) Version_$(VERSION) . lib BSpline ooyama ;\
	fi

release::
	@if [ x$(VERSION) = xdev ]; then \
	echo "Need a VERSION string."; exit 1 ; else \
	rm -rf $(TARFILE) ; \
	cvs co $(CVSOPTS) -r Version_$(VERSION) -d $(TARFILE) bspline;\
	(cd $(TARFILE) && $(MAKE) VERSION=$(VERSION) TARDEST=.. tarfile);\
	fi

publish:: release
	mv -f $(TARFILE).tar.gz $(PUBDIR)

