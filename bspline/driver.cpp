// -*- mode: c++; c-basic-offset: 4; -*-
//
// $Id$
//
// Simple test driver for the bspline interface.
//

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

using namespace std;

#include "BSpline.h"

static void 
DumpSpline (vector<float> &x, vector<float> &y, BSpline &spline, ostream &out);
static void
EvalSpline (BSpline &spline, ostream &out);
#if OOYAMA
static bool
vic (float *xt, int nxp, float wl, int bc, float *y);
extern "C" 
{
    void vicsetup_ (float *xt, int *nxp, float *ydcwl, int *nx, 
		    float *ynb, float *ynt, float *fmin, int *ierr, int *echo);
    void vicspl_ (float *xt, float *xd, float *y, int *nxp, int *nxp, 
		  int *kdat, float *ynb, float *ynt, int *nx,
		  float *ydcwl, int *kybc, int *kybc, 
		  float *ydcwl, float *ydcwl, int *ierr);
    void spotval_ (float *xi, int *kdat, float *fout, float *foutd);
}
#endif /* OOYAMA */


int
main (int argc, char *argv[])
{
    // Some basics first
    cout << "BSpline interface version: "
	 << BSplineBase::IfaceVersion() << endl;
    cout << "BSpline implementation version: "
	 << BSplineBase::ImplVersion() << endl;
	
	// Sub-sampling and cutoff wavelength come from command-line
    if (argc != 3)
    {
	cerr << "Usage: " << argv[0] << " <step> <cutoff>" << endl;
	exit (1);
    }
	
    int step = (int) atof (argv[1]);
    float wl = atof (argv[2]);

    // Read the x and y pairs from stdin
    vector<float> x;
    vector<float> y;
    float f, g;
    float base = 0;
    int i = 0;

    while (cin >> f)
    {
	if (++i == 1)
	{
	    base = f;
	}
	f -= base;

	if (cin >> g)
	{
	    x.push_back (f);
	    y.push_back (g);
	    // cout << x.back() << " " << y.back() << endl;
	}
	else
	    break;
    }

    // Subsample the arrays
    vector<float>::iterator xi = x.begin(), yi = y.begin();
    vector<float>::iterator xo = x.begin(), yo = y.begin();
    if (step > 1)
    {
	while (xo < x.end() && yo < y.end())
	{
	    //cout << *xo;
	    //cout << " ";
	    //cout << *yo;
	    //cout << endl;
	    *xi++ = *xo;
	    *yi++ = *yo;
	    xo += step;
	    yo += step;
	}
	x.resize (xi - x.begin());
	y.resize (yi - y.begin());
    }
		
    // Create our bspline base on the X vector with a simple 
    // wavelength.
    int bc = BSplineBase::BC_ZERO_SECOND;
    BSpline::Debug = true;
#if 0
    //float wl = (x.back() - x.front()) / (x.size()/3 + 1);
    //float wl = 30.0; /*secs*/
    BSplineBase::Debug = true;
    BSplineBase bb (x.begin(), x.size(), wl, bc);
    if (bb.ok())
    {
	// Now apply our y values to get the smoothed curve.
	BSpline spline(bb, y.begin());
    }
#endif
    BSpline spline (x.begin(), x.size(), y.begin(), wl, bc);
    if (spline.ok())
    {
	// And finally write the curve to a file
	ofstream fspline("spline.out");
	ofstream fcurve("smooth.out");
	DumpSpline (x, y, spline, fspline);
	EvalSpline (spline, fcurve);
    }
    else
	cerr << "Spline setup failed." << endl;

#if OOYAMA
    cerr << "Computing Ooyama FORTRAN results." << endl;
    if (vic (x.begin(), x.size(), wl, bc, y.begin()))
    {
	cerr << "Done." << endl;
	ofstream vspline("ooyama.out");
	ostream_iterator<float> of(vspline, "\t ");

	float fout, foutd;
	int kdat = 1;
	float xi = x.front();
	float xs = (x.back() - x.front())/2000.0;
	for (int i = 0; i < 2000; ++i, xi += xs)
	{
	    spotval_ (&xi, &kdat, &fout, &foutd);
	    *of++ = xi;
	    *of++ = fout;
	    vspline << endl;
	}
    }
    else
    {
	cerr << "vic() failed." << endl;
    }
#endif /* OOYAMA */

    return 0;
}


void
DumpSpline (vector<float> &x, vector<float> &y, BSpline &spline, ostream &out)
{
#ifdef notdef
    int nx;
    const float *x = spline.nodes (&nx);
    const float *y2 = spline.curve ();
    const float *end = x + nx;
#endif
    ostream_iterator<float> of(out, "\t ");

    for (unsigned int i = 0; i < x.size(); ++i)
    {
	*of++ = x[i];
	*of++ = y[i];
	*of++ = spline.evaluate (x[i]);
	out << endl;
    }
}



void
EvalSpline (BSpline &spline, ostream &out)
{
    ostream_iterator<float> of(out, "\t ");

    float x = spline.Xmin();
    float xs = (spline.Xmax() - x) / 2000.0;
    for (; x <= spline.Xmax(); x += xs)
    {
	*of++ = x;
	*of++ = spline.evaluate (x);
	out << endl;
    }
}


/*
 * This is the FORTRAN code which computes the spline and evaluates it.
 */
#if 0
      CALL VICSETUP(XT,NDIMX,YDCWL,NX,YNB,YNT,2.0,IERR,ECHO)
      IF (IERR.NE.0) GOTO 900
      CALL VICSPL(XT,XW,X,NDIMX,MXRC,KDAT,YNB,YNT,NX,YDCWL,
     *            KYBC,KYBC,YBCWL,YBCWL,IERR)
      IF (IERR.NE.0) GOTO 900
      DO 500 I = 3,NDIMX-2
         IF (X(I).LE.-999.) GOTO 500
	 CALL SPOTVAL(XT(I),KDAT,FOUT,FOUTD)
         IF (ABS(X(I)-FOUT).GT.DEVX) THEN 
              NBAD = NBAD+1 
              IBAD(NBAD) = I
              NBADX = NBADX+1 
              ENDIF 
500      CONTINUE 
#endif


#if OOYAMA
static bool
vic (float *xt, int nxp, float wl, int bc, float *y)
{
    static float fmin = 2.0;
    int ierr;
    int nx;
    float ynb, ynt;
    int echo = 1;
    float ydcwl = wl;
    int kybc = bc;			// endpoint boundary conditions
    int kdat = 1;
    vector<float> xw(nxp, 1.0);		// relative weights all set to 1

    ierr = 1;
    vicsetup_ (xt, &nxp, &ydcwl, &nx, &ynb, &ynt, &fmin, &ierr, &echo);
    if (ierr)
	return (false);

    ierr = 1;
    vicspl_ (xt, xw.begin(), /*xdat*/ y, &nxp, &nxp, &kdat,
	     &ynb, &ynt, &nx,
	     &ydcwl, &kybc, &kybc, &ydcwl, &ydcwl, &ierr);
    if (ierr)
	return (false);

    return true;
}

#endif /* OOYAMA */
