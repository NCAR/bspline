// -*- mode: c++; c-basic-offset: 4; -*-
//
// $Id$
//
// Simple test driver for the bspline interface.
//

#include "BSpline.cxx"

template class BSpline<double>;
template class BSpline<float>;

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

using namespace std;

//typedef float datum;
typedef double datum;
typedef BSpline<datum> SplineT;
typedef BSplineBase<datum> SplineBase;
typedef BSpline<double> SplineD;

void 
DumpSpline (vector<datum> &x, vector<datum> &y, SplineT &spline, ostream &out);
static void
EvalSpline (SplineT &spline, ostream &out);
#if OOYAMA
static bool
vic (float *xt, int nxp, double wl, int bc, float *y);
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
	 << SplineBase::IfaceVersion() << endl;
    cout << "BSpline implementation version: "
	 << SplineBase::ImplVersion() << endl;
	
	// Sub-sampling and cutoff wavelength come from command-line
    if (argc != 3)
    {
	cerr << "Usage: " << argv[0] << " <step> <cutoff>" << endl;
	exit (1);
    }
	
    int step = (int) atof (argv[1]);
    double wl = atof (argv[2]);

    // Read the x and y pairs from stdin
    vector<datum> x;
    vector<datum> y;
    datum f, g;
    datum base = 0;
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
    vector<datum>::iterator xi = x.begin(), yi = y.begin();
    vector<datum>::iterator xo = x.begin(), yo = y.begin();
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
    int bc = SplineBase::BC_ZERO_SECOND;
    SplineT::Debug = true;
    SplineT spline (x.begin(), x.size(), y.begin(), wl, bc);
    if (spline.ok())
    {
	// And finally write the curve to a file
	ofstream fspline("input.out");
	DumpSpline (x, y, spline, fspline);
	ofstream fcurve("spline.out");
	EvalSpline (spline, fcurve);
    }
    else
	cerr << "Spline setup failed." << endl;

#if OOYAMA
    // Need to copy the data into float arrays for vic().
    vector<float> fx;
    vector<float> fy;
    fx.resize (x.size());
    fy.resize (y.size());
    std::copy (x.begin(), x.end(), fx.begin());
    std::copy (y.begin(), y.end(), fy.begin());
    cerr << "Computing Ooyama FORTRAN results." << endl;
    if (vic (fx.begin(), fx.size(), wl, bc, fy.begin()))
    {
	cerr << "Done." << endl;
	ofstream vspline("ooyama.out");
	ostream_iterator<float> of(vspline, "\t ");

	float fout, foutd;
	int kdat = 1;
	float xi = fx.front();
	float xs = (fx.back() - fx.front())/2000.0;
	for (int i = 0; i < 2000; ++i, xi += xs)
	{
	    spotval_ (&xi, &kdat, &fout, &foutd);
	    *of++ = xi;
	    *of++ = fout;
	    *of++ = foutd;
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
DumpSpline (vector<datum> &x, vector<datum> &y, SplineT &spline, ostream &out)
{
    ostream_iterator<datum> of(out, "\t ");

    for (unsigned int i = 0; i < x.size(); ++i)
    {
	*of++ = x[i];
	*of++ = y[i];
	*of++ = spline.evaluate (x[i]);
	out << endl;
    }
}



void
EvalSpline (SplineT &spline, ostream &out)
{
    ostream_iterator<datum> of(out, "\t ");

    datum x = spline.Xmin();
    double xs = (spline.Xmax() - x) / 2000.0;
    for (; x <= spline.Xmax(); x += xs)
    {
	*of++ = x;
	*of++ = spline.evaluate (x);
	*of++ = spline.slope (x);
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
vic (float *xt, int nxp, double wl, int bc, float *y)
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
