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
    float f, g, base;
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
    //float wl = (x.back() - x.front()) / (x.size()/3 + 1);
    //float wl = 30.0; /*secs*/
    BSplineBase::Debug = true;
    BSplineBase bb (x.begin(), x.size(), wl, 2 /*bc*/);
    if (bb.ok())
    {
	// Now apply our y values to get the smoothed curve.
	BSpline spline(bb, y.begin());

		// And finally write the curve to a file
	ofstream fspline("spline.out");
	ofstream fcurve("curve.out");
	DumpSpline (x, y, spline, fspline);
	EvalSpline (spline, fcurve);
    }
    else
	cerr << "Spline setup failed." << endl;

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
