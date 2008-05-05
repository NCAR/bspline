// -*- mode: c++; c-basic-offset: 4; -*-
//
// $Id: driver.cpp 4168 2006-01-08 18:39:01Z martinc $
//
// Simple test driver for the bspline interface.
//
/*
 * Copyright (c) 1998,1999,2008
 * University Corporation for Atmospheric Research, UCAR
 *
 * Permission to use, copy, modify, distribute and sell this software and
 * its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and that
 * both that copyright notice and this permission notice appear in
 * supporting documentation.  UCAR makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 * 
 * Note from the author:
 *
 * Where possible, you are encouraged to follow the GNU General Public
 * License, or at least the spirit of the license, for the distribution and
 * licensing of this software and any derived works.  See
 * http://www.gnu.org/copyleft/gpl.html.
 */

#include "BSpline/BSpline.h"
#include "options.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

typedef double datum;
typedef BSplinePlus<datum> SplineT;
typedef BSplineBase<datum> SplineBase;
typedef BSpline<double> SplineD;

void DumpSpline(vector<datum> &x,
                vector<datum> &y,
                SplineT &spline,
                ostream &out);
static void EvalSpline(SplineT &spline,
                       ostream &out);
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

///////////////////////////////////////////////////////////////////////////////
static const char
        * optv[] =
            { "s:step        <step interval>",
                    "w:wavelength  <spline wavelength (required)>",
                    "b:bcdegree    <bc derivative degree (0,1,2)>",
                    "n:nodes       <specify number of nodes (n)>",
                    "m:blendmode   <blending mode (none, start, end, both)> (requires blendlength)",
                    "l:blendlength <length> (requires blendmode)",
                    "d|debug       <enable diagnostic output>"
                    "h|help        <print this help>",
                    NULL };

///////////////////////////////////////////////////////////////////////////////
void parseCommandLine(int argc,
                      char * const argv[],
                      int* step,
                      double* wavelength,
                      int* bc,
                      int* num_nodes,
                      BSplinePlus<double>::BLENDMODE* blendmode,
                      double* blendlength,
                      bool* debug) {

    // initialize the optional parameters
    *step = 0;
    *bc = SplineBase::BC_ZERO_SECOND;
    *num_nodes = 0;
    *blendmode = BSplinePlus<double>::BLENDNONE;
    *blendlength = -1;
    *debug = false;
    
    // indicate that the wavelength has not been set
    *wavelength = -1.0;

    bool err = false;
    const char *optarg;
    char optchar;
    Options opts(*argv, optv);
    OptArgvIter iter(--argc, ++argv);

    while (optchar = opts(iter, optarg) ) {
        switch (optchar) {
        case 'h': {
            opts.usage(std::cout, "");
            exit(0);
            break;
        }
        case 's': {
            if (optarg)
                *step = atoi(optarg);
            else
                err++;
            break;
        }
        case 'w': {
            if (optarg)
                *wavelength = atof(optarg);
            else
                err++;
            break;
        }
        case 'b': {
            if (optarg) {
                int degree = atoi(optarg);
                switch (degree) {
                case 0:
                    *bc = SplineBase::BC_ZERO_ENDPOINTS;
                    break;
                case 1:
                    *bc = SplineBase::BC_ZERO_FIRST;
                    break;
                case 2:
                default:
                    *bc = SplineBase::BC_ZERO_SECOND;
                    break;
                }
            } else
                err++;
            break;
        }
        case 'n': {
            if (optarg)
                *num_nodes = atoi(optarg);
            else
                err++;
            break;
        }
        case 'm': {
            if (optarg)
                if (!strcmp(optarg, "none"))
                    *blendmode = BSplinePlus<double>::BLENDNONE;
                else if (!strcmp(optarg, "start"))
                    *blendmode = BSplinePlus<double>::BLENDSTART;
                else if (!strcmp(optarg, "finish"))
                    *blendmode = BSplinePlus<double>::BLENDFINISH;
                else if (!strcmp(optarg, "both"))
                    *blendmode = BSplinePlus<double>::BLENDBOTH;
                else
                    err++;
            else
                err++;
            break;
        }
        case 'l': {
            if (optarg)
                *blendlength = atof(optarg);
            else
                err++;
            break;
        }
        case 'd': {
            *debug = true;
            break;
        }
        default: {
            ++err;
        }
            break;
        }
    }

    // wavelength must be supplied
    if (*wavelength < 0)
        err++;
    
    // if blending is requested, there must be a blendlength
    if (*blendmode != BSplinePlus<double>::BLENDNONE && *blendlength == -1)
        err++;

    // if blendlength is ot specified, blending must not be requested
    if (*blendmode == BSplinePlus<double>::BLENDNONE && *blendlength != -1)
        err++;

    if (err) {
        opts.usage(std::cerr, "");
        exit(1);
    }

}

///////////////////////////////////////////////////////////////////////////////
int main(int argc,
         char *argv[]) {

    // Some basics first
    cout << "BSpline interface version: " << SplineBase::IfaceVersion() << endl;
    cout << "BSpline implementation version: " << SplineBase::ImplVersion()
            << endl;

    int step;
    double wavelength;
    int bc;
    int num_nodes;
    BSplinePlus<double>::BLENDMODE blendmode;
    double blendlength;
    bool debug;

    parseCommandLine(argc,
                     argv,
                      &step,
                      &wavelength,
                      &bc,
                      &num_nodes,
                      &blendmode,
                      &blendlength,
                      &debug);

    if (debug) {
        cout << "Using step interval " << step << ", cutoff frequency "
            << wavelength << ", number of nodes " << num_nodes
            << ", and boundary condition type " << bc << "\n";
    }

    // Read the x and y pairs from stdin
    vector<datum> x;
    vector<datum> y;
    datum f, g;
    datum base = 0;
    int i = 0;

    while (cin >> f) {
        if (++i == 1) {
            base = f;
        }
        f -= base;

        if (cin >> g) {
            x.push_back(f);
            y.push_back(g);
            // cout << x.back() << " " << y.back() << endl;
        } else
            break;
    }

    // Subsample the arrays
    vector<datum>::iterator xi = x.begin(), yi = y.begin();
    vector<datum>::iterator xo = x.begin(), yo = y.begin();
    if (step > 1) {
        while (xo < x.end() && yo < y.end()) {
            //cout << *xo;
            //cout << " ";
            //cout << *yo;
            //cout << endl;
            *xi++ = *xo;
            *yi++ = *yo;
            xo += step;
            yo += step;
        }
        x.resize(xi - x.begin());
        y.resize(yi - y.begin());
    }

    // Create our bspline base on the X vector with a simple 
    // wavelength.
    if (debug)
        SplineT::Debug(1);
    SplineT spline(&x[0],
                   x.size(),
                   &y[0],
                   wavelength,
                   bc,
                   num_nodes,
                   blendmode,
                   blendlength);
    if (spline.ok()) {
        // And finally write the curve to a file
        ofstream fspline("input.out");
        DumpSpline(x, y, spline, fspline);
        ofstream fcurve("spline.out");
        EvalSpline(spline, fcurve);
    } else
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

void DumpSpline(vector<datum> &x,
                vector<datum> &y,
                SplineT &spline,
                ostream &out) {
    ostream_iterator<datum> of(out, "\t ");
    datum variance = 0;

    for (unsigned int i = 0; i < x.size(); ++i) {
        *of++ = x[i];
        *of++ = y[i];
        datum ys = spline.evaluate(x[i]);
        *of++ = ys;
        datum slope = spline.slope(x[i]);
        *of++ = slope;
        variance += (ys - y[i])*(ys - y[i]);
        out << endl;
    }
    variance /= (datum)x.size();
    cerr << "Variance: " << variance << endl;
}

void EvalSpline(SplineT &spline,
                ostream &out) {
    ostream_iterator<datum> of(out, "\t ");

    datum x = spline.Xmin();
    double xs = (spline.Xmax() - x) / 2000.0;
    for (; x <= spline.Xmax(); x += xs) {
        *of++ = x;
        *of++ = spline.evaluate(x);
        *of++ = spline.slope(x);
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
        * KYBC,KYBC,YBCWL,YBCWL,IERR)
IF (IERR.NE.0) GOTO 900
DO 500 I = 3,NDIMX-2
IF (X(I).LE.-999.) GOTO 500
CALL SPOTVAL(XT(I),KDAT,FOUT,FOUTD)
IF (ABS(X(I)-FOUT).GT.DEVX) THEN
NBAD = NBAD+1
IBAD(NBAD) = I
NBADX = NBADX+1
ENDIF
500 CONTINUE
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
    int kybc = bc; // endpoint boundary conditions
    int kdat = 1;
    vector<float> xw(nxp, 1.0); // relative weights all set to 1

    ierr = 1;
    vicsetup_ (xt, &nxp, &ydcwl, &nx, &ynb, &ynt, &fmin, &ierr, &echo);
    if (ierr)
    return (false);

    ierr = 1;
    vicspl_ (xt, xw.begin(), /*xdat*/y, &nxp, &nxp, &kdat,
            &ynb, &ynt, &nx,
            &ydcwl, &kybc, &kybc, &ydcwl, &ydcwl, &ierr);
    if (ierr)
    return (false);

    return true;
}

#endif /* OOYAMA */
