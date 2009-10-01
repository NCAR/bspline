/************************************************************************
 * Copyright 2009 University Corporation for Atmospheric Research.
 * All rights reserved.
 *
 * Use of this code is subject to the standard BSD license:
 *
 *  http://www.opensource.org/licenses/bsd-license.html
 *
 * See the COPYRIGHT file in the source distribution for the license text,
 * or see this web page:
 *
 *  http://www.eol.ucar.edu/homes/granger/bspline/doc/
 *
 *************************************************************************/

#include "BSpline.h"

#include "options.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

typedef double datum;
typedef BSpline<datum> SplineT;
typedef BSplineBase<datum> SplineBase;

void DumpSpline(vector<datum> &x,
                vector<datum> &y,
                SplineT &spline,
                ostream* out,
                bool debug);

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
            {
                    "i:input       <input file>",
                    "o:output       <output file>",
                    "w:wavelength  <spline wavelength> (required)",
                    "s:step        <step interval>",
                    "b:bcdegree    <bc derivative degree (0,1,2)>",
                    "n:nodes       <specify number of nodes (n)>",
                    "d|debug       <enable diagnostic output>"
		    "h|help        <print this help>",
                    NULL };

///////////////////////////////////////////////////////////////////////////////
void parseCommandLine(int argc,
                      char * const argv[],
                      std::string& infile,
                      std::string& outfile,
                      int& step,
                      double& wavelength,
                      int& bc,
                      int& num_nodes,
                      bool& debug)
{

    // initialize the optional parameters
    step = 0;
    bc = SplineBase::BC_ZERO_SECOND;
    num_nodes = 0;
    debug = false;

    // indicate that the wavelength has not been set
    wavelength = -1.0;
    bool err = false;
    const char *optarg;
    char optchar;
    Options opts(*argv, optv);
    OptArgvIter iter(--argc, ++argv);

    while (optchar = opts(iter, optarg) ) {
        switch (optchar)
        {
        case 'h':
            {
                // Some basics first
                cout << "BSpline interface version: "
                        << SplineBase::IfaceVersion() << endl;
                cout << "BSpline implementation version: "
                        << SplineBase::ImplVersion() << endl;
                opts.usage(std::cout, "");
                exit(0);
                break;
            }
        case 'i':
            {
                if (optarg)
                    infile = std::string(optarg);
                else
                    err++;
                break;
            }
        case 'o':
            {
                if (optarg)
                    outfile = std::string(optarg);
                else
                    err++;
                break;
            }
        case 's':
            {
                if (optarg)
                    step = atoi(optarg);
                else
                    err++;
                break;
            }
        case 'w':
            {
                if (optarg)
                    wavelength = atof(optarg);
                else
                    err++;
                break;
            }
        case 'b':
            {
                if (optarg) {
                    int degree = atoi(optarg);
                    switch (degree)
                    {
                    case 0:
                        bc = SplineBase::BC_ZERO_ENDPOINTS;
                        break;
                    case 1:
                        bc = SplineBase::BC_ZERO_FIRST;
                        break;
                    case 2:
                    default:
                        bc = SplineBase::BC_ZERO_SECOND;
                        break;
                    }
                } else
                    err++;
                break;
            }
        case 'n':
            {
                if (optarg)
                    num_nodes = atoi(optarg);
                else
                    err++;
                break;
            }
        case 'd':
            {
                debug = true;
                break;
            }
        default:
            {
                ++err;
            }
            break;
        }
    }

    // wavelength must be supplied
    if (wavelength < 0)
        err++;

    if (err) {
        opts.usage(std::cerr, "");
        exit(1);
    }
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc,
         char *argv[])
{
    std::string infile;
    std::string outfile;
    int step;
    double wavelength;
    int bc;
    int num_nodes;
    bool debug;

    parseCommandLine(argc,
                     argv,
                     infile,
                     outfile,
                     step,
                     wavelength,
                     bc,
                     num_nodes,
                     debug);

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

    // input file
    std::istream* instream;
    if (infile.size() > 0) {
        instream = new std::ifstream(infile.c_str());
        if (!*instream) {
            std::cerr << "Unable to open " << infile << "\n";
            exit(1);
        }
    }
    else
        instream = &std::cin;

    // output file
    std::ostream* outstream;
    if (outfile.size() > 0) {
        outstream = new std::ofstream(outfile.c_str());
        if (!*outstream) {
            std::cerr << "Unable to open " << outfile << "\n";
            exit(1);
        }
    }
    else
        outstream = &std::cout;

    // read data
    while (*instream >> f) {
        if (++i == 1) {
            base = f;
        }
        f -= base;

        if (*instream >> g) {
            x.push_back(f);
            y.push_back(g);
        } else
            break;
    }

    // Subsample the arrays
    vector<datum>::iterator xi = x.begin(), yi = y.begin();
    vector<datum>::iterator xo = x.begin(), yo = y.begin();
    if (step > 1) {
        while (xo < x.end() && yo < y.end()) {
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
                   num_nodes);
    if (spline.ok()) {
        // And finally write the curve to a file
        DumpSpline(x, y, spline, outstream, debug);
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
        ostream_iterator<float> of(vspline, "  ");

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

    if (infile.size())
        delete instream;
    if (outfile.size())
        delete outstream;
    return 0;
}

void DumpSpline(vector<datum> &x,
                vector<datum> &y,
                SplineT &spline,
                ostream* out,
                bool debug)
{
    int width = 15;
    
    // write column headings
    *out << "x, y, spline(x), slope(spline(x))\n";
    
    datum variance = 0;

    for (unsigned int i = 0; i < x.size(); ++i) {
        *out << x[i] << ", ";
        *out << y[i] << ", ";
        datum ys = spline.evaluate(x[i]);
        *out << ys << ", ";
        datum slope = spline.slope(x[i]);
        *out << slope;
        *out << endl;
        variance += (ys - y[i])*(ys - y[i]);
    }
    variance /= (datum)x.size();
    if (debug)
        cerr << "Variance: " << variance << endl;
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
