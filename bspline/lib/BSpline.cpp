// -*- mode: c++; c-basic-offset: 4; -*-
//
// BSpline.cxx: implementation of the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////

#include <vector>
#include <algorithm>
#include <iterator>

#if WIN32
#include <iostream>
#endif

#include <assert.h>

using namespace std;

/*
 * These conflict with Windows if not in my namespace, but egcs does not yet
 * support namespaces.  And I can't use Windows' abs() because
 * it truncs the arg to an int.
 */
#if WIN32
namespace my {
#endif
template <class T> 
inline T abs(const T t) { return (t < 0) ? -t : t; }

#if WIN32
template <class T>
inline const T& min (const T& a, const T& b) { return (a < b) ? a : b; }

template <class T>
inline const T& max (const T& a, const T& b) { return (a > b) ? a : b; }
#endif

#if WIN32
}

using my::min;
using my::max;
using my::abs;
#endif

#include "BandedMatrix.h"
#include "BSplineLU.h"
#include "BSpline.h"

template <class T> 
void setup (T &matrix, int n);

void setup<> (BandedMatrix<float> &matrix, int n) { matrix.setup (n, 3); }

// Our private state structure, which hides our use of some matrix
// template classes.

typedef BandedMatrix<float> MatrixT;

struct BSplineBaseP 
{
    MatrixT Q;				// Holds P+Q and its factorization
    std::vector<float> X;
    std::vector<float> Nodes;
};


// For now, hardcoding type 1 boundary conditions, 
// which constrains the derivative to zero at the endpoints.
const float BSplineBase::BoundaryConditions[3][4] = 
{ 
    //	0		1		M-1		M
    {	-4,		-1,		-1,		-4 },
    {	0,		1,		1,		0 },
    {	2,		-1,		-1,		2 }
};

const double BSplineBase::PI = 3.1415927;

bool BSplineBase::Debug = false;

const char *
BSplineBase::ImplVersion()
{
    return ("$Id$");
}

const char *
BSplineBase::IfaceVersion()
{
    return (_BSPLINEBASE_IFACE_ID);
}

	
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


BSplineBase::~BSplineBase()
{
    delete base;
}


// This is a member-wise copy except for replacing our
// private base structure with the source's, rather than just copying
// the pointer.  But we use the compiler's default copy constructor for
// constructing our BSplineBaseP.
BSplineBase::BSplineBase (const BSplineBase &bb) : 
    K(bb.K), BC(bb.BC), OK(bb.OK), base(new BSplineBaseP(*bb.base))
{
    xmin = bb.xmin;
    xmax = bb.xmax;
    alpha = bb.alpha;
    waveLength = bb.waveLength;
    DX = bb.DX;
    M = bb.M;
    NX = base->X.size();
}


BSplineBase::BSplineBase (const float *x, int nx, float wl, int bc) : 
    K(2), OK(false), base(new BSplineBaseP)
{
    setDomain (x, nx, wl, bc);
}


// Methods


bool
BSplineBase::setDomain (const float *x, int nx, float wl, int bc)
{
    if (nx <= 0 || x == 0 || wl < 0 || bc < 0 || bc > 2)
    {
	return false;
    }
    OK = false;
    waveLength = wl;
    BC = bc;
		
    // Copy the x array into our storage.
    base->X.resize (nx);
    std::copy (x, x+nx, base->X.begin());
    NX = base->X.size();

    // The Setup() method determines the number and size of node intervals.
    if (Setup())
    {
	if (Debug) 
	{
	    cerr << "Using M node intervals: " << M << " of length DX: "
		 << DX << endl;
	    cerr << "X min: " << xmin << " ; X max: " << xmax << endl;
	    cerr << "Data points per interval: " << (float)NX/(float)M << endl;
	    cerr << "Derivative constraint degree: " << K << endl;
	}

	// Now we can calculate alpha and our Q matrix
	alpha = Alpha (waveLength);
	if (Debug)
	{
	    cerr << "Cutoff wavelength: " << waveLength << " ; "
		 << "Alpha: " << alpha << endl;
	    cerr << "Calculating Q..." << endl;
	}
	calculateQ ();
	if (Debug && M < 30)
	{
	    cerr.fill(' ');
	    cerr.precision(2);
	    cerr.width(5);
	    cerr << base->Q << endl;
	}
	
	if (Debug) cerr << "Calculating P..." << endl;
	addP ();
	if (Debug)
	{
	    cerr << "Done." << endl;
	    if (M < 30)
	    {
		cerr << "Array Q after addition of P." << endl;
		cerr << base->Q;
	    }
	}

	// Now perform the LU factorization on Q
	if (Debug) cerr << "Beginning LU factoring of P+Q..." << endl;
	if (! factor ())
	{
	    if (Debug) cerr << "Factoring failed." << endl;
	}
	else
	{
	    if (Debug) cerr << "Done." << endl;
	    OK = true;
	}
    }
    return OK;
}



/*
 * Calculate the alpha parameter given a wavelength.
 */
float
BSplineBase::Alpha (float wl)
{
    // K is the degree of the derivative constraint: 1, 2, or 3
    float a = (float) (wl / (2 * PI * DX));
    a *= a;			// a^2
    if (K == 2)
	a = a * a;		// a^4
    else if (K == 3)
	a = a * a * a;		// a^6
    return a;
}


/*
 * Return the correct beta value given the node index.  The value depends
 * on the node index and the current boundary condition type.
 */
inline float
BSplineBase::Beta (int m)
{
    if (m > 1 && m < M-1)
	return 0.0;
    if (m >= M-1)
	m -= M-3;
    assert (0 <= BC && BC <= 2);
    assert (0 <= m && m <= 3);
    return BoundaryConditions[BC][m];
}



/*
 * Given an array of y data points defined over the domain
 * of x data points in this BSplineBase, create a BSpline
 * object which contains the smoothed curve for the y array.
 */
BSpline *
BSplineBase::apply (const float *y)
{
    BSpline *spline = new BSpline (*this, y);

    return (spline);
}


/*
 * Evaluate the closed basis function at node m for value x,
 * using the parameters for the current boundary conditions.
 */
float
BSplineBase::Basis (int m, float x)
{
    float y = 0;
    float xm = xmin + (m * DX);
    float z = /*my::*/abs((float)(x - xm) / (float)DX);
    if (z < 2.0)
    {
	z = 2 - z;
	y = 0.25 * (z*z*z);
	z -= 1.0;
	if (z > 0)
	    y -= (z*z*z);
    }

    // Boundary conditions, if any, are an additional addend.
    if (m == 0 || m == 1)
	y += Beta(m) * Basis (-1, x);
    else if (m == M-1 || m == M)
	y += Beta(m) * Basis (M+1, x);

    return y;
}




MatrixT &operator += (MatrixT &A, const MatrixT &B)
{
    MatrixT::size_type M = A.num_rows();
    MatrixT::size_type N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    MatrixT::size_type i,j;
    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            A[i][j] += B[i][j];
    return A;
}



float
BSplineBase::qDelta (int m1, int m2)
/*
 * Return the integral of the product of the basis function derivative
 * restricted to the node domain, 0 to M.
 */
{
    // At present Q is hardcoded for the first derivative
    // filter constraint and the type 1 boundary constraint.

    // These are the products of the first derivative of the
    // normalized basis functions
    // given a distance m nodes apart, qparts[m], 0 <= m <= 3
    // Each column is the integral over each unit domain, -2 to 2
    static const float qparts[3][4][4] = 
    {
	{
	    { 0.11250f,   0.63750f,   0.63750f,   0.11250f },
	    { 0.00000f,   0.13125f,  -0.54375f,   0.13125f },
	    { 0.00000f,   0.00000f,  -0.22500f,  -0.22500f },
	    { 0.00000f,   0.00000f,   0.00000f,  -0.01875f }
	},
	{
	    { 0.75000f,   2.25000f,   2.25000f,   0.75000f },
	    { 0.00000f,  -1.12500f,  -1.12500f,  -1.12500f },
	    { 0.00000f,   0.00000f,   0.00000f,   0.00000f },
	    { 0.00000f,   0.00000f,   0.00000f,   0.37500f }
	},
	{
	    { 2.25000f,  20.25000f,  20.25000f,   2.25000f },
	    { 0.00000f,  -6.75000f, -20.25000f,  -6.75000f },
	    { 0.00000f,   0.00000f,   6.75000f,   6.75000f },
	    { 0.00000f,   0.00000f,   0.00000f,  -2.25000f }
	}
    };

    if (m1 > m2)
	std::swap (m1, m2);

    if (m2 - m1 > 3)
	return 0.0;

    float q = 0.0;
    for (int m = /*my::*/max (m1-2,0); m < /*my::*/min (m1+2, M); ++m)
	q += qparts[K-1][m2-m1][m-m1+2];
    return q * alpha;
}



void
BSplineBase::calculateQ ()
{
    MatrixT &Q = base->Q;
    setup (Q, M+1);
    Q = 0.0;
    if (alpha == 0)
	return;

    // First fill in the q values without the boundary constraints.
    int i;
    for (i = 0; i <= M; ++i)
    {
	Q[i][i] = qDelta(i,i);
	for (int j = 1; j < 4 && i+j <= M; ++j)
	{
	    Q[i][i+j] = Q[i+j][i] = qDelta (i, i+j);
	}
    }

    // Now add the boundary constraints:
    // First the upper left corner.
    float b1, b2, q;
    for (i = 0; i <= 1; ++i)
    {
	b1 = Beta(i);
	for (int j = i; j < i+4; ++j)
	{
	    b2 = Beta(j);
	    assert (j-i >= 0 && j - i < 4);
	    q = 0.0;
	    if (i+1 < 4)
		q += b2*qDelta(-1,i);
	    if (j+1 < 4)
		q += b1*qDelta(-1,j);
	    q += b1*b2*qDelta(-1,-1);
	    Q[j][i] = (Q[i][j] += q);
	}
    }

    // Then the lower right
    for (i = M-1; i <= M; ++i)
    {
	b1 = Beta(i);
	for (int j = i - 3; j <= i; ++j)
	{
	    b2 = Beta(j);
	    q = 0.0;
	    if (M+1-i < 4)
		q += b2*qDelta(i,M+1);
	    if (M+1-j < 4)
		q += b1*qDelta(j,M+1);
	    q += b1*b2*qDelta(M+1,M+1);
	    Q[j][i] = (Q[i][j] += q);
	}
    }
}




void
BSplineBase::addP ()
{
    // Add directly to Q's elements
    MatrixT &P = base->Q;
    std::vector<float> &X = base->X;

    // For each data point, sum the product of the nearest, non-zero Basis
    // nodes
    int mx, m, n, i;
    for (i = 0; i < NX; ++i)
    {
	// Which node does this put us in?
	float &x = X[i];
	mx = (int)((x - xmin) / DX);

	// Loop over the upper triangle of nonzero basis functions,
	// and add in the products on each side of the diagonal.
	for (m = /*my::*/max(0, mx-1); m <= /*my::*/min(M, mx+2); ++m)
	{
	    float pn;
	    float pm = Basis (m, x);
	    float sum = pm * pm;
	    P[m][m] += sum;
	    for (n = m+1; n <= /*my::*/min(M, mx+2); ++n)
	    {
		pn = Basis (n, x);
		sum = pm * pn;
		P[m][n] += sum;
		P[n][m] += sum;
	    }
	}
    }
}



bool
BSplineBase::factor ()
{	
    MatrixT &LU = base->Q;

    if (LU_factor_banded (LU, 3) != 0)
    {
        if (Debug) cerr << "LU_factor() failed." << endl;
	return false;
    }
    if (Debug && M < 30)
	cerr << "LU decomposition: " << endl << LU << endl;
    return true;
}

	

inline int 
BSplineBase::Ratio (int &ni, float &deltax, float &ratiof,
		    float *ratiod)
{
    deltax = (xmax - xmin) / ni;
    ratiof = waveLength / deltax;
    float rd = (float) NX / (float) (ni + 1);
    if (ratiod)
	*ratiod = rd;
    return (rd >= 1.0);
}


/*
 * Return zero if this fails, non-zero otherwise.
 */
bool BSplineBase::Setup()
{
    std::vector<float> &X = base->X;
	
    // Find the min and max of the x domain
    xmin = X[0];
    xmax = X[0];

    int i;
    for (i = 1; i < NX; ++i)
    {
	if (X[i] < xmin)
	    xmin = X[i];
	else if (X[i] > xmax)
	    xmax = X[i];
    }

    if (waveLength > xmax - xmin)
    {
	return (false);
    }

    int ni = 9;		// Number of node intervals (NX - 1)
    float deltax;

    if (waveLength == 0)	// Allows turning off frequency constraint
    {
	ni = NX;
	deltax = (xmax - xmin) / (float)NX;
    }
    else
    {
	// Minimum acceptable number of node intervals per cutoff wavelength.
	static const float fmin = 2.0;

	float ratiof;	// Nodes per wavelength for current deltax
	float ratiod;	// Points per node interval

	do {
	    if (! Ratio (++ni, deltax, ratiof))
		return false;
	}
	while (ratiof < fmin);

	// Tweak the estimates obtained above
	do {
	    if (! Ratio (++ni, deltax, ratiof, &ratiod) || 
		ratiof > 15.0)
	    {
		Ratio (--ni, deltax, ratiof);
		break;
	    }
	}
	while (ratiof < 4 || ratiod > 2.0);
    }

    // Store the calculations in our state
    M = ni;
    DX = deltax;

    return (true);
}


const float *
BSplineBase::nodes (int *nn)
{
    if (base->Nodes.size() == 0)
    {
	base->Nodes.reserve (M+1);
	for (int i = 0; i <= M; ++i)
	{
	    base->Nodes.push_back ( xmin + (i * DX) );
	}
    }

    if (nn)
	*nn = base->Nodes.size();

    assert (base->Nodes.size() == (unsigned)(M+1));
    return base->Nodes.begin();
}



ostream &operator<< (ostream &out, const vector<float> &c)
{
    for (vector<float>::const_iterator it = c.begin(); it < c.end(); ++it)
	out << *it << ", ";
    out << endl;
    return out;
}



//////////////////////////////////////////////////////////////////////
// BSpline Class
//////////////////////////////////////////////////////////////////////

struct BSplineP
{
    std::vector<float> spline;
    std::vector<float> A;
};


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


/*
 * This BSpline constructor constructs and sets up a new base and 
 * solves for the spline curve coeffiecients all at once.
 */
BSpline::BSpline (const float *x, int nx, const float *y,
		  float wl, int bc_type = BC_ZERO_SECOND) :
    BSplineBase(x, nx, wl, bc_type), s(new BSplineP)
{
    solve (y);
}



/*
 * Create a new spline given a BSplineBase.
 */
BSpline::BSpline (BSplineBase &bb, const float *y) :
    BSplineBase(bb), s(new BSplineP)
{
    solve (y);
}



/*
 * (Re)calculate the spline for the given set of y values.
 */
bool
BSpline::solve (const float *y)
{
    if (! OK)
	return false;

    // Any previously calculated curve is now invalid.
    s->spline.clear ();
    OK = false;

    // Given an array of data points over x and its precalculated
    // P+Q matrix, calculate the b vector and solve for the coefficients.
    std::vector<float> &B = s->A;
    std::vector<float> &A = s->A;
    A.clear ();
    A.resize (M+1);

    if (Debug) cerr << "Solving for B..." << endl;

    // Find the mean of these data
    mean = 0.0;
    int i;
    for (i = 0; i < NX; ++i)
    {
	mean += y[i];
    }
    mean = mean / (float)NX;
    if (Debug)
	cerr << "Mean for y: " << mean << endl;

    int mx, m, j;
    for (j = 0; j < NX; ++j)
    {
	// Which node does this put us in?
	float &xj = base->X[j];
	float yj = y[j] - mean;
	mx = (int)((xj - xmin) / DX);

	for (m = max(0,mx-1); m <= min(mx+2,M); ++m)
	{
	    B[m] += yj * Basis (m, xj);
	}
    }

    if (Debug && M < 30)
    {
	cerr << "Solution a for (P+Q)a = b" << endl;
	cerr << " b: " << B << endl;
    }

    // Now solve for the A vector in place.
    if (LU_solve_banded (base->Q, A) != 0)
    {
	if (Debug)
	    cerr << "LU_Solve() failed." << endl;
    }
    else
    {
	OK = true;
	if (Debug) cerr << "Done." << endl;
	if (Debug && M < 30)
	{
	    cerr << " a: " << A << endl;
	    cerr << "LU factor of (P+Q) = " << endl << base->Q << endl;
	}
    }
    return (OK);
}



BSpline::~BSpline()
{
    delete s;
}



float BSpline::coefficient (int n)
{
    if (OK)
	if (0 <= n && n <= M)
	    return s->A[n];
    return 0;
}



float BSpline::evaluate (float x)
{
    float y = 0;
    if (OK)
    {
	int n = (int)((x - xmin)/DX);
	for (int i = max(0,n-1); i <= min(M,n+2); ++i)
	{
	    y += s->A[i] * Basis (i, x);
	}
	y += mean;
    }
    return y;
}



const float *BSpline::curve (int *nx)
{
    if (! OK)
	return 0;

    // If we already have the curve calculated, don't do it again.
    std::vector<float> &spline = s->spline;
    if (spline.size() == 0)
    {
	spline.reserve (M+1);
	for (int n = 0; n <= M; ++n)
	{
	    float x = xmin + (n * DX);
	    spline.push_back (evaluate (x));
	}
    }

    if (nx)
	*nx = spline.size();
    return spline.begin();
}
