//
// BSpline.cxx: implementation of the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////

//#include <math.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <ios>
#include <iomanip>

#include <assert.h>

//using namespace std;

//#include <tnt.h>
//#include <cmat.h>
//#include <lu.h>

#include "BandedMatrix.h"

namespace my {
template <class T> 
inline T abs(const T t) { return (t < 0) ? -t : t; }

template <class T>
inline const T& min (const T& a, const T& b) { return (a < b) ? a : b; }

template <class T>
inline const T& max (const T& a, const T& b) { return (a > b) ? a : b; }
}


/*
 * This is a modified version of R. Pozo's LU_factor template procedure from
 * the Template Numerical Toolkit.  It was modified to limit pivot searching
 * in the case of banded diagonal matrices.  The extra parameter BANDS
 * is the number of bands below the diagonal.
 */
template <class Matrix, class VectorSubscript>
int LU_factor_banded ( Matrix &A, VectorSubscript &indx, int bands)
{
    // assert(A.lbound() == 1);                // currently for 1-offset
    // assert(indx.lbound() == 1);             // vectors and matrices

	Matrix::size_type M = A.num_rows();
	Matrix::size_type N = A.num_cols();

    if (M == 0 || N==0) return 0;
	indx.assign (M);

	Matrix::size_type i=0,j=0,k=0,jp=0;
    Matrix::element_type t;
	Matrix::size_type minMN =  (M < N ? M : N) ;        // min(M,N);

    for (j=1; j<= minMN; j++)
    {

        // find pivot in column j and  test for singularity.

        jp = j;
        t = abs(A(j,j));
        for (i=j+1; i<=my::min(j+bands,M); i++)
            if ( abs(A(i,j)) > t)
            {
                jp = i;
                t = abs(A(i,j));
            }

        indx[j-1] = jp;

        // jp now has the index of maximum element 
        // of column j, below the diagonal

        if ( A(jp,j) == 0 )                 
            return 1;       // factorization failed because of zero pivot


        if (jp != j)            // swap rows j and jp
            for (k=1; k<=N; k++)
            {
                t = A(j,k);
                A(j,k) = A(jp,k);
                A(jp,k) =t;
            }

        if (j<M)                // compute elements j+1:M of jth column
        {
            // note A(j,j), was A(jp,p) previously which was
            // guarranteed not to be zero (Label #1)
            //
            Matrix::element_type recp =  1.0 / A(j,j);

            for (k=j+1; k<=M; k++)
                A(k,j) *= recp;
        }


        if (j < minMN)
        {
            // rank-1 update to trailing submatrix:   E = E - x*y;
            //
            // E is the region A(j+1:M, j+1:N)
            // x is the column vector A(j+1:M,j)
            // y is row vector A(j,j+1:N)

			Matrix::size_type ii,jj;

            for (ii=j+1; ii<=M; ii++)
                for (jj=j+1; jj<=N; jj++)
                    A(ii,jj) -= A(ii,j)*A(j,jj);
        }
    }

    return 0;
}   



template <class Matrix, class Vector, class VectorSubscripts>
int LU_solve_banded(const Matrix &A, const VectorSubscripts &indx, Vector &b)
{
    //assert(A.lbound() == 1);                // currently for 1-offset
    //assert(indx.lbound() == 1);             // vectors and matrices
    //assert(b.lbound() == 1);

	Matrix::size_type i,ii=0,ip,j;
	Matrix::size_type n = A.num_rows();
    Matrix::element_type sum = 0.0;

    for (i=1;i<=n;i++) 
    {
        ip=indx[i-1];
        sum=b[ip-1];
        b[i-1]=b[i-1];
        if (ii)
            for (j=ii;j<=i-1;j++) 
                sum -= A(i,j)*b[j-1];
        else if (sum) ii=i;
            b[i-1]=sum;
    }
    for (i=n;i>=1;i--) 
    {
        sum=b[i-1];
        for (j=i+1;j<=n;j++) 
            sum -= A(i,j)*b[j-1];
        b[i-1]=sum/A(i,i);
    }

    return 0;
}


#include "BSpline.h"
#include "BSplineSolver.h"

// Our private state structure, which also hides our use
// of TNT for matrices.

//typedef C_matrix<float> MatrixT;
typedef BandedMatrix<float> MatrixT;

struct BSplineBaseP 
{
	MatrixT Q;					// Holds P+Q
	BSplineSolver<MatrixT> solver;
	MatrixT LU;					// LU factorization of PQ
	std::vector<MatrixT::element_type> index;
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
	//Copy (bb.base->X.begin(), bb.NX, bb.waveLength);
	xmin = bb.xmin;
	xmax = bb.xmax;
	alpha = bb.alpha;
	waveLength = bb.waveLength;
	DX = bb.DX;
	M = bb.M;
	NX = base->X.size();
}


BSplineBase::BSplineBase (const float *x, int nx, float wl, int bc) : 
	K(1), base(new BSplineBaseP)
{
	setDomain (x, nx, wl, bc);
}


// Methods


bool
BSplineBase::setDomain (const float *x, int nx, float wl, int bc)
{
	// If called while we have an existing array, release it.
	//Copy (x, nx, wl);
	OK = false;
	BC = bc;
	waveLength = wl;
	if (nx <= 0 || x == 0)
	{
		return false;
	}
		
	// Copy the x array into our storage.
	base->X.resize (nx);
	std::copy (x, x+nx, base->X.begin());
	NX = base->X.size();

	// The Setup() method determines the number and size of node intervals.
	if (Setup())
	{
		if (Debug) 
			cerr << "Using M node intervals: " << M << " of length DX: "
				 << DX << endl;

		// Now we can calculate alpha and our Q matrix
		alpha = Alpha (waveLength);
		if (Debug)
		{
			cerr << "Alpha: " << alpha << endl;
			cerr << "Calculating Q..." << endl;
		}
		calculateQ ();
		if (Debug && M < 30)
		{
			cerr.fill(' ');
			cerr.precision(2);
			cerr.setf(std::ios_base::right, std::ios_base::adjustfield);
			cerr.width(5);
			cerr << base->Q << endl;
		}
	
		if (Debug) cerr << "Calculating P..." << endl;
		addP ();
		if (Debug && M < 30)
		{
			cerr << "Array Q after addition of P." << endl;
			cerr << base->Q;
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
	float a = (float) (wl / (2 * PI));
	a *= a;				// a^2
	if (K == 2)
		a *= a;			// a^4
	else if (K == 3)
		a *= a * a;		// a^6
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
	float z = abs((x - xm) / DX);
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



#if 0
template <class T>
Vector<T>& operator *= (Vector<T> &v, const T mult)
{
	Subscript N = v.dim();
	Vector<T>::iterator vi;
	
	for (vi = v.begin(); vi != v.end(); ++vi)
	{
		*vi *= mult;
	}
	return v;
}
#endif


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
 * Return the integral of the product of the basis function derivative restricted
 * to the node domain, 0 to M.
 */
{
	// At present Q is hardcoded for the first derivative
	// filter constraint and the type 1 boundary constraint.

	// These are the products of the first derivative of the
	// normalized basis functions
	// given a distance m nodes apart, qparts[m], 0 <= m <= 3
	// Each column is the integral over each unit domain, -2 to 2
	static const float qparts[4][4] = 
	{
		0.11250f,   0.63750f,   0.63750f,   0.11250f,
		0.00000f,   0.13125f,  -0.54375f,   0.13125f,
		0.00000f,   0.00000f,  -0.22500f,  -0.22500f,
		0.00000f,   0.00000f,   0.00000f,  -0.01875f
	};

	// Simply the sums of each row above, for interior nodes
	// where the integral is not restricted by the domain.
	static const float qinterior[4] =
	{
		1.5f,	-0.28125f,	-0.450f,	-0.01875f
	};

	if (m1 > m2)
		std::swap (m1, m2);

	if (m2 - m1 > 3)
		return 0.0;

	float q = 0.0;
	for (int m = my::max (m1-2,0); m < my::min (m1+2, M); ++m)
		q += qparts[m2-m1][m-m1+2];
	return q * DX * alpha;
}


template <class T> 
void setup (T &matrix, int n);

//void setup<> (C_matrix<float> &matrix, int n) { matrix.newsize (n, n); }

void setup<> (BandedMatrix<float> &matrix, int n) { matrix.setup (n, 3); }

void
BSplineBase::calculateQ ()
{
	MatrixT &Q = base->Q;
	setup (Q, M+1);
	// Q.newsize (M+1,M+1);
	Q = 0.0;

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
	MatrixT P;
	setup (P, M+1);
	P = 0.0;
	std::vector<float> &X = base->X;

	// For each data point, sum the product of the nearest, non-zero Basis nodes
	int m, n, i;
	for (i = 0; i < NX; ++i)
	{
		// Which node does this put us in?
		float x = X[i];
		m = (x - xmin) / DX;

		// Loop over the upper triangle of the 5x5 array of basis functions,
		// and add in the products on each side of diagonal.
		for (m = my::max(0, m-2); m <= my::min(M, m+2); ++m)
		{
			float pn;
			float pm = Basis (m, x);
			float sum = pm * pm * DX;
			P[m][m] += sum;
			for (n = m+1; n <= my::min(M, m+3); ++n)
			{
				pm = Basis (m, x);
				pn = Basis (n, x);
				sum = pm * pn * DX;
				P[m][n] += sum;
				P[n][m] += sum;
			}
		}
	}

	base->Q += P;
}



bool
BSplineBase::factor ()
{	
//#if 0
	base->index.assign (M+1);
	base->LU = base->Q;

    if (LU_factor_banded (base->LU, base->index, 3) != 0)
    {
        if (Debug) cerr << "LU_factor() failed." << endl;
		return false;
    }
//#endif

	if (! base->solver.upper (base->Q))
	{
		if (Debug) cerr << "BSplineSolver::upper failed." << endl;
		return false;
	}
	if (Debug && M < 30)
		cerr << *base->solver.matrix();
	return true;
}

	

inline int 
BSplineBase::Ratio (int &ni, float &deltax, float &ratiof,
					float *ratiod)
{
		deltax = (xmax - xmin) / ni;
		ratiof = deltax / waveLength;
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

	// Minimum acceptable number of nodes per cutoff wavelength
	static const float fmin = 2.0;

	int ni = 9;		// Number of node intervals
	float deltax;
	float ratiof;	// Nodes per wavelength for current deltax
	float ratiod;	// Points per node interval

	do {
		if (! Ratio (++ni, deltax, ratiof))
			return false;
	}
	while (ratiof > fmin);

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

	assert (base->Nodes.size() == M+1);
	return base->Nodes.begin();
}


#include <iterator>
#include <iostream>

template <class C>
ostream &operator<< (ostream &out, C &c)
{
	for (C::iterator it = c.begin(); it < c.end(); ++it)
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
	std::vector<float> A2;
};


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BSpline::BSpline (BSplineBase &bb, const float *y) :
	BSplineBase(bb), s(new BSplineP)
{
	if (! OK)
		return;

	// Given an array of data points over x and its precalculated
	// P+Q matrix, calculate the b vector and solve for the coefficients.

	std::vector<float> B(M+1);

	// Find the mean of these data
	mean = 0.0;
	int i;
	for (i = 0; i < NX; ++i)
	{
		mean += y[i];
	}
	mean = mean / (float)NX;

	int m, j;
	for (m = 0; m < M+1; ++m)
	{
		float sum = 0.0;
		for (j = 0; j < NX; ++j)
		{
			sum += (y[j] - mean) * Basis (m, base->X[j]);
		}
		B[m] = sum * DX;
	}

	std::vector<float> &luA = s->A;
	std::vector<float> &sA = s->A2;

	// Now solve for the A vector.
	if (! base->solver.solve (B, sA))
	{
		cerr << "Solver failed." << endl;
		exit(1);
	}
//#if 0
	luA = B;
    if (LU_solve_banded (base->LU, base->index, luA) != 0)
    {
        cerr << "LU_Solve() failed." << endl;
        exit(1);
    }
//#endif
	if (Debug && M < 30)
	{
	    cerr << "Solution a for (P+Q)a = b" << endl;
		//cerr << " b: " << B << endl;
		cerr << "solver a: " << sA << endl;
		cerr << "    lu a: " << luA << endl;

	    //cerr << "A*x should be the vector b, ";
		//cerr << "residual [s->A*x - b]: " << endl;
		//cerr << matmult(base->Q, s->A) - B << endl;
	}
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
		for (int i = 0; i <= M; ++i)
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
