/* -*- mode: c++; c-basic-offset: 4; -*- */
//
// BSpline.h: interface for the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _BSPLINEBASE_IFACE_ID
#define _BSPLINEBASE_IFACE_ID "$Id$"

/*
 * If including this file without the implementation (BSpline.cxx) on WIN32,
 * it assumes the implementation will come from a DLL and thus declare dllimport.
 * To explicitly instantiate the implementation in a source file, include the 
 * _IMPLEMENTATION_ file BSpline.cxx and _NOT_ this file.
 */
#if WIN32
# ifndef BSPLINE_DLL_
#  define BSPLINE_DLL_ __declspec(dllimport)
# endif
#else
# define BSPLINE_DLL_
#endif /* WIN32 */


template <class T> class BSpline;

// Opaque structure to hide our matrix implementation, ala
// Cheshire cat.
template <class T> struct BSplineBaseP;

/**
 * Base class for a spline object containing the nodes for a given domain,
 * cutoff wavelength, and boundary condition.  To smooth a single curve,
 * the BSpline interface contains a constructor which both sets up the
 * domain and solves for the spline.  Subsequent curves over the same
 * domain can be created by apply()ing them to the BSpline object, where
 * apply() is a BSplineBase method.  [See \Ref{apply}.]  New curves can
 * also be smoothed within the same BSpline object by calling solve() with
 * the new set of y values.  [See \Ref{BSpline}.]  A BSplineBase can be
 * created on its own, in which case all of the computations dependent on
 * the x values, boundary conditions, and cutoff wavelength have already
 * been completed.
 *
 * The solution of the cubic b-spline is divided into two parts.  The first
 * is the setup of the domain given the x values, boundary conditions, and
 * wavelength.  The second is the solution of the spline for a set of y
 * values corresponding to the x values in the domain.  The first part is
 * done in the creation of the BSplineBase object (or when calling the
 * setDomain method).  The second part is done when creating a BSpline
 * object (or calling solve() on a BSpline object).
 *
 * A BSpline object can be created with either one of its constructors, or
 * by calling apply() on an existing BSplineBase object.  Once a spline has
 * been solved, it can be evaluated at any x value.  The following example
 * creates a spline curve and evaluates it over the domain:
 *
\begin{verbatim}

    vector<float> x;
    vector<float> y;
    { ... }
    int bc = BSplineBase::BC_ZERO_SECOND;
    BSpline::Debug = true;
    BSpline spline (x.begin(), x.size(), y.begin(), wl, bc);
    if (spline.ok())
    {
        ostream_iterator<float> of(cout, "\t ");
    	float xi = spline.Xmin();
	float xs = (spline.Xmax() - xi) / 2000.0;
	for (; x <= spline.Xmax(); x += xs)
	{
	    *of++ = spline.evaluate (xi);
	}
    }
 
\end{verbatim}
 *
 * The algorithm is based on the cubic spline described by Katsuyuki Ooyama
 * in Montly Weather Review, Vol 115, October 1987 and on a previous
 * FORTRAN implementation.  The cubic b-spline is formulated as the sum of
 * some multiple of the basis function centered at each node in the domain.
 * The number of nodes is determined by the desired cutoff wavelength and a
 * desirable number of x values per node.  The basis function is continuous
 * and differentiable up to the second degree.  A derivative constraint is
 * included in the solution to achieve the effect of a low-pass frequency
 * filter with the given cutoff wavelength.  The domain nodes, boundary
 * constraints, and wavelength determine a linear system of equations,
 * Qa=b, where a is the vector of basis function coefficients at each node.
 * The coefficient vector is solved by first LU factoring along the
 * diagonally banded matrix Q in BSplineBase.  The BSpline object then
 * computes the B vector for a set of y values and solves for the
 * coefficient vector with the LU matrix.  Only the diagonal bands are
 * stored in memory and calculated during LU factoring and back
 * substitution, and the basis function is evaluated as few times as
 * possible in computing the diagonal matrix and B vector.
 *
 * {\small Interface version: $Id$}
 *
 * @author \URL[Gary Granger]{mailto:granger@atd.ucar.edu}
 * @see BSpline

 */
template <class T> 
class BSPLINE_DLL_ BSplineBase  
{
public:
    // Datum type
    typedef T datum_type;

    /// Return a string describing the implementation version.
    static const char *ImplVersion();

    /// Return a string describing the interface version.
    static const char *IfaceVersion();
	
    /**
     * Call this class method with a value greater than zero to enable
     * debug messages, or with zero to disable messages.  Calling with
     * no arguments returns true if debugging enabled, else false.
     */
    static bool Debug (int on = -1);

    /**
     * Boundary condition types.
     *
     * \begin{description}
     * \item[BC_ZERO_ENDPOINTS]	Set the endpoints of the spline to zero.
     * \item[BC_ZERO_FIRST]	Set the first derivative of the spline
     *				to zero at the endpoints.
     * \item[BC_ZERO_SECOND]	Set the second derivative to zero.
     * \end{description}
     */
    enum BoundaryConditionTypes
    {
	BC_ZERO_ENDPOINTS = 0,
	BC_ZERO_FIRST = 1,
	BC_ZERO_SECOND = 2
    };

public:

    /**
     * Construct a spline domain for the given set of x values, cutoff
     * wavelength, and boundary condition type.  The parameters are the
     * same as for setDomain().  Call ok() to check whether domain
     * setup succeeded after construction.
     *
     * @see setDomain 
     * @see ok
     */
    BSplineBase (const T *x, int nx, 
		 double wl, int bc_type = BC_ZERO_SECOND);

    /// Copy constructor
    BSplineBase (const BSplineBase &);

    /**
     * Change the domain of this base.  [If this is part of a BSpline
     * object, this method {\em will not} change the existing curve or
     * re-apply the smoothing to any set of y values.]
     *
     * The x values can be in any order, but they must be of sufficient
     * density to support the requested cutoff wavelength.  The setup of
     * the domain may fail because of either inconsistency between the x
     * density and the cutoff wavelength, or because the resulting matrix
     * could not be factored.  If setup fails, the method returns false.
     *
     * @param x		The array of x values in the domain.
     * @param nx	The number of values in the {\em x} array.
     * @param wl	The cutoff wavelength, in the same units as the
     *			{\em x} values.
     * @param bc_type	The enumerated boundary condition type.  If
     *			omitted it defaults to BC_ZERO_SECOND.
     *
     * @see ok
     */
    bool setDomain (const T *x, int nx, double wl, 
		    int bc_type = BC_ZERO_SECOND);

    /**
     * Create a BSpline smoothed curve for the given set of NX y values.
     * The returned object will need to be deleted by the caller.
     * @param y The array of y values corresponding to each of the nX()
     *		x values in the domain.
     * @see ok
     */
    BSpline<T> *apply (const T *y);

    /**
     * Return array of the node coordinates.  Returns 0 if not ok().  The
     * array of nodes returned by nodes() belongs to the object and should
     * not be deleted; it will also be invalid if the object is destroyed.
     */
    const T *nodes (int *nnodes);

    /** 
     * Return the number of nodes (one more than the number of intervals).
     */
    int nNodes () { return M+1; }

    /**
     * Number of original x values.
     */
    int nX () { return NX; }

    /// Minimum x value found.
    T Xmin () { return xmin; }

    /// Maximum x value found.
    T Xmax () { return xmin + (M * DX); }

    /** 
     * Return the Alpha value for a given wavelength.  Note that this
     * depends on the current node interval length (DX).
     */
    double Alpha (double wavelength);

    /**
     * Return alpha currently in use by this domain.
     */
    double Alpha () { return alpha; }

    /**
     * Return the current state of the object, either ok or not ok.
     * Use this method to test for valid state after construction or after
     * a call to setDomain().  ok() will return false if either fail, such
     * as when an appropriate number of nodes and node interval cannot be
     * found for a given wavelength, or when the linear equation for the
     * coefficients cannot be solved.
     */
    bool ok () { return OK; }

    virtual ~BSplineBase();

protected:

    typedef BSplineBaseP<T> Base;

    // Provided
    double waveLength;	// Cutoff wavelength (l sub c)
    int NX;
    int K;	// Degree of derivative constraint (currently fixed at 1)
    int BC;			// Boundary conditions type (0,1,2)

    // Derived
    T xmax;
    T xmin;
    int M;			// Number of intervals (M+1 nodes)
    double DX;			// Interval length in same units as X
    double alpha;
    bool OK;
    Base *base;			// Hide more complicated state members
    				// from the public interface.

    bool Setup ();
    void calculateQ ();
    double qDelta (int m1, int m2);
    double Beta (int m);
    void addP ();
    bool factor ();
    double Basis (int m, T x);
    double DBasis (int m, T x);

    static const double BoundaryConditions[3][4];
    static const double PI;

private:

    Ratio (int&, double &, double &, double *rd = 0);

};


template <class T> struct BSplineP;


/**
 * Inherit the BSplineBase domain information and interface and add
 * smoothing.  See the \Ref{BSplineBase} documentation for a summary of the
 * BSpline interface.
 *
 * @author \URL[Gary Granger]{mailto:granger@atd.ucar.edu}

 */
template <class T>
class BSPLINE_DLL_ BSpline : public BSplineBase<T>
{
public:
    /**
     * Create a single spline with the parameters required to set up
     * the domain and subsequently smooth the given set of y values.
     * The y values must correspond to each of the values in the x array.
     * If either the domain setup fails or the spline cannot be solved,
     * the state will be set to not ok.
     *
     * @see ok
     *
     * @param x		The array of x values in the domain.
     * @param nx	The number of values in the {\em x} array.
     * @param y		The array of y values corresponding to each of the
     *			nX() x values in the domain.
     * @param wl	The cutoff wavelength, in the same units as the
     *			{\em x} values.
     * @param bc_type	The enumerated boundary condition type.  If
     *			omitted it defaults to BC_ZERO_SECOND.
     */
    BSpline (const T *x, int nx, 		/* independent variable */
	     const T *y,			/* dependent values @ ea X */
	     double wl,				/* cutoff wavelength */
	     int bc_type = BC_ZERO_SECOND);

    /**
     * A BSpline curve can be derived from a separate Base and a set
     * of data points over that base.
     */
    BSpline (BSplineBase<T> &base, const T *y);

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the nX()
     *		x values in the domain.
     */
    bool solve (const T *y);

    /**
     * Return the entire curve evaluated at each of the nodes.
     * The array is held by the object, and thus should not be freed and
     * is only valid while the object exists.
     * If the current state is not ok(), the method returns zero.
     *
     * @param nx  If non-zero, returns the number of points in the curve.
     */
    const T *curve (int *nx = 0);

    /**
     * Return the evaluation of the smoothed curve 
     * at a particular x value.  If current state is not ok(), returns 0.
     */
    T evaluate (T x);

    /** 
     * Return the first derivative of the spline curve at the given x.
     * Returns zero if the current state is not ok().
     */
    T slope (T x);

    /**
     * Return the n-th basis coefficient, from 0 to M.  If the current
     * state is not ok(), or n is out of range, the method returns zero.
     */
    T coefficient (int n);

    virtual ~BSpline();

protected:

    // Our hidden state structure
    BSplineP<T> *s;
    T mean;			// Fit without mean and add it in later

};

#endif // !defined _BSPLINEBASE_IFACE_ID
