//
// BSpline.h: interface for the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _BSPLINEBASE_IFACE_ID
#define _BSPLINEBASE_IFACE_ID "$Id$"

class BSpline;

// Opaque structure to hide our matrix implementation, ala
// Cheshire cat.
struct BSplineBaseP;

class BSplineBase  
{
public:
	// Class members
	static const double PI;
	static bool Debug;						/* True enables debug messages */
	static const char *ImplVersion();		/* Return implementation version string */
	static const char *IfaceVersion();
	
	// Boundary condition types
	enum
	{
		BC_ZERO_ENDPOINTS = 0,			/* zero function value at endpoints */
		BC_ZERO_FIRST = 1,				/* zero first derivative */
		BC_ZERO_SECOND = 2				/* zero second derivative */
	};

public:

	BSplineBase (const float *x, int nx, 
				 float wl,							/* cutoff wavelength */
				 int bc_type = BC_ZERO_SECOND		/* boundary condition type */
				 );

	// Copy constructor
	BSplineBase (const BSplineBase &);

	// Set this spline to a whole new domain.  This does not yet work to
	// re-apply the smoothing to a BSpline curve.
	bool setDomain (const float *x, int nx, float wl, 
				    int bc_type = BC_ZERO_SECOND);

	// Create a BSpline smoothed curve for the given set of NX y values.  The
	// returned object will need to be deleted.
	BSpline *apply (const float *y);

	// These methods return info about the spline domain.  The array of nodes returned
	// by nodes() belongs to the object and should not be deleted; it will also be
	// invalid if the object is destroyed.
	const float *nodes (int *nnodes);
	int nNodes () { return M+1; }

	// These methods return information about the original X domain.
	int nX () { return NX; }
	float Xmin () { return xmin; }
	float Xmax () { return xmin + (M * DX); }

	// Return the Alpha value for a given wavelength.  With no argument, return
	// the object's current alpha.
	float Alpha (float wavelength);
	float Alpha () { return alpha; }

	// Use this method to test for valid state after construction or after a call to
	// setDomain().  ok() will return false if either fail, such as when an
	// appropriate number of nodes and node interval cannot be found for a given
	// wavelength, or when the linear equation for the coefficients
	// cannot be solved.
	bool ok () { return OK; }

	virtual ~BSplineBase();

protected:

	// Provided
	float waveLength;	// Cutoff wavelength (l sub c)
	int NX;
	int K;				// Degree of derivative constraint (currently fixed at 1)
	int BC;				// Boundary conditions type (0,1,2)

	// Derived
	float xmax;
	float xmin;
	int M;				// Number of intervals (M+1 nodes)
	float DX;			// Interval length in same units as X
	float alpha;
	bool OK;
	BSplineBaseP *base;	// Hide more complicated state members
						// from the public interface.

	//void Copy (const float *x, int nx, float wl);
	//void Reset ();
	bool Setup ();
	void calculateQ ();
	float qDelta (int m1, int m2);
	float Beta (int m);
	void addP ();
	bool factor ();
	float Basis (int m, float x);

	static const float BoundaryConditions[3][4];

private:

	Ratio (int&, float &, float &, float *rd = 0);

};


struct BSplineP;


class BSpline : public BSplineBase
{
public:
	// Return the entire curve evaluated at each of the nodes.
	// The array is held by the object, and thus should not be freed and
	// is only valid while the object exists.
	const float *curve (int *nx = 0);

	// Return the evaluation of the smoothed curve 
	// at a particular x value.
	float evaluate (float x);

	// Return the n-th basis coefficient, from 0 to M
	float coefficient (int n);

	// A BSpline curve must be derived from a Base and a set
	// of data points over that base.
	BSpline (BSplineBase &, const float *y);

	virtual ~BSpline();

protected:

	// Our hidden state structure
	BSplineP *s;
	float mean;			// Fit without mean and add it in later

};

#endif // !defined _BSPLINEBASE_IFACE_ID
