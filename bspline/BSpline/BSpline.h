//
// $Id$
//
// BSpline.h: interface for the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BSPLINE_H__97C5C270_06EE_11D2_BF96_00C04FA30D33__INCLUDED_)
#define AFX_BSPLINE_H__97C5C270_06EE_11D2_BF96_00C04FA30D33__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

class BSpline;

// Opaque structure to hide our matrix implementation, ala
// Cheshire cat.
struct BSplineBaseP;

class BSplineBase  
{
public:
	// Class members
	static const double PI;

public:

	BSplineBase (const float *x, int nx, float wl /*cutoff wavelength*/);
	/* BSplineBase (); */

	// The copy constructor, which is especially used by the
	// BSpline subclass constructor.
	BSplineBase (const BSplineBase &);

	virtual ~BSplineBase();

	int nX () { return NX; }
	void setDomain (const float *x, int nx, float wl);

	BSpline *apply (const float *y);
	const float *nodes (int *nnodes);
	float Xmin () { return xmin; }
	float Xmax () { return xmin + (M * DX); }
	float Alpha (float wavelength);

protected:

	// Provided
	float waveLength;	// Cutoff wavelength (l sub c)
	int NX;
	int K;				// Degree of derivative constraint
	int BC;				// Boundary conditions type (0,1,2)

	// Derived
	float xmax;
	float xmin;
	int M;				// Number of intervals (M+1 nodes)
	float DX;			// Interval length in same units as X
	float alpha;
	BSplineBaseP *base;	// Hide more complicated state members
						// from the public interface.

	void Copy (const float *x, int nx, float wl);
	void Reset ();
	int Setup ();
	void calculateQ ();
	float qDelta (int m1, int m2);
	float Beta (int m);
	void addP ();
	void factor ();
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
	float mean;			// Fit without mean but add it in later

};

#endif // !defined(AFX_BSPLINE_H__97C5C270_06EE_11D2_BF96_00C04FA30D33__INCLUDED_)
