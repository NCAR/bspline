
#ifndef _BSPLINESOLVER_ID
#define _BSPLINESOLVER_ID "$Id$"

#include <vector>

template <class Matrix> class BSplineSolver
{
public:

	typedef Matrix::element_type element_type;
	typedef Matrix::size_type size_type;
	typedef std::vector<element_type> vector_type;

	BSplineSolver () : mp(0), N(0), ok(false)
	{
		// We are invalid until upper() is called with a reference matrix!
	}

	// Solve vector a for M*a = b
	bool solve (const vector_type &b, vector_type &a)
	{
		if (!ok)
			return false;

		// Traverse backwards to solve for each a, calculating fm along the way
		a.assign (N);
		size_type i = N;
		Matrix &M = *mp;
#if 0
		if (--i >= 0)
			a[i] = b[i] / M[i][i];
		if (--i >= 0)
			a[i] = (e1[i]*a[i+1]) + (b[i]/M[i][i]);
		if (--i >= 0)
			a[i] = (e1[i]*a[i+1]) + (e2[i]*a[i+2]) + (b[i]/M[i][i]);
		while (--i >= 0)
			a[i] = (e1[i]*a[i+1]) + (e2[i]*a[i+2]) + (e3[i]*a[i+3]) + (b[i]/M[i][i]);
#endif
		if (--i >= 0)
			a[i] = b[i] / M[i][i];
		if (--i >= 0)
			a[i] = (b[i] - M[i][i+1]*a[i+1]) / M[i][i];
		if (--i >= 0)
			a[i] = (b[i] - M[i][i+1]*a[i+1] - M[i][i+2]*a[i+2]) / M[i][i];
		while (--i >= 0)
			a[i] = (b[i] - M[i][i+1]*a[i+1] - M[i][i+2]*a[i+2] - M[i][i+3]*a[i+3]) / M[i][i];

		return true;	
	}

	bool OK() { return ok; }

	const Matrix *matrix () { return mp; }

	// Upper triangularize and generate our auxiliary e vectors
	bool upper (Matrix &_M)
	{
		mp = &_M;
		Matrix &M = *mp;
		N = M.num_rows();
		
		// Use Gaussian elimination to convert to upper diagonal matrix.  NOTE this
		// does no pivoting, relying heavily on the assumption of a spline matrix
		// having exactly 7 bands ceneterd on a dominant diagonal.
		e1.assign(N);
		e2.assign(N);
		e3.assign(N);

		// For each diagonal element, zero the elements in its column on the three
		// rows beneath it.  We only need to update to the third column right of the 
		// diagonal on the row we're multiplying/adding; everything right of that 
		// in that row is zero.
		size_type i, j;
		for (j = 0; j < N; ++j)
		{
			for (i = j+1; (i <= j+3) && (i < N); ++i)
			{
				if (M[j][j] == 0)
					return (ok = false);
				element_type r = - (M[i][j] / M[j][j]);
				M[i][j] = 0;
				for (int jj = j+1; (jj <= j + 3) && (jj < N); ++jj)
				{
					M[i][jj] += M[j][jj] * r;
				}
			}
					
			// Now calculate the e vector elements for this row.
		}

		return (ok = true);
	}

private:

	Matrix *mp;
	Matrix::size_type N;
	bool ok;

	vector_type e1, e2, e3;

};



#endif /* _BSPLINESOLVER_ID */
