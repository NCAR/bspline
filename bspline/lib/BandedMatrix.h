/*
 * Template for a diagonally banded matrix.
 */

#ifndef _BANDEDMATRIX_ID
#define _BANDEDMATRIX_ID "$Id$"

#include <iostream.h>
#include <vector>


template <class T> class BandedMatrixRow;


template <class T> class BandedMatrix
{
public:
	typedef int size_type;
	typedef T element_type;

	// Create a banded matrix with the same number of bands above and below the
	// diagonal.
	BandedMatrix (int _N = 1, int nbands_off_diagonal = 0) : bands(0)
	{
		if (! setup (_N, nbands_off_diagonal))
			setup ();
	}

	// Create a banded matrix by naming the first and last non-zero bands, where
	// the diagonal is at zero, and bands below the diagonal are negative, bands
	// above the diagonal are positive.
	BandedMatrix (int _N, int first, int last) : bands(0)
	{
		if (! setup (_N, first, last))
			setup ();
	}

	// Copy constructor
	BandedMatrix (const BandedMatrix &b) : bands(0)
	{
		Copy (*this, b);
	}

	inline bool setup (int _N = 1, int noff = 0)
	{
		return setup (_N, -noff, noff);
	}

	bool setup (int _N, int first, int last)
	{
		// Check our limits first and make sure they make sense.
		// Don't change anything until we know it will work.
		if (first > last || _N <= 0)
			return false;

		// Need at least as many _N as columns and as rows in the bands.
		if (_N < abs(first) || _N < abs(last))
			return false;

		top = last;
		bot = first;
		N = _N;
		out_of_bounds = T();

		// Finally setup the diagonal vectors
		nbands = last - first + 1;
		if (bands) delete[] bands;
		bands = new std::vector<T>[nbands];
		int i;
		for (i = 0; i < nbands; ++i)
		{
			// The length of each array varies with its distance from the diagonal
			int len = N - (abs(bot + i));
			bands[i].assign (len);
		}
		return true;
	}

	BandedMatrix<T> & operator= (const BandedMatrix<T> &b) 
	{
		return Copy (*this, b);
	}

	BandedMatrix<T> & operator= (const T &e)
	{
		int i;
		for (i = 0; i < nbands; ++i)
		{
			bands[i].assign (bands[i].size(), e);
		}
		out_of_bounds = e;
		return (*this);
	}

	~BandedMatrix ()
	{
		if (bands)
			delete[] bands;
	}

private:
	// Return false if coordinates are out of bounds
	inline bool check_bounds (int i, int j, int &v, int &e) const
	{
		v = (j - i) - bot;
		e = (i >= j) ? j : i;
		return !(v < 0 || v >= nbands || e < 0 || e >= bands[v].size());
	}

	static BandedMatrix & Copy (BandedMatrix &a, const BandedMatrix &b)
	{
		if (a.bands) delete[] a.bands;
		a.top = b.top;
		a.bot = b.bot;
		a.N = b.N;
		a.out_of_bounds = b.out_of_bounds;
		a.nbands = a.top - a.bot + 1;
		a.bands = new std::vector<T>[a.nbands];
		int i;
		for (i = 0; i < a.nbands; ++i)
		{
			a.bands[i] = b.bands[i];
		}
		return a;
	}

public:
	T &element (int i, int j)
	{
		int v, e;
		if (check_bounds(i, j, v, e))
			return (bands[v][e]);
		else
			return out_of_bounds;
	}

	const T &element (int i, int j) const
	{
		int v, e;
		if (check_bounds(i, j, v, e))
			return (bands[v][e]);
		else
			return out_of_bounds;
	}

	inline T & operator() (int i, int j) 
	{
		return element (i-1,j-1);
	}

	inline const T & operator() (int i, int j) const
	{
		return element (i-1,j-1);
	}

	size_type num_rows() const { return N; }

	size_type num_cols() const { return N; }

	const BandedMatrixRow<T> operator[] (int row) const
	{
		return BandedMatrixRow<T>(*this, row);
	}

	BandedMatrixRow<T> operator[] (int row)
	{
		return BandedMatrixRow<T>(*this, row);
	}


private:

	//typedef std::vector<T> band_type;	// These parse under MSVC, not tested
	//std::vector<band_type> b;

	int top;
	int bot;
	int nbands;
	std::vector<T> *bands;
	//std::vector<vector<T>> bands;		// MSVC++ won't take this
	int N;
	T out_of_bounds;

};


template <class T>
ostream &operator<< (ostream &out, const BandedMatrix<T> &m)
{
	int i, j;
	for (i = 0; i < m.num_rows(); ++i)
	{
		for (j = 0; j < m.num_cols(); ++j)
		{
			out << m.element (i, j) << " ";
		}
		out << endl;
	}
	return out;
}



/*
 * Helper class for the intermediate in the [][] operation.
 */
template <class T> class BandedMatrixRow
{
public:
	BandedMatrixRow (const BandedMatrix<T> &_m, int _row) : bm(_m), i(_row)
	{ }

	BandedMatrixRow (BandedMatrix<T> &_m, int _row) : bm(_m), i(_row)
	{ }

	~BandedMatrixRow () {}

	BandedMatrix<T>::element_type & operator[] (int j)
	{
		return const_cast<BandedMatrix<T> &>(bm).element (i, j);
	}

	const BandedMatrix<T>::element_type & operator[] (int j) const
	{
		return bm.element (i, j);
	}

private:
	const BandedMatrix<T> &bm;
	int i;
};



#endif /* _BANDEDMATRIX_ID */

