/*
 * Template for a diagonally banded matrix.
 */

#ifndef _BANDEDMATRIX_ID
#define _BANDEDMATRIX_ID "$Id$"


#include <vector>


template <class T> class BandedMatrixRow;


template <class T> class BandedMatrix
{
public:
	typedef int size_type;
	typedef T element_type;

	// Create a banded matrix with the same number of bands above and below the
	// diagonal.
	BandedMatrix (int _N = 1, int nbands_off_diagonal = 0)
	{
		setup (_N, nbands_off_diagonal);
	}

	// Create a banded matrix by naming the first and last non-zero bands, where
	// the diagonal is at zero, and bands below the diagonal are negative, bands
	// above the diagonal are positive.
	BandedMatrix (int _N, int first, int last)
	{
		setup (_N, first, last);
	}

	inline bool setup (int _N, int noff)
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
		int nbands = last - first + 1;
		if (bands) delete[] bands;
		bands = new vector<T>[nbands];
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
		if (bands) delete[] bands;
		top = b.top;
		bot = b.bot;
		N = b.N;
		out_of_bounds = b.out_of_bounds;
		int nbands = top - bot + 1;
		bands = new vector<T>[top - bot + 1];
		int i;
		for (i = 0; i < nbands; ++i)
		{
			bands[i] = b.bands[i];
		}
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
		return !(v < 0 || v >= bands.size() || e < 0 || e >= bands[v].size());
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

	size_type num_rows() const { return N; }

	size_type num_cols() const { return N; }

	BandedMatrixRow<T> operator[] (int row) const
	{
		return BandedMatrixRow<T>(*this, row);
	}


private:

	int top;
	int bot;
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
	BandedMatrixRow (BandedMatrix<T> &_m, int _row) : bm(_m), i(_row)
	{ }

	BandedMatrix<T>::element_type & operator[] (int j)
	{
		return bm.element (i, j);
	}

private:
	BandedMatrix<T> bm;
	int i;
};


#endif /* _BANDEDMATRIX_ID */

