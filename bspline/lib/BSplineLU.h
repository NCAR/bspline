/* -*- mode: c++; c-basic-offset: 4; -*- */

/*
 * LU factor a diagonally banded matrix using Crout's algorithm, but
 * limiting the trailing sub-matrix multiplication to the non-zero
 * elements in the diagonal bands.  Return nonzero if a problem occurs.
 */
template <class Matrix>
int LU_factor_banded (Matrix &A, unsigned int bands)
{
    typename Matrix::size_type M = A.num_rows();
    typename Matrix::size_type N = A.num_cols();
    if (M != N)
	return 1;

    typename Matrix::size_type i,j,k;
    typename Matrix::element_type sum;

    for (j = 1; j <= N; ++j)
    {
	// Check for zero pivot
        if ( A(j,j) == 0 )                 
            return 1;

	// Calculate rows above and on diagonal. A(1,j) remains as A(1,j).
	for (i = (j > bands) ? j-bands : 1; i <= j; ++i)
	{	
	    sum = 0;
	    for (k = (j > bands) ? j-bands : 1; k < i; ++k)
	    {
		sum += A(i,k)*A(k,j);
	    }
	    A(i,j) -= sum;
	}

	// Calculate rows below the diagonal.
	for (i = j+1; (i <= M) && (i <= j+bands); ++i)
	{
	    sum = 0;
	    for (k = (i > bands) ? i-bands : 1; k < j; ++k)
	    {
		sum += A(i,k)*A(k,j);
	    }
	    A(i,j) = (A(i,j) - sum) / A(j,j);
	}
    }

    return 0;
}   



/*
 * Solving (LU)x = B.  First forward substitute to solve for y, Ly = B.
 * Then backwards substitute to find x, Ux = y.  Return nonzero if a
 * problem occurs.  Limit the substitution sums to the elements on the
 * bands above and below the diagonal.
 */
template <class Matrix, class Vector>
int LU_solve_banded(const Matrix &A, Vector &b, unsigned int bands)
{
    typename Matrix::size_type i,j;
    typename Matrix::size_type M = A.num_rows();
    typename Matrix::size_type N = A.num_cols();
    typename Matrix::element_type sum;

    if (M != N || M == 0)
	return 1;

    // Forward substitution to find y.  The diagonals of the lower
    // triangular matrix are taken to be 1.
    for (i = 2; i <= M; ++i)
    {
	sum = b[i-1];
	for (j = (i > bands) ? i-bands : 1; j < i; ++j)
	{
	    sum -= A(i,j)*b[j-1];
	}
	b[i-1] = sum;
    }

    // Now for the backward substitution
    b[M-1] /= A(M,M);
    for (i = M-1; i >= 1; --i)
    {
	if (A(i,i) == 0)	// oops!
	    return 1;
	sum = b[i-1];
	for (j = i+1; (j <= N) && (j <= i+bands); ++j)
	{
	    sum -= A(i,j)*b[j-1];
	}
	b[i-1] = sum / A(i,i);
    }

    return 0;
}

