/* -*- mode: c++; c-basic-offset: 4; -*- */

/*
 * This is a modified version of R. Pozo's LU_factor template procedure from
 * the Template Numerical Toolkit.  It was modified to limit pivot searching
 * in the case of banded diagonal matrices.  The extra parameter BANDS
 * is the number of bands below the diagonal.  Also, this routine now does
 * NO pivoting, so the index vector has been eliminated.
 */
template <class Matrix>
int LU_factor_banded (Matrix &A, int bands)
{
    typename Matrix::size_type M = A.num_rows();
    typename Matrix::size_type N = A.num_cols();

    typename Matrix::size_type j,k,jp;
    typename Matrix::size_type minMN = min(M,N);

    for (j=1; j<= minMN; j++)
    {
        jp = j;

        // jp now has the index of maximum element 
        // of column j, below the diagonal

        if ( A(jp,j) == 0 )                 
            return 1;       // factorization failed because of zero pivot

        if (j<M)                // compute elements j+1:M of jth column
        {
            // note A(j,j), was A(jp,p) previously which was
            // guarranteed not to be zero (Label #1)
            //
            typename Matrix::element_type recp =  1.0 / A(j,j);

            for (k=j+1; (k <= j+bands) && (k<=M); k++)
                A(k,j) *= recp;
        }

        if (j < minMN)
        {
            // rank-1 update to trailing submatrix:   E = E - x*y;
            //
            // E is the region A(j+1:M, j+1:N)
            // x is the column vector A(j+1:M,j)
            // y is row vector A(j,j+1:N)

	    typename Matrix::size_type ii,jj;

            for (ii=j+1; (ii <= j+bands) && (ii<=M); ii++)
                for (jj=j+1; (jj <= j+bands) && (jj<=N); jj++)
                    A(ii,jj) -= A(ii,j)*A(j,jj);
        }
    }

    return 0;
}   



template <class Matrix, class Vector>
int LU_solve_banded(const Matrix &A, Vector &b)
{
    typename Matrix::size_type i,ii=0,ip,j;
    typename Matrix::size_type n = A.num_rows();
    typename Matrix::element_type sum = 0.0;

    for (i=1;i<=n;i++) 
    {
        ip=i;
        sum=b[ip-1];
        b[ip-1]=b[i-1];
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

