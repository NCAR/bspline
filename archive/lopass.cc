// $Id$
//
// Converted from a FORTRAN original by someone else
//
// "Complete symmetric filter with multiplicative adjustment"
//

#include <stdlib.h>
#include <math.h>
#include <iostream.h>

// FORTRAN version
extern "C" {
void flopass_ (float *gin, float *gout, int *nmax, float *frac, int *ntrm);
}

/*
 * The input values are in 'gin', the filtered values are returned in
 * 'gout'.  Both arrays have 'nmax' elements.  'ntrm' is the maximum radius
 * of influence from surrounding points, reduced at the edges and without
 * wraparound to the opposite edge.  'ntrm' must be greater than 0.  'frac'
 * is the multiplicative adjustment.
 */
void
lopass (float *gin, float *gout, int nmax, float frac, int ntrm)
{
	const int mterms = 2*ntrm+2;
	float *wt = new float[mterms];		// Use 1 as base index
	float omcut = M_PI * frac;
	int n;

	// ---- "Endpoints" ----
	gout[0] = gin[0];
	gout[nmax - 1] = gin[nmax - 1];

	/*
	 * Loop over each interior point in the input array.
	 */
	for (int l = 1; l < nmax - 1; ++l)
	{
		// mtrm is the number of terms to include on either side
		// of the current point.  It equals ntrm if there are ntrm
		// points on each side of the current point, otherwise it
		// is reduced to maximum number of points which fit on both
		// sides of the point.  So points near the edge have a 
		// smaller weighting window.
		int mtrm = ntrm;

		if (l < mtrm)
			mtrm = l;
		if (l >= nmax - mtrm)
			mtrm = nmax - l - 1;

		// ---- "Calculate weights" ----
		float wfsum = 0.0;
		for (n = 1; n <= mtrm; ++n)
		{
			float fac = 2.0 * M_PI * n / (2.0 * mtrm + 1);
			wt[n] = sin(omcut * n)/(M_PI * n)*sin(fac)/fac;
			// Sum the weight for both(2) sides
			wfsum += 2.0 * wt[n];
		}

		// Add the weight for the center point
		wfsum += frac;
		float wto = frac / wfsum;	// wt0, or "weight not"

		// Normalize the weights
		for (n = 1; n <= mtrm; ++n)
		{
			wt[n] /= wfsum;
		}

		// Finally sum the weighted surrounding points for the output
		float sum = 0.0;
		for (n = 1; n <= mtrm; ++n)
		{
			sum += (gin[l-n] + gin[l+n]) * wt[n];
		}
		gout[l] = sum + gin[l] * wto;
	}

	delete[] wt;
}



int main (int argc, char *argv[])
{
	// Build some test arrays and run them through the filter

	if (argc > 3)
	{
		cerr << "Usage: " << argv[0] << " [fraction] [nterm]\n";
		exit (1);
	}

	float frac = 0.5;
	int nterm = 20;

	if (argc == 3)
		frac = atof(argv[2]);
	if (argc >= 2)
		nterm = atoi(argv[1]);

	const int N = 5000;

	float *x = new float[N];
	float *y = new float[N];
	float *y2 = new float[N];

	int n = 0;
	while (cin >> x[n])
	{
		// cout << x[n] << endl;
		++n;
	}

	lopass (x, y, n, frac, nterm);
	flopass_ (x, y2, &n, &frac, &nterm);
	for (int i = 0; i < n; ++i)
	{
		
		cout << " " << y[i] << " " << y2[i];
		cout << " " << y2[i] - y[i];
		if (fabs(y[i] - y2[i]) > 0.00001)
			cout << " *** ";
		cout << endl;
	}

	delete[] x;
	delete[] y;
	delete[] y2;
}

