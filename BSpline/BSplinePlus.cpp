#include "BSpline.h"
#include "BandedMatrix.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <map>
#include <assert.h>

//////////////////////////////////////////////////////////////////////
// BSplinePlus Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
/*
 * This BSpline constructor constructs and sets up a new base and 
 * solves for the spline curve coeffiecients all at once.
 */
template<class T> BSplinePlus<T>::BSplinePlus(const T *x,
                                              int nx,
                                              const T *y,
                                              double wl,
                                              int bc_type,
                                              int num_nodes,
                                              BLENDMODE blendMode,
                                              T blendSpan) :
    BSpline<T>(x, nx, y, wl, bc_type, num_nodes), _blendMode(blendMode),
            _blendSpan(blendSpan), _blendXleftIndex(0), _blendXrightIndex(nx-1) {

    // verify that x values are in ascending order.
    for (int i = 1; i < nx; i++) {
        assert(x[i-1] < x[i]);
    }

    // determine the x indices of the data points that lie outside of the 
    // blending range. This will be used to select the points that the
    // blending will be applied to. They are computed and saved now, 
    // for convenience later on.
    for (int i = 1; i < nx; i++) {
        if (_blendMode == BSplinePlus<T>::BLENDSTART || _blendMode
                == BSplinePlus<T>::BLENDBOTH) {
            if ((x[i-1]-BSplineBase<T>::xmin) < _blendSpan)
                _blendXleftIndex = i;
        }
        if (_blendMode == BSplinePlus<T>::BLENDFINISH || _blendMode
                == BSplinePlus<T>::BLENDBOTH) {
            if ((x[i-1] + _blendSpan) < x[nx-1])
                _blendXrightIndex = i;
        }
    }

    // Create the blending coefficients. Once again computed here, 
    // for convenience later on. The blend vector is initialized to
    // 1.0, and the ends will be adjusted according to the blend mode.
    _blendWeights.assign(nx, 1.0);

    // The left end:
    for (unsigned int i = 0; i < _blendXleftIndex; i++) {
        double deltax = (x[i]-BSplineBase<T>::xmin);
        double weight = deltax/_blendSpan;
        _blendWeights[i] = weight;
    }
    // The right end:
    for (unsigned int i = _blendXrightIndex; i < nx; i++) {
        double deltax = (BSplineBase<T>::xmax - x[i]);
        double weight = deltax/_blendSpan;
        _blendWeights[i] = weight;
    }

    if (BSplineBase<T>::Debug()) {
        // save all of our diagnostic prints for here.

        std::map<BLENDMODE, std::string > modeNames;
        modeNames[BLENDNONE] = "none";
        modeNames[BLENDSTART] = "start";
        modeNames[BLENDFINISH] = "finish";
        modeNames[BLENDBOTH] = "both";

        std::cout << "Blending mode: " << modeNames[_blendMode] << "\n";
        std::cout << "Blending length: " << _blendSpan
                << "(in same units as x)\n";
        std::cout << "xmin:" << BSplineBase<T>::xmin << ", xmax:" << BSplineBase<T>::xmax <<"\n";
        std::cout << "\n";
        for (unsigned int i = 0; i < _blendWeights.size(); i++) {
            if (_blendWeights[i] != 1.0) {
                std::cout << "Spline weight[" << i << "] for x("
                        << std::setprecision(4) << x[i] << "):"
                        << std::setprecision(2) << _blendWeights[i] << "\n";
            }
        }
        std::cout << "\n";
    }

}

//////////////////////////////////////////////////////////////////////
/*
 * Create a new spline given a BSplineBase.
 */
template<class T> BSplinePlus<T>::BSplinePlus(BSplineBase<T> &bb,
                                              const T *y) :
    BSpline<T>(bb, y) {
}
//////////////////////////////////////////////////////////////////////
/*
 * (Re)calculate the spline for the given set of y values.
 */
template<class T> bool BSplinePlus<T>::solve(const T *y) {
    bool ok = BSpline<T>::solve(y);

    if (!ok)
        return false;

    return (true);
}
//////////////////////////////////////////////////////////////////////
template<class T> BSplinePlus<T>::~BSplinePlus() {
}
//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::evaluate(T x) {
    T y = BSpline<T>::evaluate(x);

    return y;
}
//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::slope(T x) {

    T dy = BSpline<T>::slope(x);

    return dy;
}
//////////////////////////////////////////////////////////////////////
template<class T> const T *BSplinePlus<T>::curve(int *nx) {
    const T* c = BSpline<T>::curve(nx);

    if (!c)
        return 0;

    return c;
}


