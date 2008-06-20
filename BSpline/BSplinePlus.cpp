/************************************************************************\
	Copyright 2008 University Corporation for Atmospheric Research.
	All rights reserved.
	Use of this code is subject to UCAR's standard Terms of Use,
	which can be found at http://www.ucar.edu/legal/terms_of_use.shtml .
	By using this source code, you agree to abide by those Terms of Use.
\*************************************************************************/

#include "BSplinePlus.h"
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
            _blendSpan(blendSpan)
{
    // setup the blending artifacts.
    initBlending(y);

}

//////////////////////////////////////////////////////////////////////
/*
 * Create a new spline given a BSplineBase.
 */
template<class T> BSplinePlus<T>::BSplinePlus(BSplineBase<T> &bb,
                                              const T *y,
                                              BLENDMODE blendmode,
                                              T blendspan) :
    BSpline<T>(bb, y), _blendMode(blendmode), _blendSpan(blendspan)
{
    if (!OK)
        return;
    
    // setup the blending artifacts.
    initBlending(y);
}

//////////////////////////////////////////////////////////////////////
template<class T> void BSplinePlus<T>::initBlending(const T* y)
{
    // verify that the x values are ordered
    for (unsigned int i = 1; i < NX; i++) {
        assert(base->X[i-1] < base->X[i]);
    }

    // save y
    _y.resize(NX);
    for (unsigned int i = 0; i < NX; i++)
        _y[i] = y[i];

    // save the boundaries of our blending areas
    _xLeft = xmin + _blendSpan;
    _xRight = xmax - _blendSpan;

    if (BSplineBase<T>::Debug()) {
        // save all of our diagnostic prints for here.

        std::map<BLENDMODE, std::string > modeNames;
        modeNames[BLENDNONE] = "none";
        modeNames[BLENDSTART] = "start";
        modeNames[BLENDFINISH] = "finish";
        modeNames[BLENDBOTH] = "both";

        std::cout << "Blending mode: " << modeNames[_blendMode].c_str() << "\n";
        std::cout << "Blending span: " << _blendSpan
                << "(in same units as x)\n";
        std::cout << "Blend start from " << xmin << " to " << _xLeft <<"\n";
        std::cout << "Blend finish from " << _xRight << " to " << xmax << "\n";
    }
}

//////////////////////////////////////////////////////////////////////
template<class T> BSplinePlus<T>::~BSplinePlus()
{
}

//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::evaluate(T x)
{
    T y = BSpline<T>::evaluate(x);

    // Perform the blending. Two interpolations are involved:
    //   The value of the original Y series, at the desired x value.
    //   The ratio of the x position in the span provides a 
    //     weighting foactor for the blending. As the x value approaches
    //     the end points of the series, the weighting goes to 100% in
    //     favor of the original data.
    
    
    if (x < _xLeft && x >= xmin) {
        // blend the left (start) side of the series
        
        // get the value of the original y,
        // interpolated at X
        T originalY = interpLeft(x);
        // calculate out the blending ratio
        T factor = (x-xmin)/_blendSpan;
        // blend the spline value and the original value
        T newY = factor*y + (1.0-factor)*originalY;
        y = newY;
    } else if (x > _xRight && x <= xmax) {
        // blend the right (Finish) side of the series
        
        // get the value of the original y,
        // interpolated at X
        T originalY = interpRight(x);
        // calculate out the blending ratio
        T factor = (xmax-x)/_blendSpan;
        // blend the spline value and the original value
        T newY = factor*y + (1.0-factor)*originalY;
        y = newY;
    }

    return y;
}

//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::slope(T x)
{

    T dy = BSpline<T>::slope(x);

    return dy;
}

//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::interpLeft(T x)
{

    T result = 0.0;
    if (x >= xmin) {
        for (int i = 0; i < NX-2; i++) {
            // look for the X interval containing the requested x
            if (x <= base->X[i+1]) {
                T deltaY = (_y[i+1]-_y[i]);
                T factor = (x - base->X[i])/(base->X[i+1]-base->X[i]);
                T increment = factor*deltaY;
                result = _y[i] + increment;
                break;
            }
        }
    }

    return result;
}

//////////////////////////////////////////////////////////////////////
template<class T> T BSplinePlus<T>::interpRight(T x)
{

    T result = 0.0;
    if (x <= xmax) {
        for (int i = NX-1; i > 0; i--) {
            // look for the X interval containing the requested x
            if (x >= base->X[i-1]) {
                T deltaY = (_y[i]-_y[i-1]);
                T factor = (base->X[i]-x)/(base->X[i]-base->X[i-1]);
                T increment = factor*deltaY;
                result = _y[i] - increment;
                break;
            }
        }
    }

    return result;

}
