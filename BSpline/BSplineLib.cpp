/************************************************************************\
	Copyright 2008 University Corporation for Atmospheric Research.
	All rights reserved.
	Use of this code is subject to UCAR's standard Terms of Use,
	which can be found at http://www.ucar.edu/legal/terms_of_use.shtml .
	By using this source code, you agree to abide by those Terms of Use.
\*************************************************************************/

// Instantiate the BSpline templates for type 

#include "BSplineBase.cpp"
#include "BSpline.cpp"
#include "BSplinePlus.cpp"

/// Instantiate BSplineBase for a library
template class  BSplineBase<double>;
template class  BSplineBase<float>;

/// Instantiate BSpline for a library
template class BSpline<double>;
template class BSpline<float>;

/// Instantiate BSplinePlus for a library
template class BSplinePlus<double>;
template class BSplinePlus<float>;
