/* -*- mode: c++; c-basic-offset: 4; -*- */
//
// $Id$
//
// Instantiate the BSpline templates for type float.  For WIN32, this
// module becomes a DLL.
//
//////////////////////////////////////////////////////////////////////

/*
 * We pre-empt the DLL macro definition in BSpline.cxx and BSpline.h
 * so that we can explicitly instantiate and export a template for
 * a specific type.
 */
#if WIN32
# define BSPLINE_DLL_ __declspec(dllexport)
#endif

#include "BSpline.cxx"

template class BSPLINE_DLL_ BSplineBase<float>;
template class BSPLINE_DLL_ BSpline<float>;

/* That's it! */
