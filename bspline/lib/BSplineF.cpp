/* -*- mode: c++; c-basic-offset: 4; -*- */
//
// $Id$
//
// Instantiate the BSpline templates for type float.  For WIN32, this
// module becomes a DLL.
//
//////////////////////////////////////////////////////////////////////
/*
 * Copyright (c) 1998,1999
 * University Corporation for Atmospheric Research, UCAR
 *
 * Permission to use, copy, modify, distribute and sell this software and
 * its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and that
 * both that copyright notice and this permission notice appear in
 * supporting documentation.  UCAR makes no representations about the
 * suitability of this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 * 
 * Note from the author:
 *
 * Where possible, you are encouraged to follow the GNU General Public
 * License, or at least the spirit of the license, for the distribution and
 * licensing of this software and any derived works.  See
 * http://www.gnu.org/copyleft/gpl.html.
 */

/*
 * We pre-empt the DLL macro definition in BSpline.cpp and BSpline.h
 * so that we can explicitly instantiate and export a template for
 * a specific type.
 */
#if WIN32
# define BSPLINE_DLL_ __declspec(dllexport)
#endif

#include "BSpline.cpp"

template class BSPLINE_DLL_ BSplineBase<float>;
template class BSPLINE_DLL_ BSpline<float>;

/* That's it! */
