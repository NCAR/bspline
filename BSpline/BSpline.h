/* -*- mode: c++; c-basic-offset: 4; -*- */
//
// BSpline.h: interface for the BSplineBase class.
//
//////////////////////////////////////////////////////////////////////
/*
 * Copyright (c) 1998-2008
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

#ifndef BSPLINE_H
#define BSPLINE_H

#include "BSplineBase.h"
#include <vector>

template <class T> struct BSplineP;


/**
 * Used to evaluate a BSpline.
 * Inherits the BSplineBase domain information and interface and adds
 * smoothing.  See the BSplineBase documentation for a summary of the
 * BSpline interface.
 */
template <class T>
class BSpline : public BSplineBase<T>
{
public:
    /**
     * Create a single spline with the parameters required to set up
     * the domain and subsequently smooth the given set of y values.
     * The y values must correspond to each of the values in the x array.
     * If either the domain setup fails or the spline cannot be solved,
     * the state will be set to not ok.
     *
     * @see ok().
     *
     * @param x		The array of x values in the domain.
     * @param nx	The number of values in the @p x array.
     * @param y		The array of y values corresponding to each of the
     *			nX() x values in the domain.
     * @param wl	The cutoff wavelength, in the same units as the
     *			@p x values.  A wavelength of zero disables
     *			the derivative constraint.
     * @param bc_type	The enumerated boundary condition type.  If
     *			omitted it defaults to BC_ZERO_SECOND.
     * @param num_nodes The number of nodes to use for the cubic b-spline.
     *			If less than 2 a "reasonable" number will be
     *			calculated automatically, taking into account
     *			the given cutoff wavelength.
     */
    BSpline (const T *x, int nx, 		/* independent variable */
	     const T *y,			/* dependent values @ ea X */
	     double wl,				/* cutoff wavelength */
	     int bc_type = BSplineBase<T>::BC_ZERO_SECOND,
	     int num_nodes = 0);

    /**
     * A BSpline curve can be derived from a separate @p base and a set
     * of data points @p y over that base.
     */
    BSpline (BSplineBase<T> &base, const T *y);

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the nX()
     *		x values in the domain.
     */
    bool solve (const T *y);

    /**
     * Return the entire curve evaluated at each of the nodes.
     * The array is held by the object, and thus should not be freed and
     * is only valid while the object exists.
     * If the current state is not ok(), the method returns zero.
     *
     * @param nx  If non-zero, returns the number of points in the curve.
     */
    const T *curve (int *nx = 0);

    /**
     * Return the evaluation of the smoothed curve 
     * at a particular @p x value.  If current state is not ok(), returns 0.
     */
    T evaluate (T x);

    /** 
     * Return the first derivative of the spline curve at the given @p x.
     * Returns zero if the current state is not ok().
     */
    T slope (T x);

    /**
     * Return the @p n-th basis coefficient, from 0 to M.  If the current
     * state is not ok(), or @p n is out of range, the method returns zero.
     */
    T coefficient (int n);

    virtual ~BSpline();

    using BSplineBase<T>::Debug;

protected:

    using BSplineBase<T>::OK;
    using BSplineBase<T>::M;
    using BSplineBase<T>::NX;
    using BSplineBase<T>::DX;
    using BSplineBase<T>::base;
    using BSplineBase<T>::xmin;
    using BSplineBase<T>::xmax;

    // Our hidden state structure
    BSplineP<T> *s;
    T mean;			// Fit without mean and add it in later

};

#endif 
