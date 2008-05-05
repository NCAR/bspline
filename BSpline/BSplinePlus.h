#ifndef BSPLINEPLUS_H_
#define BSPLINEPLUS_H_

#include "BSpline.h"

/**
 * Add modified behavior to the standard BSpline.
 * Inherits the BSpline. Allows the BSpline boundary 
 * conditions to be overriden with a blending scheme.
 * See the BSpline documentation for a summary of the
 * BSpline interface. Note that for BSplinePluse, @em the 
 * x values must be in ascending order.
 */
template <class T>
class BSplinePlus : public BSpline<T>
{
public:
    /// Blending modes
    enum BLENDMODE {BLENDNONE, BLENDSTART, BLENDFINISH, BLENDBOTH};
    /**
     * Create a single spline with the parameters required to set up
     * the domain and subsequently smooth the given set of y values.
     * The y values must correspond to each of the values in the x array.
     * If either the domain setup fails or the spline cannot be solved,
     * the state will be set to not ok.
     * 
     * The blending pan must be smaller than half the total span of
     * the spline.
     *
     * @see ok().
     *
     * @param x     The array of x values in the domain.
     * @param nx    The number of values in the @p x array.
     * @param y     The array of y values corresponding to each of the
     *          nX() x values in the domain.
     * @param wl    The cutoff wavelength, in the same units as the
     *          @p x values.  A wavelength of zero disables
     *          the derivative constraint.
     * @param bc_type   The enumerated boundary condition type.  If
     *          omitted it defaults to BC_ZERO_SECOND.
     * @param num_nodes The number of nodes to use for the cubic b-spline.
     *          If less than 2 a "reasonable" number will be
     *          calculated automatically, taking into account
     *          the given cutoff wavelength.
     * @param blendmode Choose blending at the endpoints.
     * @param blendspan The distance over which the blending is applied.
     */
    BSplinePlus (const T *x, 
                 int nx,        
                 const T *y,            
                 double wl,             
                 int bc_type = BSplineBase<T>::BC_ZERO_SECOND,
                 int num_nodes = 0,
                 BLENDMODE blendmode = BSplinePlus<T>::BLENDNONE,
                 T blendspan = 0.0);

    /**
     * A BSpline curve can be derived from a separate @p base and a set
     * of data points @p y over that base.
     */
    BSplinePlus (BSplineBase<T> &base, const T *y);

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the nX()
     *      x values in the domain.
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

    virtual ~BSplinePlus();

protected:
    /// The blending mode for the end points
    BLENDMODE _blendMode;
    /// The span over which the blending occurs
    T _blendSpan;
    /// The index of the first x value that is past the start (left side)
    /// blending interval. Will be undefined for blend mode == BLENDNONE os
    /// blend mode == BLENDFINISH
    int _blendXleftIndex;
    /// The index of the first x value that is before the finish (right side)
    /// blending interval. Will be undefined for blend mode == BLENDNONE os
    /// blend mode == BLENDSTART
    int _blendXrightIndex;
    /// weighting factors to be applied to the points in the series
    /// Weighting of one preserves the spline calculation, a zero 
    /// gives full weight to the input data
    std::vector<double> _blendWeights;    
};

#endif /*BSPLINEPLUS_H_*/
