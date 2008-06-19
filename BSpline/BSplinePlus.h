#ifndef BSPLINEPLUS_H_
#define BSPLINEPLUS_H_

#include "BSpline.h"

/**
 * Add modified behavior to the standard BSpline.
 * Inherits the BSpline. Allows the BSpline boundary 
 * conditions to be overriden with a blending scheme.
 * See the BSpline documentation for a summary of the
 * BSpline interface. 
 */
template<class T> class BSplinePlus : public BSpline<T> {
    public:
        /// Blending modes
        enum BLENDMODE {BLENDNONE, ///< No blending 
            BLENDSTART, ///< Blend the begininng of the series
            BLENDFINISH, ///< Blend the end of the series
            BLENDBOTH ///< Blend both ends of the series.
        };
        /**
         * Create a single spline with the parameters required to set up
         * the domain and subsequently smooth the given set of y values.
         * The y values must correspond to each of the values in the x array.
         * If either the domain setup fails or the spline cannot be solved,
         * the state will be set to not ok.
         * 
         * The blending span must be smaller than half the total span of
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
        BSplinePlus(const T *x,
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
         * @param blendmode Choose blending at the endpoints.
         * @param blendspan The distance over which the blending is applied.
         */
        BSplinePlus(BSplineBase<T> &base,
                    const T *y,
                    BLENDMODE blendmode = BSplinePlus<T>::BLENDNONE,
                    T blendspan = 0.0);

        /**
         * Return the evaluation of the smoothed curve 
         * at a particular @p x value.  If current state is not ok(), returns 0.
         */
        T evaluate(T x);

        /** 
         * Return the first derivative of the spline curve at the given @p x.
         * Returns zero if the current state is not ok().
         */
        T slope(T x);

        virtual ~BSplinePlus();

    protected:
        /// Initialize the blending data.
        /// @param y The y input vector
        void initBlending(const T* y);
        /// Find the interpolated value of the original Y for the specified X value
        /// starting from the left of the series. If X is not within the series, return 
        /// 0.0
        /// @param x Location at which to do the interpolation
        /// @return The interpolated value.
        T interpLeft(T x);
        /// Find the interpolated value of the original Y for the specified X value
        /// starting from the right of the series.If X is not within the series, return 
        /// 0.0
        /// @param x Location at which to do the interpolation
        /// 2retun The interpolated value.
        T interpRight(T x);
        /// A copy of the input y vector, needed for the blending operation
        std::vector<T> _y;
        /// The blending mode.
        BLENDMODE _blendMode;
        /// The span over which the blending occurs
        T _blendSpan;
        /// The left x value where blending begins
        T _xLeft;
        /// The right x value where blending begins
        T _xRight;
        
        using BSplineBase<T>::OK;
        using BSplineBase<T>::NX;
        using BSplineBase<T>::base;
        using BSplineBase<T>::xmin;
        using BSplineBase<T>::xmax;


};

#endif /*BSPLINEPLUS_H_*/