/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Azra Bano

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file letsberational.hpp
    \brief Normalised Black function and its inverse following Jäckel's
           "Let's Be Rational"

    These are the building blocks of
    blackFormulaImpliedStdDevLetsBeRational(); most users will want
    that function instead.
*/

#ifndef quantlib_lets_be_rational_hpp
#define quantlib_lets_be_rational_hpp

#include <ql/option.hpp>

namespace QuantLib {

    namespace detail {

        /*! Normalised undiscounted Black call price

            \f[ b(x,s) = \Phi(x/s + s/2)\,e^{x/2} - \Phi(x/s - s/2)\,e^{-x/2} \f]

            with log-moneyness \f$ x = \ln(F/K) \f$ and total standard
            deviation \f$ s = \sigma\sqrt{T} \f$; the undiscounted Black
            price is \f$ \sqrt{FK}\,b(x,s) \f$.
        */
        Real normalisedBlackCall(Real x, Real s);

        //! Derivative of normalisedBlackCall() with respect to \f$ s \f$
        Real normalisedVega(Real x, Real s);

        /*! Total standard deviation \f$ s \f$ such that the normalised
            Black price of the given option type equals \f$ \beta \f$,
            using P. Jäckel, "Let's Be Rational", Wilmott (2015), pp. 40-53.

            Prices below the normalised intrinsic value give 0; prices at
            or above the upper bound \f$ e^{\theta x/2} \f$ give QL_MAX_REAL.
            No other checks are made.

            \param iterations if not null, receives the number of
                              Householder iterations that were taken.
        */
        Real normalisedImpliedStdDevLetsBeRational(Real beta,
                                                   Real x,
                                                   Option::Type optionType,
                                                   Size maxIterations = 2,
                                                   Size* iterations = nullptr);

    }

}

#endif
