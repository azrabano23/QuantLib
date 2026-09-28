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

/*
 Most of this file is a port of the reference implementation that
 accompanies P. Jäckel, "Let's Be Rational", Wilmott (2015), pp. 40-53
 (www.jaeckel.org/LetsBeRational.7z, revision 990), which carries the
 following notice:

 ======================================================================================
 Copyright © 2013-2014 Peter Jäckel.

 Permission to use, copy, modify, and distribute this software is freely granted,
 provided that this notice is preserved.

 WARRANTY DISCLAIMER
 The Software is provided "as is" without warranty of any kind, either express or implied,
 including without limitation any implied warranties of condition, uninterrupted use,
 merchantability, fitness for a particular purpose, or non-infringement.
 ======================================================================================

 The complementary error functions are a port of W. J. Cody's CALERF
 (http://www.netlib.org/specfun/erf), as used by the reference
 implementation.
*/

#include <ql/errors.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/pricingengines/letsberational.hpp>
#include <algorithm>
#include <cmath>

namespace QuantLib {

    namespace detail {

        namespace {

            const Real oneOverSqrtTwo = 0.7071067811865475244008443621048490392848359376887;
            const Real oneOverSqrtTwoPi = 0.3989422804014326779399460599343818684758586311649;
            const Real sqrtTwoPi = 2.506628274631000502415765284811045253006986740610;
            const Real twoPi = 6.283185307179586476925286766559005768394338798750;
            const Real sqrtPiOverTwo = 1.253314137315500251207882642405522626503493370305;
            const Real sqrtThree = 1.732050807568877293527446341505872366942805253810;
            const Real sqrtOneOverThree = 0.577350269189625764509148780501957455647601751270;
            const Real twoPiOverSqrtTwentySeven =
                1.209199576156145233729385505094770488189377498728;
            const Real piOverSix = 0.523598775598298873077107230546583814032861566563;

            const Real sqrtEpsilon = std::sqrt(QL_EPSILON);
            const Real fourthRootEpsilon = std::sqrt(sqrtEpsilon);
            const Real eighthRootEpsilon = std::sqrt(fourthRootEpsilon);
            const Real sixteenthRootEpsilon = std::sqrt(eighthRootEpsilon);
            const Real sqrtMinPositive = std::sqrt(QL_MIN_POSITIVE_REAL);
            const Real sqrtMax = std::sqrt(QL_MAX_REAL);

            // η: below h = x/s = η the asymptotic expansion is used
            const Real asymptoticExpansionThreshold = -10.0;
            // τ: below t = s/2 = τ the small-t expansion is used
            const Real smallTExpansionThreshold = 2.0 * sixteenthRootEpsilon;

            Real square(Real x) {
                return x * x;
            }

            // truncation towards zero, i.e., Fortran's AINT
            Real aint(Real x) {
                return x > 0.0 ? Real(std::floor(x)) : Real(-std::floor(-x));
            }

            /* W. J. Cody, "Rational Chebyshev approximations for the error
               function", Math. Comp. 23 (1969), pp. 631-638.  Returns
               erfc(x) or, if scaled is true, erfcx(x) = exp(x^2) erfc(x). */
            Real calerf(Real x, bool scaled) {
                static const Real a[5] = {3.1611237438705656, 113.864154151050156,
                                          377.485237685302021, 3209.37758913846947,
                                          .185777706184603153};
                static const Real b[4] = {23.6012909523441209, 244.024637934444173,
                                          1282.61652607737228, 2844.23683343917062};
                static const Real c[9] = {
                    .564188496988670089, 8.88314979438837594, 66.1191906371416295,
                    298.635138197400131, 881.95222124176909,  1712.04761263407058,
                    2051.07837782607147, 1230.33935479799725, 2.15311535474403846e-8};
                static const Real d[8] = {15.7449261107098347, 117.693950891312499,
                                          537.181101862009858, 1621.38957456669019,
                                          3290.79923573345963, 4362.61909014324716,
                                          3439.36767414372164, 1230.33935480374942};
                static const Real p[6] = {.305326634961232344,    .360344899949804439,
                                          .125781726111229246,    .0160837851487422766,
                                          6.58749161529837803e-4, .0163153871373020978};
                static const Real q[5] = {2.56852019228982242, 1.87295284992346047,
                                          .527905102951428412, .0605183413124413191,
                                          .00233520497626869185};
                const Real sqrpi = 0.56418958354775628695; // 1/sqrt(pi)
                const Real thresh = 0.46875;
                // machine-dependent constants for IEEE double precision
                const Real xinf = 1.79e308, xneg = -26.628, xsmall = 1.11e-16, xbig = 26.543,
                           xhuge = 6.71e7, xmax = 2.53e307;

                const Real y = std::fabs(x);
                Real result;
                if (y <= thresh) {
                    // erfc(x) = 1 - erf(x) for |x| <= 0.46875
                    const Real ysq = y > xsmall ? Real(y * y) : Real(0.0);
                    Real xnum = a[4] * ysq, xden = ysq;
                    for (Size i = 0; i < 3; ++i) {
                        xnum = (xnum + a[i]) * ysq;
                        xden = (xden + b[i]) * ysq;
                    }
                    result = 1.0 - x * (xnum + a[3]) / (xden + b[3]);
                    return scaled ? Real(std::exp(ysq) * result) : result;
                } else if (y <= 4.0) {
                    // erfc for 0.46875 <= |x| <= 4
                    Real xnum = c[8] * y, xden = y;
                    for (Size i = 0; i < 7; ++i) {
                        xnum = (xnum + c[i]) * y;
                        xden = (xden + d[i]) * y;
                    }
                    result = (xnum + c[7]) / (xden + d[7]);
                    if (!scaled) {
                        // exp(-y^2) evaluated as a product to limit round-off
                        const Real ysq = aint(y * 16.0) / 16.0, del = (y - ysq) * (y + ysq);
                        result *= std::exp(-ysq * ysq) * std::exp(-del);
                    }
                } else {
                    // erfc for |x| > 4
                    result = 0.0;
                    if (y >= xbig && (!scaled || y >= xmax)) {
                        // erfc underflows, or erfcx is below the smallest number
                    } else if (y >= xbig && y >= xhuge) {
                        result = sqrpi / y;
                    } else {
                        const Real ysq = 1.0 / (y * y);
                        Real xnum = p[5] * ysq, xden = ysq;
                        for (Size i = 0; i < 4; ++i) {
                            xnum = (xnum + p[i]) * ysq;
                            xden = (xden + q[i]) * ysq;
                        }
                        result = ysq * (xnum + p[4]) / (xden + q[4]);
                        result = (sqrpi - result) / y;
                        if (!scaled) {
                            const Real ysq2 = aint(y * 16.0) / 16.0, del = (y - ysq2) * (y + ysq2);
                            result *= std::exp(-ysq2 * ysq2) * std::exp(-del);
                        }
                    }
                }
                // fix up for negative argument
                if (x < 0.0) {
                    if (!scaled) {
                        result = 2.0 - result;
                    } else if (x < xneg) {
                        result = xinf;
                    } else {
                        const Real ysq = aint(x * 16.0) / 16.0, del = (x - ysq) * (x + ysq);
                        const Real e = std::exp(ysq * ysq) * std::exp(del);
                        result = (e + e) - result;
                    }
                }
                return result;
            }

            Real erfc(Real x) {
                return calerf(x, false);
            }
            Real erfcx(Real x) {
                return calerf(x, true);
            }

            Real normPdf(Real x) {
                return oneOverSqrtTwoPi * std::exp(-0.5 * x * x);
            }

            Real normCdf(Real z) {
                /* For z <= -10, the asymptotic expansion
                   Phi(z) = phi(z)/|z| [1 - 1/z^2 + ...], Abramowitz & Stegun (26.2.12),
                   has a relative accuracy of 1.64e-16 with at most 17 terms. */
                if (z <= -10.0) {
                    Real sum = 1.0;
                    if (z >= -1.0 / sqrtEpsilon) {
                        const Real zsqr = z * z;
                        Real i = 1.0, g = 1.0, a = QL_MAX_REAL, lasta;
                        do {
                            lasta = a;
                            const Real x = (4.0 * i - 3.0) / zsqr;
                            const Real y = x * ((4.0 * i - 1.0) / zsqr);
                            a = g * (x - y);
                            sum -= a;
                            g *= y;
                            ++i;
                            a = std::fabs(a);
                        } while (lasta > a && a >= std::fabs(sum * QL_EPSILON));
                    }
                    return -normPdf(z) * sum / z;
                }
                return 0.5 * erfc(-z * oneOverSqrtTwo);
            }

            Real inverseNormCdf(Real u) {
                // only used for the initial guess; the Householder
                // iterations take care of its ~1e-9 relative error.
                return InverseCumulativeNormal::standard_value(u);
            }

            /* normalised intrinsic value max(theta (e^{x/2} - e^{-x/2}), 0);
               a Taylor expansion avoids cancellation for small |x|. */
            Real normalisedIntrinsic(Real x, Real theta) {
                if (theta * x <= 0.0)
                    return 0.0;
                const Real x2 = x * x;
                // the factor 98 is the fourth root of the last coefficient, 92897280
                if (x2 < 98.0 * fourthRootEpsilon)
                    return std::fabs(std::max(
                        Real((theta < 0.0 ? -1.0 : 1.0) * x *
                             (1.0 + x2 * ((1.0 / 24.0) +
                                          x2 * ((1.0 / 1920.0) + x2 * ((1.0 / 322560.0) +
                                                                       (1.0 / 92897280.0) * x2))))),
                        Real(0.0)));
                const Real bMax = std::exp(0.5 * x), oneOverBMax = 1.0 / bMax;
                return std::fabs(
                    std::max(Real((theta < 0.0 ? -1.0 : 1.0) * (bMax - oneOverBMax)), Real(0.0)));
            }

            /* Asymptotic expansion of
                   b = Phi(h+t) exp(h t) - Phi(h-t) exp(-h t)
                     = exp(-(h^2+t^2)/2) / sqrt(2 pi) [Y(h+t) - Y(h-t)],  Y(z) = Phi(z)/phi(z)
               for large negative h - t (region 1), with h = x/s and t = s/2,
               based on Abramowitz & Stegun (26.2.12).
               With r = (h+t)(h-t), A(h,t) = r/t [Y(h+t) - Y(h-t)] is
               expanded to 17th order in q = (h/r)^2 as
                   A = 2 + q (P_1 + 3 q (P_2 + 5 q (P_3 + ... + 33 q P_17)))
               where P_k(e) = (-1)^k 2 sum_{j=0..k} C(2k+1, 2j+1) e^j and
               e = (t/h)^2.  This is the same polynomial that the reference
               implementation spells out term by term. */
            Real asymptoticExpansionOfNormalisedBlackCall(Real h, Real t) {
                const Real e = (t / h) * (t / h), r = ((h + t) * (h - t)), q = (h / r) * (h / r);
                const Integer order = 17;
                Real sum = 0.0;
                for (Integer k = order; k >= 1; --k) {
                    // P_k(e) by Horner's scheme, from the highest power of e;
                    // the binomial coefficients are exact integers in Real.
                    const Real n = 2 * k + 1;
                    Real binomial = 1.0; // C(2k+1, 2k+1)
                    Real p = binomial;
                    for (Integer m = 2 * k + 1; m > 1; m -= 2) {
                        // C(n, m-2), exact: the division leaves no remainder
                        binomial = binomial * m * (m - 1) / ((n - m + 2) * (n - m + 1));
                        p = binomial + e * p;
                    }
                    p *= (k % 2 == 0) ? 2.0 : -2.0;
                    sum = (k == order) ? p : Real(p + (2 * k + 1) * q * sum);
                }
                sum = 2.0 + q * sum;
                const Real b = oneOverSqrtTwoPi * std::exp(-0.5 * (h * h + t * t)) * (t / r) * sum;
                return std::max(b, Real(0.0));
            }

            /* Expansion of Y(h+t) - Y(h-t) to twelfth order in t (region 2);
               accurate to machine precision for h <= 0 and t < τ = 2 ε^(1/16) ≈ 0.21.
               With Y(h) = sqrt(pi/2) erfcx(-h/sqrt(2)), the coefficient
               a = 1 + h Y(h) is positive and goes to 0 as h -> -inf; it is
               the main source of round-off when |h| > 1. */
            Real smallTExpansionOfNormalisedBlackCall(Real h, Real t) {
                const Real a = 1.0 + h * (0.5 * sqrtTwoPi) * erfcx(-oneOverSqrtTwo * h), w = t * t,
                           h2 = h * h;
                // coefficients of w^k, k = 1...6
                const Real c1 = (-1 + 3 * a + a * h2) / 6;
                const Real c2 = (-7 + 15 * a + h2 * (-1 + 10 * a + a * h2)) / 120;
                const Real c3 =
                    (-57 + 105 * a + h2 * (-18 + 105 * a + h2 * (-1 + 21 * a + a * h2))) / 5040;
                const Real c4 =
                    (-561 + 945 * a +
                     h2 * (-285 + 1260 * a + h2 * (-33 + 378 * a + h2 * (-1 + 36 * a + a * h2)))) /
                    362880;
                const Real c5 = (-6555 + 10395 * a +
                                 h2 * (-4680 + 17325 * a +
                                       h2 * (-840 + 6930 * a +
                                             h2 * (-52 + 990 * a + h2 * (-1 + 55 * a + a * h2))))) /
                                39916800;
                const Real n6 =
                    -89055 + 135135 * a +
                    h2 * (-82845 + 270270 * a +
                          h2 * (-20370 + 135135 * a +
                                h2 * (-1926 + 25740 * a +
                                      h2 * (-75 + 2145 * a + h2 * (-1 + 78 * a + a * h2)))));
                const Real expansion =
                    2 * t *
                    (a +
                     w * (c1 + w * (c2 + w * (c3 + w * (c4 + w * (c5 + n6 * w / 6227020800.0))))));
                const Real b = oneOverSqrtTwoPi * std::exp(-0.5 * (h * h + t * t)) * expansion;
                return std::max(b, Real(0.0));
            }

            /* Region 3: b is dominated by its first term, so the direct
               formula loses nothing. */
            Real normalisedBlackCallUsingNormCdf(Real x, Real s) {
                const Real h = x / s, t = 0.5 * s, bMax = std::exp(0.5 * x);
                const Real b = normCdf(h + t) * bMax - normCdf(h - t) / bMax;
                return std::max(b, Real(0.0));
            }

            /* Region 4: with Phi(z) = exp(-z^2/2) erfcx(-z/sqrt(2)) / 2,
                   b = exp(-(h^2+t^2)/2) [erfcx(-(h+t)/sqrt(2)) - erfcx(-(h-t)/sqrt(2))] / 2,
               which needs a single exponential instead of four. */
            Real normalisedBlackCallUsingErfcx(Real h, Real t) {
                const Real b =
                    0.5 * std::exp(-0.5 * (h * h + t * t)) *
                    (erfcx(-oneOverSqrtTwo * (h + t)) - erfcx(-oneOverSqrtTwo * (h - t)));
                return std::max(b, Real(0.0));
            }

            /* Rational cubic interpolation, following R. Delbourgo and
               J. A. Gregory, "Shape preserving piecewise rational
               interpolation", SIAM J. Sci. Stat. Comput. 6 (1985). */

            const Real minimumRationalCubicControlParameter = -(1.0 - sqrtEpsilon);
            const Real maximumRationalCubicControlParameter = 2.0 / (QL_EPSILON * QL_EPSILON);

            bool isZero(Real x) {
                return std::fabs(x) < QL_MIN_POSITIVE_REAL;
            }

            Real rationalCubicInterpolation(
                Real x, Real xL, Real xR, Real yL, Real yR, Real dL, Real dR, Real r) {
                const Real h = xR - xL;
                if (std::fabs(h) <= 0.0)
                    return 0.5 * (yL + yR);
                const Real t = (x - xL) / h;
                // r >= maximum means linear interpolation (this also lets NaN propagate)
                if (!(r >= maximumRationalCubicControlParameter)) {
                    const Real omt = 1.0 - t, t2 = t * t, omt2 = omt * omt;
                    // formula (2.4) divided by formula (2.5)
                    return (yR * t2 * t + (r * yR - h * dR) * t2 * omt +
                            (r * yL + h * dL) * t * omt2 + yL * omt2 * omt) /
                           (1.0 + (r - 3.0) * t * omt);
                }
                return yR * t + yL * (1.0 - t);
            }

            Real rationalCubicControlParameterToFitSecondDerivativeAtLeftSide(
                Real xL, Real xR, Real yL, Real yR, Real dL, Real dR, Real secondDerivativeL) {
                const Real h = xR - xL, numerator = 0.5 * h * secondDerivativeL + (dR - dL);
                if (isZero(numerator))
                    return 0.0;
                const Real denominator = (yR - yL) / h - dL;
                if (isZero(denominator))
                    return numerator > 0.0 ? maximumRationalCubicControlParameter :
                                             minimumRationalCubicControlParameter;
                return numerator / denominator;
            }

            Real rationalCubicControlParameterToFitSecondDerivativeAtRightSide(
                Real xL, Real xR, Real yL, Real yR, Real dL, Real dR, Real secondDerivativeR) {
                const Real h = xR - xL, numerator = 0.5 * h * secondDerivativeR + (dR - dL);
                if (isZero(numerator))
                    return 0.0;
                const Real denominator = dR - (yR - yL) / h;
                if (isZero(denominator))
                    return numerator > 0.0 ? maximumRationalCubicControlParameter :
                                             minimumRationalCubicControlParameter;
                return numerator / denominator;
            }

            // smallest control parameter that preserves monotonicity and convexity
            Real minimumRationalCubicControlParameterValue(Real dL,
                                                           Real dR,
                                                           Real s,
                                                           bool preferShapePreservation) {
                const bool monotonic = dL * s >= 0.0 && dR * s >= 0.0, convex = dL <= s && s <= dR,
                           concave = dL >= s && s >= dR;
                if (!monotonic && !convex && !concave)
                    return minimumRationalCubicControlParameter;
                const Real dRmdL = dR - dL, dRms = dR - s, smdL = s - dL;
                Real r1 = -QL_MAX_REAL, r2 = r1;
                if (monotonic) {
                    if (!isZero(s))
                        r1 = (dR + dL) / s; // (3.8)
                    else if (preferShapePreservation)
                        r1 = maximumRationalCubicControlParameter;
                }
                if (convex || concave) {
                    if (!(isZero(smdL) || isZero(dRms))) // (3.18)
                        r2 = std::max(std::fabs(dRmdL / dRms), std::fabs(dRmdL / smdL));
                    else if (preferShapePreservation)
                        r2 = maximumRationalCubicControlParameter;
                } else if (monotonic && preferShapePreservation) {
                    r2 = maximumRationalCubicControlParameter;
                }
                return std::max(minimumRationalCubicControlParameter, std::max(r1, r2));
            }

            Real convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide(
                Real xL,
                Real xR,
                Real yL,
                Real yR,
                Real dL,
                Real dR,
                Real secondDerivativeL,
                bool preferShapePreservation) {
                const Real r = rationalCubicControlParameterToFitSecondDerivativeAtLeftSide(
                    xL, xR, yL, yR, dL, dR, secondDerivativeL);
                const Real rMin = minimumRationalCubicControlParameterValue(
                    dL, dR, (yR - yL) / (xR - xL), preferShapePreservation);
                return std::max(r, rMin);
            }

            Real convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide(
                Real xL,
                Real xR,
                Real yL,
                Real yR,
                Real dL,
                Real dR,
                Real secondDerivativeR,
                bool preferShapePreservation) {
                const Real r = rationalCubicControlParameterToFitSecondDerivativeAtRightSide(
                    xL, xR, yL, yR, dL, dR, secondDerivativeR);
                const Real rMin = minimumRationalCubicControlParameterValue(
                    dL, dR, (yR - yL) / (xR - xL), preferShapePreservation);
                return std::max(r, rMin);
            }

            /* Lower map f(beta) = 2 pi |x| / sqrt(27) Phi(-|x|/(sqrt(3) s))^3
               and its first two derivatives with respect to beta, at s. */
            void lowerMapAndDerivatives(Real x, Real s, Real& f, Real& fp, Real& fpp) {
                const Real ax = std::fabs(x), z = sqrtOneOverThree * ax / s, y = z * z, s2 = s * s,
                           Phi = normCdf(-z), phi = normPdf(z);
                fpp =
                    piOverSix * y / (s2 * s) * Phi *
                    (8.0 * sqrtThree * s * ax + (3.0 * s2 * (s2 - 8.0) - 8.0 * x * x) * Phi / phi) *
                    std::exp(2.0 * y + 0.25 * s2);
                const Real Phi2 = Phi * Phi;
                fp = twoPi * y * Phi2 * std::exp(y + 0.125 * s * s);
                f = twoPiOverSqrtTwentySeven * ax * (Phi2 * Phi);
            }

            Real inverseLowerMap(Real x, Real f) {
                return std::fabs(x / (sqrtThree * inverseNormCdf(std::pow(
                                                      f / (twoPiOverSqrtTwentySeven * std::fabs(x)),
                                                      1.0 / 3.0))));
            }

            /* Upper map f(beta) = Phi(-s/2) and its first two derivatives
               with respect to beta, at s. */
            void upperMapAndDerivatives(Real x, Real s, Real& f, Real& fp, Real& fpp) {
                f = normCdf(-0.5 * s);
                const Real w = square(x / s);
                fp = -0.5 * std::exp(0.5 * w);
                fpp = sqrtPiOverTwo * std::exp(w + 0.125 * s * s) * w / s;
            }

            Real inverseUpperMap(Real f) {
                return -2.0 * inverseNormCdf(f);
            }

            /* Third-order Householder step: with newton = -g/g',
               halley = g''/g' and hh3 = g'''/g', the step is
               newton * (1 + halley newton / 2) / (1 + newton (halley + hh3 newton / 6)). */
            Real householderFactor(Real newton, Real halley, Real hh3) {
                return (1.0 + 0.5 * halley * newton) /
                       (1.0 + newton * (halley + hh3 * newton / 6.0));
            }

            /* Safeguarded iteration shared by the three objective functions.
               The bracket [sLeft, sRight] is tightened at every step; if the
               iteration leaves it, or changes direction three times, it falls
               back to bisection.  This only happens for extreme inputs such as
               |x| > 500. */
            template <class Step>
            Real householderIteration(Real x,
                                      Real beta,
                                      Real s,
                                      Real sLeft,
                                      Real sRight,
                                      Size maxIterations,
                                      Size& iterations,
                                      const Step& step) {
                Real ds = QL_MAX_REAL, dsPrevious = 0.0;
                Size directionReversals = 0;
                for (iterations = 0; iterations < maxIterations && std::fabs(ds) > QL_EPSILON * s;
                     ++iterations) {
                    if (ds * dsPrevious < 0.0)
                        ++directionReversals;
                    if (iterations > 0 && (directionReversals == 3 || !(s > sLeft && s < sRight))) {
                        s = 0.5 * (sLeft + sRight);
                        if (sRight - sLeft <= QL_EPSILON * s)
                            break;
                        directionReversals = 0;
                        ds = 0.0;
                    }
                    dsPrevious = ds;
                    const Real b = normalisedBlackCall(x, s), bp = normalisedVega(x, s);
                    if (b > beta && s < sRight)
                        sRight = s;
                    else if (b < beta && s > sLeft)
                        sLeft = s;
                    ds = std::max(Real(-0.5 * s), step(s, b, bp, sLeft, sRight));
                    s += ds;
                }
                return s;
            }

        }

        Real normalisedBlackCall(Real x, Real s) {
            // in-the-money: intrinsic value plus out-of-the-money value
            if (x > 0.0)
                return normalisedIntrinsic(x, 1.0) + normalisedBlackCall(-x, s);
            if (s <= 0.0)
                return 0.0;
            // With h = x/s and t = s/2, the conditions below are written so as
            // to avoid divisions by s.
            // Region 1: h < η and t < τ + |h| - |η|
            if (x < s * asymptoticExpansionThreshold &&
                0.5 * s * s + x < s * (smallTExpansionThreshold + asymptoticExpansionThreshold))
                return asymptoticExpansionOfNormalisedBlackCall(x / s, 0.5 * s);
            // Region 2: t < τ
            if (0.5 * s < smallTExpansionThreshold)
                return smallTExpansionOfNormalisedBlackCall(x / s, 0.5 * s);
            // Region 3: h + t > 0.85, i.e., b is more than about 85% of exp(x/2)
            if (x + 0.5 * s * s > s * 0.85)
                return normalisedBlackCallUsingNormCdf(x, s);
            // Region 4
            return normalisedBlackCallUsingErfcx(x / s, 0.5 * s);
        }

        Real normalisedVega(Real x, Real s) {
            const Real ax = std::fabs(x);
            if (ax <= 0.0)
                return oneOverSqrtTwoPi * std::exp(-0.125 * s * s);
            if (s <= 0.0 || s <= ax * sqrtMinPositive)
                return 0.0;
            return oneOverSqrtTwoPi * std::exp(-0.5 * (square(x / s) + square(0.5 * s)));
        }

        Real normalisedImpliedStdDevLetsBeRational(
            Real beta, Real x, Option::Type optionType, Size maxIterations, Size* iterations) {
            Size n = 0;
            if (iterations != nullptr)
                *iterations = 0;
            Real theta = optionType == Option::Call ? 1.0 : -1.0;
            // subtract the intrinsic value: in-the-money to out-of-the-money
            if (theta * x > 0.0) {
                beta = std::fabs(std::max(Real(beta - normalisedIntrinsic(x, theta)), Real(0.0)));
                theta = -theta;
            }
            // map puts to calls: b_put(x) = b_call(-x)
            if (theta < 0.0)
                x = -x;
            // from here on, x <= 0 and beta is the price of an out-of-the-money call
            if (beta <= 0.0)
                return 0.0;
            const Real bMax = std::exp(0.5 * x);
            if (beta >= bMax)
                return QL_MAX_REAL;

            // central point: the inflexion point of b(s) is at s_c = sqrt(2|x|)
            const Real sC = std::sqrt(std::fabs(2.0 * x)), bC = normalisedBlackCall(x, sC),
                       vC = normalisedVega(x, sC);

            Real s, sLeft = QL_MIN_POSITIVE_REAL, sRight = QL_MAX_REAL;

            if (beta < bC) {
                // s_l: where the tangent at s_c crosses b = 0
                const Real sL = sC - bC / vC, bL = normalisedBlackCall(x, sL);
                if (beta < bL) {
                    // lowest segment: rational cubic in the lower map f(beta)
                    Real fL, dfL, d2fL;
                    lowerMapAndDerivatives(x, sL, fL, dfL, d2fL);
                    const Real rLL =
                        convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide(
                            0.0, bL, 0.0, fL, 1.0, dfL, d2fL, true);
                    Real f = rationalCubicInterpolation(beta, 0.0, bL, 0.0, fL, 1.0, dfL, rLL);
                    if (!(f > 0.0)) {
                        // round-off for extreme |x|: quadratic interpolation
                        // through f(0) = 0, f'(0) = 1 and f(b_l)
                        const Real t = beta / bL;
                        f = (fL * t + bL * (1.0 - t)) * t;
                    }
                    s = inverseLowerMap(x, f);
                    sRight = sL;
                    /* objective g(s) = 1/ln(b(s)) - 1/ln(beta), which is
                       close to linear in s where b is exponentially small:
                         newton = -g/g'  = (ln(beta)-ln(b)) ln(b)/ln(beta) b/b'
                         halley = g''/g' = b''/b' - b'/b (1+2/ln(b))
                         hh3 = g'''/g'   = b'''/b' + 2(b'/b)^2 (1+3/ln(b)(1+1/ln(b)))
                                           - 3 (b''/b') (b'/b) (1+2/ln(b)) */
                    const Real lnBeta = std::log(beta);
                    s = householderIteration(
                        x, beta, s, sLeft, sRight, maxIterations, n,
                        [x, lnBeta](Real s, Real b, Real bp, Real sLeft, Real sRight) -> Real {
                            if (b <= 0.0 || bp <= 0.0) // underflow: bisection
                                return 0.5 * (sLeft + sRight) - s;
                            const Real lnB = std::log(b), bpob = bp / b, h = x / s,
                                       bHalley = h * h / s - s / 4.0,
                                       newton = (lnBeta - lnB) * lnB / lnBeta / bpob,
                                       halley = bHalley - bpob * (1.0 + 2.0 / lnB);
                            const Real bHH3 = bHalley * bHalley - 3.0 * square(h / s) - 0.25,
                                       hh3 = bHH3 +
                                             2.0 * square(bpob) *
                                                 (1.0 + 3.0 / lnB * (1.0 + 1.0 / lnB)) -
                                             3.0 * bHalley * bpob * (1.0 + 2.0 / lnB);
                            return newton * householderFactor(newton, halley, hh3);
                        });
                    if (iterations != nullptr)
                        *iterations = n;
                    return s;
                } else {
                    // lower middle segment: rational cubic in s(beta)
                    const Real vL = normalisedVega(x, sL);
                    const Real rLM =
                        convexRationalCubicControlParameterToFitSecondDerivativeAtRightSide(
                            bL, bC, sL, sC, 1.0 / vL, 1.0 / vC, 0.0, false);
                    s = rationalCubicInterpolation(beta, bL, bC, sL, sC, 1.0 / vL, 1.0 / vC, rLM);
                    sLeft = sL;
                    sRight = sC;
                }
            } else {
                // s_h: where the tangent at s_c crosses b = b_max
                const Real sH = vC > QL_MIN_POSITIVE_REAL ? Real(sC + (bMax - bC) / vC) : sC,
                           bH = normalisedBlackCall(x, sH);
                if (beta <= bH) {
                    // upper middle segment: rational cubic in s(beta)
                    const Real vH = normalisedVega(x, sH);
                    const Real rHM =
                        convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide(
                            bC, bH, sC, sH, 1.0 / vC, 1.0 / vH, 0.0, false);
                    s = rationalCubicInterpolation(beta, bC, bH, sC, sH, 1.0 / vC, 1.0 / vH, rHM);
                    sLeft = sC;
                    sRight = sH;
                } else {
                    // highest segment: rational cubic in the upper map f(beta)
                    Real fH, dfH, d2fH, f = -QL_MAX_REAL;
                    upperMapAndDerivatives(x, sH, fH, dfH, d2fH);
                    if (d2fH > -sqrtMax && d2fH < sqrtMax) {
                        const Real rHH =
                            convexRationalCubicControlParameterToFitSecondDerivativeAtLeftSide(
                                bH, bMax, fH, 0.0, dfH, -0.5, d2fH, true);
                        f = rationalCubicInterpolation(beta, bH, bMax, fH, 0.0, dfH, -0.5, rHH);
                    }
                    if (f <= 0.0) {
                        // quadratic interpolation through f(b_h), f(b_max) = 0
                        // and f'(b_max) = -1/2
                        const Real h = bMax - bH, t = (beta - bH) / h;
                        f = (fH * (1.0 - t) + 0.5 * h * t) * (1.0 - t);
                    }
                    s = inverseUpperMap(f);
                    sLeft = sH;
                    if (beta > 0.5 * bMax) {
                        /* objective g(s) = ln(b_max-beta) - ln(b_max-b(s)):
                             newton = -g/g'  = ln((b_max-b)/(b_max-beta)) (b_max-b)/b'
                             halley = g''/g' = b''/b' + g'
                             hh3 = g'''/g'   = b'''/b' + g' (2g' + 3b''/b') */
                        s = householderIteration(
                            x, beta, s, sLeft, sRight, maxIterations, n,
                            [x, beta, bMax](Real s, Real b, Real bp, Real sLeft,
                                            Real sRight) -> Real {
                                if (b >= bMax || bp <= QL_MIN_POSITIVE_REAL) // bisection
                                    return 0.5 * (sLeft + sRight) - s;
                                const Real bMaxMinusB = bMax - b,
                                           g = std::log((bMax - beta) / bMaxMinusB),
                                           gp = bp / bMaxMinusB;
                                const Real bHalley = square(x / s) / s - s / 4.0,
                                           bHH3 =
                                               bHalley * bHalley - 3.0 * square(x / (s * s)) - 0.25;
                                const Real newton = -g / gp, halley = bHalley + gp,
                                           hh3 = bHH3 + gp * (2.0 * gp + 3.0 * bHalley);
                                return newton * householderFactor(newton, halley, hh3);
                            });
                        if (iterations != nullptr)
                            *iterations = n;
                        return s;
                    }
                }
            }
            /* middle segments: objective g(s) = b(s) - beta
                 newton = -g/g'  = (beta-b)/b'
                 halley = g''/g' = b''/b' = x^2/s^3 - s/4
                 hh3 = g'''/g'   = b'''/b' = halley^2 - 3 (x/s^2)^2 - 1/4 */
            s = householderIteration(x, beta, s, sLeft, sRight, maxIterations, n,
                                     [x, beta](Real s, Real b, Real bp, Real, Real) -> Real {
                                         const Real newton = (beta - b) / bp,
                                                    halley = square(x / s) / s - s / 4.0,
                                                    hh3 = halley * halley -
                                                          3.0 * square(x / (s * s)) - 0.25;
                                         return newton * householderFactor(newton, halley, hh3);
                                     });
            if (iterations != nullptr)
                *iterations = n;
            return s;
        }

    }

}
