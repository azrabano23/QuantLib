/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2013 Gary Kennedy
 Copyright (C) 2015, 2024 Peter Caspers
 Copyright (C) 2017 Klaus Spanderen
 Copyright (C) 2020 Marcin Rybacki
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


#include "toplevelfixture.hpp"
#include "utilities.hpp"
#include <ql/pricingengines/blackformula.hpp>
#include <ql/pricingengines/letsberational.hpp>
#include <cmath>
#include <iomanip>

using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)

BOOST_AUTO_TEST_SUITE(BlackFormulaTests)

BOOST_AUTO_TEST_CASE(testBachelierImpliedVol) {

    BOOST_TEST_MESSAGE("Testing Bachelier implied vol...");

    Real forward = 1.0;
    Real bpvol = 0.01;
    Real tte = 10.0;
    Real stdDev = bpvol * std::sqrt(tte);
    Option::Type optionType = Option::Call;
    Real discount = 0.95;

    Real d[] = {-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0};
    for (Real i : d) {


        Real strike = forward - i * bpvol * std::sqrt(tte);

        Real callPrem = bachelierBlackFormula(optionType, strike, forward, stdDev, discount);

        Real impliedBpVol = bachelierBlackFormulaImpliedVolChoi(optionType, strike, forward, tte,
                                                                callPrem, discount);

        if (std::fabs(bpvol - impliedBpVol) > 1.0e-12) {
            BOOST_ERROR("Failed, expected " << bpvol << " realised " << impliedBpVol);
        }

        Real impliedBpVolExact =
            bachelierBlackFormulaImpliedVol(optionType, strike, forward, tte, callPrem, discount);

        if (std::fabs(bpvol - impliedBpVolExact) > 1.0e-15) {
            BOOST_ERROR("Failed, expected " << bpvol << " realised " << impliedBpVolExact);
        }
    }
}

BOOST_AUTO_TEST_CASE(testChambersImpliedVol) {

    BOOST_TEST_MESSAGE("Testing Chambers-Nawalkha implied vol approximation...");

    Option::Type types[] = {Option::Call, Option::Put};
    Real displacements[] = {0.0000, 0.0010, 0.0050, 0.0100, 0.0200};
    Real forwards[] = {-0.0010, 0.0000, 0.0050, 0.0100, 0.0200, 0.0500};
    Real strikes[] = {-0.0100, -0.0050, -0.0010, 0.0000, 0.0010, 0.0050,
                      0.0100,  0.0200,  0.0500,  0.1000};
    Real stdDevs[] = {0.10, 0.15, 0.20, 0.30, 0.50, 0.60, 0.70,
                      0.80, 1.00, 1.50, 2.00};
    Real discounts[] = {1.00, 0.95, 0.80, 1.10};

    Real tol = 5.0E-4;

    for (auto& type : types) {
        for (Real& displacement : displacements) {
            for (Real& forward : forwards) {
                for (Real& strike : strikes) {
                    for (Real& stdDev : stdDevs) {
                        for (Real& discount : discounts) {
                            if (forward + displacement > 0.0 && strike + displacement > 0.0) {
                                Real premium = blackFormula(type, strike, forward, stdDev, discount,
                                                            displacement);
                                Real atmPremium = blackFormula(type, forward, forward, stdDev,
                                                               discount, displacement);
                                Real iStdDev = blackFormulaImpliedStdDevChambers(
                                    type, strike, forward, premium, atmPremium, discount,
                                    displacement);
                                Real moneyness = (strike + displacement) / (forward + displacement);
                                if(moneyness > 1.0) moneyness = 1.0 / moneyness;
                                Real error = (iStdDev - stdDev) / stdDev * moneyness;
                                if(error > tol)
                                    BOOST_ERROR("Failed to verify Chambers-Nawalkha "
                                                "approximation for "
                                                << type << " displacement=" << displacement
                                                << " forward=" << forward << " strike=" << strike
                                                << " discount=" << discount << " stddev=" << stdDev
                                                << " result=" << iStdDev
                                                << " exceeds maximum error tolerance");
                            }
                        }
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testRadoicicStefanicaImpliedVol) {

    BOOST_TEST_MESSAGE(
        "Testing Radoicic-Stefanica implied vol approximation...");

    const Time T = 1.7;
    const Rate r = 0.1;
    const DiscountFactor df = std::exp(-r*T);

    const Real forward = 100;

    const Volatility vol = 0.3;
    const Real stdDev = vol * std::sqrt(T);

    const Option::Type types[] = { Option::Call, Option::Put };
    const Real strikes[] = {
        50, 60, 70, 80, 90, 100, 110, 125, 150, 200, 300 };

    const Real tol = 0.02;

    for (Real strike : strikes) {
        for (auto type : types) {
            const ext::shared_ptr<PlainVanillaPayoff> payoff(
                ext::make_shared<PlainVanillaPayoff>(type, strike));

            const Real marketValue = blackFormula(payoff, forward, stdDev, df);

            const Real estVol = blackFormulaImpliedStdDevApproximationRS(
                payoff, forward, marketValue, df) / std::sqrt(T);

            const Real error = std::fabs(estVol - vol);
            if (error > tol) {
                BOOST_ERROR("Failed to verify Radoicic-Stefanica"
                    "approximation for "
                    << type
                    << "\n forward     :" << forward
                    << "\n strike      :" << strike
                    << "\n discount    :" << df
                    << "\n implied vol :" << vol
                    << "\n result      :" << estVol
                    << "\n error       :" << error
                    << "\n tolerance   :" << tol);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testRadoicicStefanicaLowerBound) {

    BOOST_TEST_MESSAGE("Testing Radoicic-Stefanica lower bound...");

    // testing lower bound plot figure 3.1 from
    // "Tighter Bounds for Implied Volatility",
    // J. Gatheral, I. Matic, R. Radoicic, D. Stefanica
    // https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2922742

    const Real forward = 1.0;
    const Real k = 1.2;

    for (Real s=0.17; s < 2.9; s+=0.01) {
        const Real strike = std::exp(k)*forward;
        const Real c = blackFormula(Option::Call, strike, forward, s);
        const Real estimate = blackFormulaImpliedStdDevApproximationRS(
            Option::Call, strike, forward, c);

        const Real error = s - estimate;
        if (std::isnan(estimate) || std::fabs(error) > 0.05) {
            BOOST_ERROR("Failed to lower bound Radoicic-Stefanica"
                "approximation for "
                << "\n forward     :" << forward
                << "\n strike      :" << k
                << "\n stdDev      :" << s
                << "\n result      :" << estimate
                << "\n error       :" << error);

        }

        if (c > 1e-6 && error < 0.0) {
            BOOST_ERROR("Failed to verify Radoicic-Stefanica is lower bound"
                    << "\n forward     :" << forward
                    << "\n strike      :" << k
                    << "\n stdDev      :" << s
                    << "\n result      :" << estimate
                    << "\n error       :" << error);
        }
    }
}

BOOST_AUTO_TEST_CASE(testImpliedVolAdaptiveSuccessiveOverRelaxation) {
    BOOST_TEST_MESSAGE("Testing implied volatility calculation via "
        "adaptive successive over-relaxation...");

    const DayCounter dc = Actual365Fixed();
    const Date today = Date(12, July, 2017);
    Settings::instance().evaluationDate() = today;

    const Date exerciseDate = today + Period(15, Months);
    const Time exerciseTime = dc.yearFraction(today, exerciseDate);

    const ext::shared_ptr<YieldTermStructure> rTS = flatRate(0.10, dc);
    const ext::shared_ptr<YieldTermStructure> qTS = flatRate(0.06, dc);

    const DiscountFactor df = rTS->discount(exerciseDate);

    const Volatility vol = 0.20;
    const Real stdDev = vol * std::sqrt(exerciseTime);

    const Real s0     = 100;
    const Real forward= s0 * qTS->discount(exerciseDate)/df;

    const Option::Type types[] = { Option::Call, Option::Put };
    const Real strikes[] = { 50, 60, 70, 80, 90, 100, 110, 125, 150, 200 };
    const Real displacements[] = { 0, 25, 50, 100};

    const Real tol = 1e-8;

    for (Real strike : strikes) {
        for (auto type : types) {
            const ext::shared_ptr<PlainVanillaPayoff> payoff(
                ext::make_shared<PlainVanillaPayoff>(type, strike));

            for (Real displacement : displacements) {

                const Real marketValue = blackFormula(payoff, forward, stdDev, df, displacement);

                const Real impliedStdDev = blackFormulaImpliedStdDevLiRS(
                    payoff, forward, marketValue, df, displacement,
                    Null<Real>(), 1.0, tol, 100);

                const Real error = std::fabs(impliedStdDev - stdDev);
                if (error > 10*tol) {
                    BOOST_ERROR("Failed to calculated implied volatility"
                                " with adaptive successive over-relaxation"
                            << "\n forward     :" << forward
                            << "\n strike      :" << strike
                            << "\n stdDev      :" << stdDev
                            << "\n displacement:" << displacement
                            << "\n result      :" << impliedStdDev
                            << "\n error       :" << error
                            << "\n tolerance   :" << tol);
                }
            }
        }
    }
}

void assertBlackFormulaForwardDerivative(
    Option::Type optionType,
    const std::vector<Real> &strikes,
    Real bpvol)
{
    Real forward = 1.0;
    Real tte = 10.0;
    Real stdDev = bpvol * std::sqrt(tte);
    Real discount = 0.95;
    Real displacement = 0.01;
    Real bump = 0.0001;
    Real epsilon = 1.e-10;
    std::string type = optionType == Option::Call ? "Call" : "Put";

    for (Real strike : strikes) {
        Real delta = blackFormulaForwardDerivative(optionType, strike, forward, stdDev, discount,
                                                   displacement);
        Real bumpedDelta = blackFormulaForwardDerivative(
            optionType, strike, forward + bump, stdDev, discount, displacement);

        Real basePremium = blackFormula(
            optionType, strike, forward, stdDev, discount, displacement);
        Real bumpedPremium = blackFormula(
            optionType, strike, forward + bump, stdDev, discount, displacement);
        Real deltaApprox = (bumpedPremium - basePremium) / bump;

        /*! Based on the Mean Value Theorem, the below inequality
            should hold for any function that is monotonic in the
            area of the bump.
         */
        bool success =
            (std::max(delta, bumpedDelta) + epsilon > deltaApprox)
            &&  (deltaApprox > std::min(delta, bumpedDelta) - epsilon);

        if (!success)
        {
            BOOST_ERROR("Failed to calculate the derivative of the"
                        " Black formula w.r.t. forward"
                        << "\n option type       :" << type
                        << "\n forward           :" << forward
                        << "\n strike            :" << strike
                        << "\n stdDev            :" << stdDev
                        << "\n displacement      :" << displacement
                        << "\n analytical delta  :" << delta
                        << "\n approximated delta:" << deltaApprox);
        }
    }
}

BOOST_AUTO_TEST_CASE(testBlackFormulaForwardDerivative) {

    BOOST_TEST_MESSAGE("Testing forward derivative of the Black formula...");

    std::vector<Real> strikes;
    strikes.push_back(0.1);
    strikes.push_back(0.5);
    strikes.push_back(1.0);
    strikes.push_back(2.0);
    strikes.push_back(3.0);
    const Real vol = 0.1;
    assertBlackFormulaForwardDerivative(Option::Call, strikes, vol);
    assertBlackFormulaForwardDerivative(Option::Put, strikes, vol);
}

BOOST_AUTO_TEST_CASE(testBlackFormulaForwardDerivativeWithZeroStrike) {

    BOOST_TEST_MESSAGE("Testing forward derivative of the Black formula "
        "with zero strike...");

    std::vector<Real> strikes;
    strikes.push_back(0.0);
    const Real vol = 0.1;
    assertBlackFormulaForwardDerivative(Option::Call, strikes, vol);
    assertBlackFormulaForwardDerivative(Option::Put, strikes, vol);
}

BOOST_AUTO_TEST_CASE(testBlackFormulaForwardDerivativeWithZeroVolatility) {

    BOOST_TEST_MESSAGE("Testing forward derivative of the Black formula "
        "with zero volatility...");

    std::vector<Real> strikes;
    strikes.push_back(0.1);
    strikes.push_back(0.5);
    strikes.push_back(1.0);
    strikes.push_back(2.0);
    strikes.push_back(3.0);
    const Real vol = 0.0;
    assertBlackFormulaForwardDerivative(Option::Call, strikes, vol);
    assertBlackFormulaForwardDerivative(Option::Put, strikes, vol);
}

void assertBachelierBlackFormulaForwardDerivative(
    Option::Type optionType,
    const std::vector<Real> &strikes,
    Real bpvol)
{
    Real forward = 1.0;
    Real tte = 10.0;
    Real stdDev = bpvol * std::sqrt(tte);
    Real discount = 0.95;
    Real bump = 0.0001;
    Real epsilon = 1.e-10;
    std::string type = optionType == Option::Call ? "Call" : "Put";

    for (Real strike : strikes) {
        Real delta =
            bachelierBlackFormulaForwardDerivative(optionType, strike, forward, stdDev, discount);
        Real bumpedDelta = bachelierBlackFormulaForwardDerivative(
            optionType, strike, forward + bump, stdDev, discount);

        Real basePremium = bachelierBlackFormula(
            optionType, strike, forward, stdDev, discount);
        Real bumpedPremium = bachelierBlackFormula(
            optionType, strike, forward + bump, stdDev, discount);
        Real deltaApprox = (bumpedPremium - basePremium) / bump;

        /*! Based on the Mean Value Theorem, the below inequality
            should hold for any function that is monotonic in the
            area of the bump.
         */
        bool success =
            (std::max(delta, bumpedDelta) + epsilon > deltaApprox)
            &&  (deltaApprox > std::min(delta, bumpedDelta) - epsilon);

        if (!success)
        {
            BOOST_ERROR("Failed to calculate the derivative of the"
                        " Bachelier Black formula w.r.t. forward"
                        << "\n option type       :" << type
                        << "\n forward           :" << forward
                        << "\n strike            :" << strike
                        << "\n stdDev            :" << stdDev
                        << "\n analytical delta  :" << delta
                        << "\n approximated delta:" << deltaApprox);
        }
    }
}

BOOST_AUTO_TEST_CASE(testBachelierBlackFormulaForwardDerivative) {

    BOOST_TEST_MESSAGE("Testing forward derivative of the "
        "Bachelier Black formula...");

    std::vector<Real> strikes;
    strikes.push_back(-3.0);
    strikes.push_back(-2.0);
    strikes.push_back(-1.0);
    strikes.push_back(-0.5);
    strikes.push_back(0.0);
    strikes.push_back(0.5);
    strikes.push_back(1.0);
    strikes.push_back(2.0);
    strikes.push_back(3.0);
    const Real vol = 0.001;
    assertBachelierBlackFormulaForwardDerivative(Option::Call, strikes, vol);
    assertBachelierBlackFormulaForwardDerivative(Option::Put, strikes, vol);
}

BOOST_AUTO_TEST_CASE(testBachelierBlackFormulaForwardDerivativeWithZeroVolatility) {

    BOOST_TEST_MESSAGE("Testing forward derivative of the Bachelier Black formula "
        "with zero volatility...");

    std::vector<Real> strikes;
    strikes.push_back(-3.0);
    strikes.push_back(-2.0);
    strikes.push_back(-1.0);
    strikes.push_back(-0.5);
    strikes.push_back(0.0);
    strikes.push_back(0.5);
    strikes.push_back(1.0);
    strikes.push_back(2.0);
    strikes.push_back(3.0);
    const Real vol = 0.0;
    assertBachelierBlackFormulaForwardDerivative(Option::Call, strikes, vol);
    assertBachelierBlackFormulaForwardDerivative(Option::Put, strikes, vol);
}

BOOST_AUTO_TEST_CASE(testBachelierImpliedVolAtIntrinsic) {

    BOOST_TEST_MESSAGE("Testing Bachelier implied volatility for prices at intrinsic value...");

    // For a deep in-the-money option the Bachelier price rounds to the
    // intrinsic value, so the time value computed by the inversion can come
    // out a few ulps below zero.  Such a price must still invert instead of
    // being rejected as arbitrageable.
    const Real tte = 0.01;
    const Real vol = 1e-5;

    for (Real forward : {-0.01, 0.005, 0.05, 1.0}) {
        for (Real offset : {0.01, 0.03, 0.07}) {
            for (auto type : {Option::Call, Option::Put}) {
                const Real strike =
                    (type == Option::Call) ? forward - offset : forward + offset;
                const Real price = bachelierBlackFormula(
                    type, strike, forward, vol * std::sqrt(tte));

                Real impliedVol = -1.0;
                BOOST_CHECK_NO_THROW(
                    impliedVol = bachelierBlackFormulaImpliedVol(
                        type, strike, forward, tte, price));
                BOOST_CHECK(impliedVol >= 0.0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(testDeepInTheMoneyPricesNotBelowIntrinsic) {

    BOOST_TEST_MESSAGE("Testing that deep in-the-money Black and Bachelier prices "
                       "are not below their intrinsic value...");

    // With a small enough standard deviation the cumulative normal
    // saturates and the price has to reduce to the discounted intrinsic
    // value exactly; a result a few ulps below it is rejected by the
    // no-arbitrage checks of the implied-volatility inversions.
    struct Case {
        Option::Type type;
        Real strike, forward;
    };
    const Case cases[] = {{Option::Call, 50.0, 100.0},    {Option::Put, 100.0, 50.0},
                          {Option::Call, 0.03, 0.05},     {Option::Put, 0.05, 0.03},
                          {Option::Call, 1000.0, 1234.5}, {Option::Put, 1234.5, 1000.0}};

    for (const Case& c : cases) {
        Real theta = (c.type == Option::Call) ? 1.0 : -1.0;
        Real moneyness = theta * (c.forward - c.strike);
        for (Real discount : {1.0, 0.97}) {
            Real intrinsic = discount * std::max(moneyness, 0.0);

            // both d1 and d2 are above 20
            Real black = blackFormula(c.type, c.strike, c.forward, 0.01, discount);
            if (black < intrinsic)
                BOOST_ERROR("Black price for strike " << c.strike << ", forward " << c.forward
                            << ", discount " << discount << " is below intrinsic by "
                            << intrinsic - black);

            // (forward - strike)/stdDev is 20
            Real bachelier =
                bachelierBlackFormula(c.type, c.strike, c.forward, moneyness / 20.0, discount);
            if (bachelier < intrinsic)
                BOOST_ERROR("Bachelier price for strike " << c.strike << ", forward "
                            << c.forward << ", discount " << discount
                            << " is below intrinsic by " << intrinsic - bachelier);
        }
    }
}

namespace {

    struct LetsBeRationalCase {
        Real price, stdDev, conditioning;
        bool resolvable;
    };

    /* Accurate Black price for a given strike, forward and total standard
       deviation, built from the normalised Black function so that the
       round trip is not limited by cancellation in the Black formula.  The
       conditioning is the factor by which a relative error in the price is
       amplified in the implied standard deviation. */
    LetsBeRationalCase letsBeRationalCase(Option::Type type,
                                          Real strike,
                                          Real forward,
                                          Real stdDev,
                                          Real discount,
                                          Real displacement) {
        const Real k = strike + displacement, f = forward + displacement;
        // accurate log-moneyness of the given strike and forward
        const Real x = (k <= 2.0 * f && f <= 2.0 * k) ? std::log1p((f - k) / k) : std::log(f / k);
        const Real theta = (type == Option::Call) ? 1.0 : -1.0;
        const Real intrinsic = std::max(theta * (f - k), Real(0.0));
        const Real otm =
            std::sqrt(f) * std::sqrt(k) * detail::normalisedBlackCall(-std::fabs(x), stdDev);
        const Real undiscounted = intrinsic + otm;
        const Real price = discount * undiscounted;
        const Real vega = std::sqrt(f) * std::sqrt(k) * detail::normalisedVega(x, stdDev);
        const Real upperBound = (type == Option::Call) ? f : k;
        // subnormal prices, prices rounding to intrinsic or to the upper
        // bound do not determine the standard deviation
        const bool resolvable = otm > QL_MIN_POSITIVE_REAL / QL_EPSILON &&
                                price / discount > intrinsic && undiscounted < upperBound;
        return {price, stdDev, std::max(Real(1.0), undiscounted / (stdDev * vega)), resolvable};
    }

}

BOOST_AUTO_TEST_CASE(testLetsBeRationalNormalisedBlack) {

    BOOST_TEST_MESSAGE("Testing the normalised Black function used by Let's Be Rational...");

    // reference values computed with 50-digit arithmetic, covering the four
    // evaluation branches (asymptotic expansion, small-t expansion, direct
    // formula, erfcx-based formula)
    struct Case {
        Real x, s, expected;
    };
    const Case cases[] = {
        {-30.0, 1.0, 1.4405938487063514186e-199},    {-12.0, 1.0, 1.2910437115451598069e-34},
        {-15.0, 1.2, 2.9623898007235195834e-37},     {-1.0, 0.1, 7.4653976115686172329e-26},
        {-0.001, 0.0001, 7.4745602454210835492e-29}, {0.0, 1e-8, 3.9894228040143268462e-9},
        {-0.1, 2.0, 0.63412712270593895706},         {0.0, 3.0, 0.86638559746228386799},
        {-30.0, 4.0, 2.3782760139263438428e-15},     {-3.0, 1.0, 0.00034299241765225832057},
        {-1.0, 0.5, 0.0041423587506785749349},       {-10.0, 5.0, 0.0041547780138823765422},
        {-5.0, 1.5, 0.00013118738389704313023},      {-4.5, 0.45, 3.2810356800176665841e-25},
        {-20.0, 2.0, 9.152756597844069956e-25},      {0.001, 0.0001, 0.0010000000416666672083},
        {1e-5, 0.1, 0.039882611876465373266}};

    // The error is a few tens of ulps at most where the result is tiny: the
    // exponential of a large argument cannot be more accurate than that, and
    // the small-t expansion has some cancellation in its leading
    // coefficient, as documented in the reference implementation. The
    // plain Black formula would lose several hundred ulps at some of the
    // points above.
    const Real tolerance = 1.0e-14;

    for (const Case& c : cases) {
        const Real b = detail::normalisedBlackCall(c.x, c.s);
        const Real error = std::fabs(b - c.expected) / c.expected;
        if (error > tolerance)
            BOOST_ERROR("failed to reproduce normalised Black price"
                        << std::setprecision(17) << "\n x:          " << c.x
                        << "\n s:          " << c.s << "\n calculated: " << b
                        << "\n expected:   " << c.expected << "\n rel. error: " << error);

        // put-call symmetry b(x) = b(-x) + intrinsic
        const Real symmetric =
            detail::normalisedBlackCall(-c.x, c.s) + (std::exp(0.5 * c.x) - std::exp(-0.5 * c.x));
        if (std::fabs(symmetric - b) >
            1.0e-15 * std::max(Real(1.0), std::exp(0.5 * std::fabs(c.x))))
            BOOST_ERROR("normalised Black function violates put-call parity"
                        << std::setprecision(17) << "\n x: " << c.x << "\n s: " << c.s
                        << "\n b(x):  " << b << "\n b(-x) + intrinsic: " << symmetric);
    }
}

BOOST_AUTO_TEST_CASE(testLetsBeRationalRoundTrip) {

    BOOST_TEST_MESSAGE("Testing Let's Be Rational implied standard deviation "
                       "round trip on a wide grid...");

    const Real forward = 100.0, discount = 0.97;
    const Real logMoneyness[] = {-30.0, -20.0, -10.0, -5.0, -3.0, -1.0, -0.5, -0.1,
                                 -0.01, -1e-3, -1e-5, 0.0,  1e-5, 1e-3, 0.01, 0.1,
                                 0.5,   1.0,   3.0,   5.0,  10.0, 20.0, 30.0};
    const Real stdDevs[] = {1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0};

    // The result is exact up to a few ulps of the price, amplified by the
    // conditioning of the inversion.
    const Real tolerance = 32 * QL_EPSILON;

    Size checked = 0;
    for (auto type : {Option::Call, Option::Put}) {
        for (Real x : logMoneyness) {
            const Real strike = forward * std::exp(-x);
            for (Real stdDev : stdDevs) {
                const LetsBeRationalCase c =
                    letsBeRationalCase(type, strike, forward, stdDev, discount, 0.0);
                if (!c.resolvable)
                    continue;
                ++checked;
                const Real implied = blackFormulaImpliedStdDevLetsBeRational(type, strike, forward,
                                                                             c.price, discount);
                const Real error = std::fabs(implied - stdDev) / stdDev;
                if (error > tolerance * c.conditioning)
                    BOOST_ERROR("failed to recover standard deviation"
                                << std::setprecision(17) << "\n type:         " << type
                                << "\n strike:       " << strike << "\n forward:      " << forward
                                << "\n price:        " << c.price << "\n stdDev:       " << stdDev
                                << "\n implied:      " << implied << "\n rel. error:   " << error
                                << "\n conditioning: " << c.conditioning);
            }
        }
    }
    // make sure that the grid is not mostly skipped
    BOOST_CHECK(checked > 300);

    // at the normalised level, two Householder iterations suffice everywhere:
    // further iterations do not change the result
    for (Real x : logMoneyness) {
        for (Real stdDev : stdDevs) {
            const Real b = detail::normalisedBlackCall(-std::fabs(x), stdDev);
            if (b <= QL_MIN_POSITIVE_REAL / QL_EPSILON)
                continue;
            Size iterations = 0;
            const Real implied = detail::normalisedImpliedStdDevLetsBeRational(
                b, -std::fabs(x), Option::Call, 100, &iterations);
            // without the limit, the iteration stops by itself once the
            // steps are at round-off level
            if (iterations >= 10)
                BOOST_ERROR("no convergence without iteration limit for x = "
                            << -std::fabs(x) << ", stdDev = " << stdDev);
            const Real limited =
                detail::normalisedImpliedStdDevLetsBeRational(b, -std::fabs(x), Option::Call);
            // a third iteration, if any, does not change the result materially
            const Real difference = std::fabs(implied - limited) / implied;
            if (difference > 64 * QL_EPSILON)
                BOOST_ERROR("two iterations are not enough for x = "
                            << -std::fabs(x) << ", stdDev = " << stdDev << ": relative difference "
                            << difference);
        }
    }
}

BOOST_AUTO_TEST_CASE(testLetsBeRationalWithDisplacement) {

    BOOST_TEST_MESSAGE("Testing Let's Be Rational implied standard deviation "
                       "with displacement...");

    struct Case {
        Real strike, forward, displacement;
    };
    // negative rates made admissible by the displacement
    const Case cases[] = {{-0.005, 0.01, 0.03}, {0.01, -0.005, 0.03}, {-0.01, -0.01, 0.02},
                          {0.02, 0.015, 0.0},   {100.0, 80.0, 25.0},  {-0.029, 0.05, 0.03}};
    const Real stdDevs[] = {0.001, 0.05, 0.2, 1.0, 2.5};
    const Real discount = 0.9;

    for (const Case& c : cases) {
        for (auto type : {Option::Call, Option::Put}) {
            for (Real stdDev : stdDevs) {
                // consistency with QuantLib's Black formula...
                const Real price =
                    blackFormula(type, c.strike, c.forward, stdDev, discount, c.displacement);
                const Real implied = blackFormulaImpliedStdDevLetsBeRational(
                    type, c.strike, c.forward, price, discount, c.displacement);
                const Real repriced =
                    blackFormula(type, c.strike, c.forward, implied, discount, c.displacement);
                // the Black formula has an absolute round-off error
                // proportional to the (shifted) forward
                if (std::fabs(repriced - price) > 1.0e-13 * (c.forward + c.displacement))
                    BOOST_ERROR("failed to reprice with displacement"
                                << std::setprecision(17) << "\n type:         " << type
                                << "\n strike:       " << c.strike << "\n forward:      "
                                << c.forward << "\n displacement: " << c.displacement
                                << "\n stdDev:       " << stdDev << "\n implied:      " << implied
                                << "\n price:        " << price << "\n repriced:     " << repriced);

                // ...and the round trip on accurate prices
                const LetsBeRationalCase lbr =
                    letsBeRationalCase(type, c.strike, c.forward, stdDev, discount, c.displacement);
                if (!lbr.resolvable)
                    continue;
                const Real exact = blackFormulaImpliedStdDevLetsBeRational(
                    type, c.strike, c.forward, lbr.price, discount, c.displacement);
                const Real error = std::fabs(exact - stdDev) / stdDev;
                if (error > 32 * QL_EPSILON * lbr.conditioning)
                    BOOST_ERROR("failed to recover standard deviation with displacement"
                                << std::setprecision(17) << "\n type:         " << type
                                << "\n strike:       " << c.strike << "\n forward:      "
                                << c.forward << "\n displacement: " << c.displacement
                                << "\n stdDev:       " << stdDev << "\n implied:      " << exact
                                << "\n rel. error:   " << error);
            }
        }
    }

    // the payoff overload agrees
    const auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Put, -0.005);
    const Real price = blackFormula(payoff, 0.01, 0.3, discount, 0.03);
    QL_CHECK_CLOSE_FRACTION(
        blackFormulaImpliedStdDevLetsBeRational(payoff, 0.01, price, discount, 0.03),
        blackFormulaImpliedStdDevLetsBeRational(Option::Put, -0.005, 0.01, price, discount, 0.03),
        1e-14);
}

BOOST_AUTO_TEST_CASE(testLetsBeRationalEdgeCases) {

    BOOST_TEST_MESSAGE("Testing Let's Be Rational implied standard deviation "
                       "in edge cases...");

    const Real forward = 100.0, discount = 0.9;

    // at the intrinsic value, including zero-priced out-of-the-money options
    BOOST_CHECK_EQUAL(blackFormulaImpliedStdDevLetsBeRational(Option::Call, 80.0, forward,
                                                              discount * 20.0, discount),
                      0.0);
    BOOST_CHECK_EQUAL(blackFormulaImpliedStdDevLetsBeRational(Option::Put, 120.0, forward,
                                                              discount * 20.0, discount),
                      0.0);
    BOOST_CHECK_EQUAL(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 120.0, forward, 0.0, discount), 0.0);
    BOOST_CHECK_EQUAL(
        blackFormulaImpliedStdDevLetsBeRational(Option::Put, forward, forward, 0.0, discount), 0.0);
    // a deep in-the-money price that rounds to the intrinsic value
    const Real deep = blackFormula(Option::Call, 50.0, forward, 0.01, discount);
    BOOST_CHECK_EQUAL(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 50.0, forward, deep, discount), 0.0);

    // no solution below the intrinsic value or at the upper bound
    BOOST_CHECK_THROW(blackFormulaImpliedStdDevLetsBeRational(
                          Option::Call, 80.0, forward, discount * 20.0 * (1.0 - 1e-10), discount),
                      Error);
    BOOST_CHECK_THROW(blackFormulaImpliedStdDevLetsBeRational(Option::Put, 120.0, forward,
                                                              discount * 19.9, discount),
                      Error);
    BOOST_CHECK_THROW(blackFormulaImpliedStdDevLetsBeRational(Option::Call, 80.0, forward,
                                                              discount * forward, discount),
                      Error);
    BOOST_CHECK_THROW(blackFormulaImpliedStdDevLetsBeRational(Option::Put, 80.0, forward,
                                                              discount * 80.0, discount),
                      Error);
    BOOST_CHECK_THROW(blackFormulaImpliedStdDevLetsBeRational(Option::Call, 80.0, forward,
                                                              discount * forward * 1.01, discount),
                      Error);

    // invalid inputs
    BOOST_CHECK_THROW(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 100.0, forward, -1.0, discount),
        Error);
    BOOST_CHECK_THROW(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 100.0, forward, 5.0, 0.0), Error);
    BOOST_CHECK_THROW(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 100.0, forward, 5.0, 1.0, -1.0),
        Error);
    BOOST_CHECK_THROW(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 100.0, -0.01, 5.0, 1.0, 0.01), Error);
    // with a zero shifted strike the call price does not depend on the volatility
    BOOST_CHECK_THROW(
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, -0.01, 0.02, 0.01, 1.0, 0.01), Error);

    // Close to the money with a tiny standard deviation, the log-moneyness
    // must be accurate to full relative precision: log(forward/strike) would
    // give a relative error around 1e-11 here. The prices were computed with
    // 50-digit arithmetic for the given (double) strike.
    struct NearTheMoney {
        Real price, stdDev;
    };
    for (const NearTheMoney& c :
         {NearTheMoney{7.474564014045043e-30, 1e-7}, NearTheMoney{8.331551225004244e-06, 1e-6}}) {
        const Real implied =
            blackFormulaImpliedStdDevLetsBeRational(Option::Call, 100.00010000005, 100.0, c.price);
        QL_CHECK_CLOSE_FRACTION(implied, c.stdDev, 1e-14);
    }

    // close to the upper bound
    for (Real stdDev : {8.0, 12.0, 14.0}) {
        for (auto type : {Option::Call, Option::Put}) {
            const Real strike = 110.0;
            const LetsBeRationalCase c =
                letsBeRationalCase(type, strike, forward, stdDev, discount, 0.0);
            BOOST_REQUIRE(c.resolvable);
            const Real implied =
                blackFormulaImpliedStdDevLetsBeRational(type, strike, forward, c.price, discount);
            const Real error = std::fabs(implied - stdDev) / stdDev;
            if (error > 32 * QL_EPSILON * c.conditioning)
                BOOST_ERROR("failed to recover standard deviation close to the upper bound"
                            << std::setprecision(17) << "\n type:       " << type
                            << "\n stdDev:     " << stdDev << "\n implied:    " << implied
                            << "\n rel. error: " << error);
        }
    }

    // deep out of the money, with prices down to about 1e-196, and tiny
    // standard deviations
    struct Case {
        Option::Type type;
        Real strike, stdDev;
    };
    const Case cases[] = {{Option::Call, forward * std::exp(30.0), 1.0},
                          {Option::Put, forward * std::exp(-30.0), 1.0},
                          {Option::Call, forward * std::exp(12.0), 0.5},
                          {Option::Put, forward * std::exp(-0.5), 0.05},
                          {Option::Call, forward, 1e-8},
                          {Option::Put, forward * (1.0 - 1e-7), 1e-6}};
    for (const Case& c : cases) {
        const LetsBeRationalCase lbr =
            letsBeRationalCase(c.type, c.strike, forward, c.stdDev, discount, 0.0);
        BOOST_REQUIRE(lbr.resolvable);
        const Real implied =
            blackFormulaImpliedStdDevLetsBeRational(c.type, c.strike, forward, lbr.price, discount);
        const Real error = std::fabs(implied - c.stdDev) / c.stdDev;
        if (error > 32 * QL_EPSILON * lbr.conditioning)
            BOOST_ERROR("failed to recover standard deviation"
                        << std::setprecision(17) << "\n type:       " << c.type
                        << "\n strike:     " << c.strike << "\n price:      " << lbr.price
                        << "\n stdDev:     " << c.stdDev << "\n implied:    " << implied
                        << "\n rel. error: " << error);
    }
}

BOOST_AUTO_TEST_CASE(testLetsBeRationalMovesWithInputs) {

    BOOST_TEST_MESSAGE("Testing that the Let's Be Rational implied standard deviation "
                       "responds to its inputs...");

    const Real strike = 105.0, forward = 100.0, discount = 0.95, price = 4.0;
    const Real base =
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, strike, forward, price, discount);

    // higher prices imply higher standard deviations
    Real previous = 0.0;
    for (Real p = 0.5; p < 90.0; p += 0.5) {
        const Real implied =
            blackFormulaImpliedStdDevLetsBeRational(Option::Call, strike, forward, p, discount);
        if (!(implied > previous))
            BOOST_ERROR("implied standard deviation not increasing with price: "
                        << implied << " at price " << p << " after " << previous);
        previous = implied;
    }

    // each input moves the result...
    const Real lowerDiscount =
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, strike, forward, price, 0.9);
    const Real higherStrike =
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, 110.0, forward, price, discount);
    const Real higherForward =
        blackFormulaImpliedStdDevLetsBeRational(Option::Call, strike, 102.0, price, discount);
    const Real displaced = blackFormulaImpliedStdDevLetsBeRational(Option::Call, strike, forward,
                                                                   price, discount, 10.0);
    BOOST_CHECK(lowerDiscount > base); // same price, larger undiscounted value
    BOOST_CHECK(higherStrike > base);  // further out of the money
    BOOST_CHECK(higherForward < base); // closer to the money
    BOOST_CHECK(displaced < base);     // same price on a larger shifted forward

    // ...consistently with the Black formula
    QL_CHECK_CLOSE_FRACTION(blackFormula(Option::Call, strike, forward, lowerDiscount, 0.9), price,
                            1e-12);
    QL_CHECK_CLOSE_FRACTION(blackFormula(Option::Call, 110.0, forward, higherStrike, discount),
                            price, 1e-12);
    QL_CHECK_CLOSE_FRACTION(blackFormula(Option::Call, strike, 102.0, higherForward, discount),
                            price, 1e-12);
    QL_CHECK_CLOSE_FRACTION(blackFormula(Option::Call, strike, forward, displaced, discount, 10.0),
                            price, 1e-12);

    // the put with the same undiscounted time value gives the same result
    const Real putPrice = price - discount * (forward - strike);
    QL_CHECK_CLOSE_FRACTION(
        blackFormulaImpliedStdDevLetsBeRational(Option::Put, strike, forward, putPrice, discount),
        base, 1e-12);
}

BOOST_AUTO_TEST_CASE(testLetsBeRationalVersusNewtonSafe) {

    BOOST_TEST_MESSAGE("Testing that Let's Be Rational is at least as accurate as "
                       "the default Black implied standard deviation...");

    const Real forward = 100.0, discount = 0.97;
    const Real logMoneyness[] = {-5.0, -3.0, -1.0, -0.5, -0.1, -0.01, 0.0,
                                 0.01, 0.1,  0.5,  1.0,  3.0,  5.0};
    const Real stdDevs[] = {1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0};

    Size compared = 0;
    // errors in units of the accuracy allowed by the price, eps * conditioning
    Real maxNormalisedErrorNewton = 0.0, maxNormalisedErrorLbr = 0.0;
    for (auto type : {Option::Call, Option::Put}) {
        for (Real x : logMoneyness) {
            const Real strike = forward * std::exp(-x);
            for (Real stdDev : stdDevs) {
                const LetsBeRationalCase c =
                    letsBeRationalCase(type, strike, forward, stdDev, discount, 0.0);
                if (!c.resolvable)
                    continue;
                Real newton;
                try {
                    newton = blackFormulaImpliedStdDev(type, strike, forward, c.price, discount);
                } catch (Error&) {
                    continue; // only compare where the existing solver converges
                }
                const Real lbr = blackFormulaImpliedStdDevLetsBeRational(type, strike, forward,
                                                                         c.price, discount);
                const Real errorNewton = std::fabs(newton - stdDev) / stdDev;
                const Real errorLbr = std::fabs(lbr - stdDev) / stdDev;
                ++compared;
                maxNormalisedErrorNewton =
                    std::max(maxNormalisedErrorNewton, errorNewton / (QL_EPSILON * c.conditioning));
                maxNormalisedErrorLbr =
                    std::max(maxNormalisedErrorLbr, errorLbr / (QL_EPSILON * c.conditioning));
                // no less accurate, up to what the price resolution allows
                if (errorLbr > std::max(errorNewton, 32 * QL_EPSILON * c.conditioning))
                    BOOST_ERROR("Let's Be Rational less accurate than the default solver"
                                << std::setprecision(17) << "\n type:          " << type
                                << "\n strike:        " << strike << "\n stdDev:        " << stdDev
                                << "\n default:       " << newton << "\n LetsBeRational: " << lbr);
            }
        }
    }
    BOOST_CHECK(compared > 150);
    // relative to the accuracy that the price allows, the new solver is
    // within a few ulps everywhere, and better than the default one
    BOOST_CHECK(maxNormalisedErrorLbr <= 32.0);
    BOOST_CHECK(maxNormalisedErrorLbr < maxNormalisedErrorNewton);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
