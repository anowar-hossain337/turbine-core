//======================================================================================================================
//
//  Unit tests for walberla::wind::WindParameters
//
//======================================================================================================================

#include "wind/WindParameters.h"

#include "core/debug/TestSubsystem.h"
#include "core/math/Math.h"

#include <cmath>

using namespace walberla;
using namespace walberla::wind;

// ---------------------------------------------------------------------------
// Test: logLawSpeed returns the reference speed at the reference height.
// ---------------------------------------------------------------------------
static void testLogLawAtReferenceHeight()
{
   WindParameters p;
   p.inflowSpeed     = real_t(1.0);
   p.referenceHeight = real_t(10.0);
   p.roughnessLength = real_t(0.1);

   const real_t speed = p.logLawSpeed(p.referenceHeight);
   WALBERLA_CHECK_FLOAT_EQUAL(speed, p.inflowSpeed,
      "logLawSpeed at reference height must equal inflowSpeed.");
}

// ---------------------------------------------------------------------------
// Test: logLawSpeed increases with height.
// ---------------------------------------------------------------------------
static void testLogLawMonotonicity()
{
   WindParameters p;
   p.inflowSpeed     = real_t(1.0);
   p.referenceHeight = real_t(10.0);
   p.roughnessLength = real_t(0.1);

   const real_t u5  = p.logLawSpeed(real_t(5.0));
   const real_t u10 = p.logLawSpeed(real_t(10.0));
   const real_t u20 = p.logLawSpeed(real_t(20.0));

   WALBERLA_CHECK_LESS(u5, u10,  "logLawSpeed must increase with height (5 < 10).");
   WALBERLA_CHECK_LESS(u10, u20, "logLawSpeed must increase with height (10 < 20).");
}

// ---------------------------------------------------------------------------
// Test: logLawSpeed value is computed correctly (manual calculation).
// ---------------------------------------------------------------------------
static void testLogLawValue()
{
   WindParameters p;
   p.inflowSpeed     = real_t(2.0);
   p.referenceHeight = real_t(10.0);
   p.roughnessLength = real_t(0.5);

   const real_t z        = real_t(5.0);
   const real_t expected = real_t(2.0)
                         * (std::log(z / p.roughnessLength)
                            / std::log(p.referenceHeight / p.roughnessLength));

   const real_t speed = p.logLawSpeed(z);
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(speed, expected, real_t(1e-10),
      "logLawSpeed value mismatch.");
}

// ---------------------------------------------------------------------------
// Test: default parameters are physically reasonable.
// ---------------------------------------------------------------------------
static void testDefaultParameters()
{
   WindParameters p;

   WALBERLA_CHECK_GREATER(p.inflowSpeed,       real_t(0));
   WALBERLA_CHECK_GREATER(p.referenceHeight,   real_t(0));
   WALBERLA_CHECK_GREATER(p.roughnessLength,   real_t(0));
   WALBERLA_CHECK_GREATER(p.vonKarmanConstant, real_t(0));

   // The inflow direction should be a non-zero vector.
   WALBERLA_CHECK_GREATER(p.inflowDirection.length(), real_t(0));
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment env(argc, argv);

   testDefaultParameters();
   testLogLawAtReferenceHeight();
   testLogLawMonotonicity();
   testLogLawValue();

   return EXIT_SUCCESS;
}
