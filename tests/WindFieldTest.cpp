//======================================================================================================================
//
//  Unit tests for walberla::wind::WindField
//
//======================================================================================================================

#include "wind/WindField.h"
#include "wind/WindParameters.h"

#include "blockforest/Initialization.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Vector3.h"

using namespace walberla;
using namespace walberla::wind;

// ---------------------------------------------------------------------------
// Helper: create a small uniform block grid for testing.
// ---------------------------------------------------------------------------
static shared_ptr< StructuredBlockStorage > createTestDomain()
{
   return blockforest::createUniformBlockGrid(
      /* blocks */ 1, 1, 1,
      /* cells  */ 8, 8, 8,
      /* dx     */ real_t(1),
      /* one block per process */ false,
      /* periodic */ false, false, false);
}

// ---------------------------------------------------------------------------
// Test: WindField can be constructed and initialized without error.
// ---------------------------------------------------------------------------
static void testConstructionAndInit()
{
   auto blocks = createTestDomain();
   WindParameters params;
   WindField wf(blocks, params);
   wf.initialize(); // must not throw / abort
}

// ---------------------------------------------------------------------------
// Test: After initialization, every cell has a non-negative wind speed and
//       velocity is aligned with the inflow direction.
// ---------------------------------------------------------------------------
static void testInitialVelocityDirection()
{
   auto blocks = createTestDomain();

   WindParameters params;
   params.inflowDirection   = Vector3< real_t >(real_t(1), real_t(0), real_t(0));
   params.inflowSpeed       = real_t(0.05);
   params.referenceHeight   = real_t(4);
   params.roughnessLength   = real_t(0.1);

   WindField wf(blocks, params);
   wf.initialize();

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      // Sample a cell well above the roughness length.
      const Vector3< real_t > vel = wf.getVelocity(*blockIt, 0, 0, 4);

      // The y and z components must be zero (pure x-direction inflow).
      WALBERLA_CHECK_FLOAT_EQUAL(vel[1], real_t(0), "y-component must be zero.");
      WALBERLA_CHECK_FLOAT_EQUAL(vel[2], real_t(0), "z-component must be zero.");

      // The x component must be positive for cells above the roughness length.
      WALBERLA_CHECK_GREATER_EQUAL(vel[0], real_t(0), "Wind speed must be non-negative.");
   }
}

// ---------------------------------------------------------------------------
// Test: getVelocity / setVelocity round-trip.
// ---------------------------------------------------------------------------
static void testGetSetVelocity()
{
   auto blocks = createTestDomain();
   WindParameters params;
   WindField wf(blocks, params);
   wf.initialize();

   const Vector3< real_t > testVel(real_t(0.1), real_t(0.2), real_t(0.3));

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      wf.setVelocity(*blockIt, 3, 3, 3, testVel);
      const Vector3< real_t > retrieved = wf.getVelocity(*blockIt, 3, 3, 3);

      WALBERLA_CHECK_FLOAT_EQUAL(retrieved[0], testVel[0]);
      WALBERLA_CHECK_FLOAT_EQUAL(retrieved[1], testVel[1]);
      WALBERLA_CHECK_FLOAT_EQUAL(retrieved[2], testVel[2]);
   }
}

// ---------------------------------------------------------------------------
// Test: update() does not corrupt the field (steady-state implementation).
// ---------------------------------------------------------------------------
static void testUpdateDoesNotCorrupt()
{
   auto blocks = createTestDomain();
   WindParameters params;
   WindField wf(blocks, params);
   wf.initialize();

   // Record velocity at a representative cell before update.
   Vector3< real_t > velBefore;
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      velBefore = wf.getVelocity(*blockIt, 2, 2, 4);

   wf.update(real_t(1.0));

   // The steady-state implementation should leave the field unchanged.
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      const Vector3< real_t > velAfter = wf.getVelocity(*blockIt, 2, 2, 4);
      WALBERLA_CHECK_FLOAT_EQUAL(velAfter[0], velBefore[0]);
      WALBERLA_CHECK_FLOAT_EQUAL(velAfter[1], velBefore[1]);
      WALBERLA_CHECK_FLOAT_EQUAL(velAfter[2], velBefore[2]);
   }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char** argv)
{
   debug::enterTestMode();
   walberla::Environment env(argc, argv);

   testConstructionAndInit();
   testInitialVelocityDirection();
   testGetSetVelocity();
   testUpdateDoesNotCorrupt();

   return EXIT_SUCCESS;
}
