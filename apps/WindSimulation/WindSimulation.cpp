//======================================================================================================================
//
//  WindSimulation – waLBerla-wind example application
//
//  Simulates incompressible wind flow over flat terrain using a D3Q19 SRT-LBM
//  on a Cartesian block-structured domain.  A logarithmic inflow profile is
//  applied at the west boundary; the east boundary uses a simple outflow
//  (zero-gradient) condition; all other faces are periodic.
//
//  Usage:  WindSimulation [<parameter-file>]
//  Default parameter file: WindSimulation.prm
//
//======================================================================================================================

#include "wind/WindField.h"
#include "wind/WindParameters.h"
#include "wind/lbm/WindLBMModel.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/SharedSweep.h"
#include "field/AddToStorage.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"
#include "timeloop/SweepTimeloop.h"

namespace walberla {

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::SRT >;
using PdfField_T     = lbm::PdfField< LatticeModel_T >;
using CommScheme_T   = blockforest::communication::UniformBufferedScheme< LatticeModel_T::CommunicationStencil >;

int main(int argc, char** argv)
{
   // ---- Environment -------------------------------------------------------
   Environment env(argc, argv);

   // ---- Parameters --------------------------------------------------------
   wind::WindParameters params;
   params.inflowDirection   = Vector3< real_t >(real_t(1), real_t(0), real_t(0));
   params.inflowSpeed       = real_t(0.05);
   params.referenceHeight   = real_t(10);
   params.roughnessLength   = real_t(0.1);
   params.vonKarmanConstant = real_t(0.41);

   const uint_t Nx = 64;
   const uint_t Ny = 32;
   const uint_t Nz = 32;

   const uint_t timeSteps = 1000;
   const real_t omega     = real_t(1.8); // SRT relaxation parameter

   // ---- Domain & Block Structure ------------------------------------------
   auto blocks = blockforest::createUniformBlockGrid(
      /* blocks  */ 1, 1, 1,
      /* cells   */ Nx, Ny, Nz,
      /* dx      */ real_t(1),
      /* one block per process */ false,
      /* periodic */ false, true, true);

   // ---- Lattice model -----------------------------------------------------
   LatticeModel_T latticeModel(lbm::collision_model::SRT(omega));

   // ---- PDF field ---------------------------------------------------------
   const BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage(blocks, "pdf", latticeModel,
                                Vector3< real_t >(real_t(0)), real_t(1),
                                uint_c(1), field::fzyx);

   // ---- Wind field --------------------------------------------------------
   auto windField = make_shared< wind::WindField >(blocks, params);
   windField->initialize();

   // ---- LBM-wind coupling -------------------------------------------------
   wind::lbm::WindLBMModel< LatticeModel_T > windLBM(windField, pdfFieldID, latticeModel);
   windLBM.initializePDFs(blocks);

   // ---- Communication scheme ----------------------------------------------
   CommScheme_T commScheme(blocks);
   commScheme.addPackInfo(make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >(pdfFieldID));

   // ---- LBM sweep (stream + collide) -------------------------------------
   // Use the makeSweep helper from lbm::SweepWrappers which does not require
   // an explicit FlagField: all cells are treated as fluid.
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField< uint32_t > >(
      pdfFieldID, /* flagFieldID */ BlockDataID(), Set< FlagUID >{ FlagUID("Fluid") });

   // ---- Time loop ---------------------------------------------------------
   SweepTimeloop timeloop(blocks->getBlockStorage(), timeSteps);

   timeloop.add() << BeforeFunction(commScheme, "Communication")
                  << Sweep(lbm::makeSharedSweep(lbmSweep), "LBM Sweep");

   WALBERLA_LOG_INFO_ON_ROOT("Starting wind simulation: "
                             << Nx << " x " << Ny << " x " << Nz
                             << " cells, " << timeSteps << " time steps.");
   timeloop.run();
   WALBERLA_LOG_INFO_ON_ROOT("Wind simulation finished.");

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv)
{
   return walberla::main(argc, argv);
}
