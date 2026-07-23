// Reason for edit start: extend the shared prm parser from the Phase 3 moving sphere to the Phase 4 moving-body setup so box dimensions and prescribed rotation can be configured without scattering prm reads through main.cu.
#pragma once

#include <core/all.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace turbine_core {
namespace moving_body {

struct MesaPDConfig
{
    bool enabled{ false };
    uint_t subcycles{ uint_t(1) };
};

struct PSMConfig
{
    bool enabled{ false };
    bool phase1ScaffoldOnly{ true };
};

struct TrajectoryPoint
{
    uint_t step{ uint_t(0) };
    walberla::Vector3< real_t > center{ real_t(0), real_t(0), real_t(0) };
};

// Reason for edit start: Phase 4e adds one piece-level config record so a moving convex-polyhedron body can be assembled from multiple mesh pieces that share one rigid-body motion but still keep separate local offsets.
struct ConvexPieceConfig
{
    std::string meshFile{};
    walberla::Vector3< real_t > offset{ real_t(0), real_t(0), real_t(0) };
};
// Reason for edit end: Phase 4e adds one piece-level config record so a moving convex-polyhedron body can be assembled from multiple mesh pieces that share one rigid-body motion but still keep separate local offsets.

struct MovingBodyConfig
{
    bool enabled{ false };
    std::string mode{ "kinematic" };
    std::string representation{ "sphere" };
    real_t radius{ real_t(1) };
    walberla::Vector3< real_t > boxEdgeLength{ real_t(1), real_t(1), real_t(1) };
    // Reason for edit start: Phase 4d adds one optional mesh-file hook so the convex-polyhedron path can load a user-provided STL/OBJ/OFF mesh while reusing the existing boxEdgeLength control as the target fit envelope.
    std::string meshFile{};
    // Reason for edit end: Phase 4d adds one optional mesh-file hook so the convex-polyhedron path can load a user-provided STL/OBJ/OFF mesh while reusing the existing boxEdgeLength control as the target fit envelope.
    walberla::Vector3< real_t > initialPosition{ real_t(0), real_t(0), real_t(0) };
    walberla::Vector3< real_t > initialVelocity{ real_t(0), real_t(0), real_t(0) };
    walberla::Vector3< real_t > initialRotation{ real_t(0), real_t(0), real_t(0) };
    walberla::Vector3< real_t > angularVelocity{ real_t(0), real_t(0), real_t(0) };
    std::string trajectoryType{ "piecewise_linear" };
    uint_t trajectoryPointCount{ uint_t(0) };
    std::vector<TrajectoryPoint> trajectoryPoints{};
    // Reason for edit start: Phase 4e stores an optional list of convex mesh pieces beside the older single-mesh hook so the app can switch between one-hull and multi-piece rigid assemblies from prm alone.
    uint_t convexPieceCount{ uint_t(0) };
    std::vector<ConvexPieceConfig> convexPieces{};
    // Reason for edit end: Phase 4e stores an optional list of convex mesh pieces beside the older single-mesh hook so the app can switch between one-hull and multi-piece rigid assemblies from prm alone.
};

struct RuntimeConfig
{
    MesaPDConfig mesaPD{};
    PSMConfig psm{};
    MovingBodyConfig body{};

    bool anyFeatureEnabled() const
    {
        return mesaPD.enabled || psm.enabled || body.enabled;
    }
};

inline RuntimeConfig loadRuntimeConfig(const std::shared_ptr< walberla::Config > & globalConfig)
{
    RuntimeConfig config{};

    if (globalConfig->getNumBlocks("MesaPD"))
    {
        const auto mesaPDConfig = globalConfig->getOneBlock("MesaPD");
        config.mesaPD.enabled = mesaPDConfig.getParameter<bool>("enabled", false);
        config.mesaPD.subcycles = mesaPDConfig.getParameter<uint_t>("subcycles", uint_t(1));
    }

    if (globalConfig->getNumBlocks("PSM"))
    {
        const auto psmConfig = globalConfig->getOneBlock("PSM");
        config.psm.enabled = psmConfig.getParameter<bool>("enabled", false);
        config.psm.phase1ScaffoldOnly = psmConfig.getParameter<bool>("phase1ScaffoldOnly", true);
    }

    if (globalConfig->getNumBlocks("MovingBody"))
    {
        const auto movingBodyConfig = globalConfig->getOneBlock("MovingBody");
        // Reason for edit start: Phase 4e centralizes the optional quoted-string normalization so both the legacy single-mesh path and the new per-piece mesh-file path accept the usual prm quoting rules.
        const auto normalizeOptionalQuotedString = [](std::string value) {
            if (value.size() >= size_t(2)) {
                const char first = value.front();
                const char last = value.back();
                if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
                    value = value.substr(size_t(1), value.size() - size_t(2));
                }
            }
            return value;
        };
        // Reason for edit end: Phase 4e centralizes the optional quoted-string normalization so both the legacy single-mesh path and the new per-piece mesh-file path accept the usual prm quoting rules.
        config.body.enabled = movingBodyConfig.getParameter<bool>("enabled", false);
        config.body.mode = movingBodyConfig.getParameter<std::string>("mode", "kinematic");
        config.body.representation = movingBodyConfig.getParameter<std::string>("representation", "sphere");
        config.body.radius = movingBodyConfig.getParameter<real_t>("radius", real_t(1));
        config.body.boxEdgeLength =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("boxEdgeLength",
                                                                          walberla::Vector3<real_t>{ real_t(1), real_t(1), real_t(1) });
        // Reason for edit start: Phase 4d and Phase 4e normalize both the older single `meshFile` setting and the new repeated `piece` mesh-file settings right after parsing so all geometry paths resolve cleanly from standard quoted prm input.
        config.body.meshFile =
                normalizeOptionalQuotedString(movingBodyConfig.getParameter<std::string>("meshFile", ""));
        walberla::Config::Blocks pieceBlocks;
        movingBodyConfig.getBlocks("piece", pieceBlocks);
        for (const auto & pieceBlock : pieceBlocks)
        {
            ConvexPieceConfig piece{};
            piece.meshFile =
                    normalizeOptionalQuotedString(pieceBlock.getParameter<std::string>("meshFile", ""));
            piece.offset =
                    pieceBlock.getParameter<walberla::Vector3<real_t>>(
                            "offset",
                            walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
            WALBERLA_CHECK(!piece.meshFile.empty(),
                           "Phase 4e MovingBody::piece blocks require a non-empty meshFile.")
            config.body.convexPieces.push_back(piece);
        }
        config.body.convexPieceCount = uint_t(config.body.convexPieces.size());
        // Reason for edit end: Phase 4d and Phase 4e normalize both the older single `meshFile` setting and the new repeated `piece` mesh-file settings right after parsing so all geometry paths resolve cleanly from standard quoted prm input.
        config.body.initialPosition =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("initialPosition",
                                                                          walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
        config.body.initialVelocity =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("initialVelocity",
                                                                          walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
        config.body.initialRotation =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("initialRotation",
                                                                          walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
        config.body.angularVelocity =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("angularVelocity",
                                                                          walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
    }

    if (globalConfig->getNumBlocks("Trajectory"))
    {
        const auto trajectoryConfig = globalConfig->getOneBlock("Trajectory");
        config.body.trajectoryType = trajectoryConfig.getParameter<std::string>("type", "piecewise_linear");
        walberla::Config::Blocks pointBlocks;
        trajectoryConfig.getBlocks("point", pointBlocks);
        std::sort(pointBlocks.begin(), pointBlocks.end(),
                  [](const walberla::Config::BlockHandle & lhs, const walberla::Config::BlockHandle & rhs) {
                      return lhs.getParameter<uint_t>("step") < rhs.getParameter<uint_t>("step");
                  });

        for (const auto & pointBlock : pointBlocks)
        {
            TrajectoryPoint point{};
            point.step = pointBlock.getParameter<uint_t>("step");
            point.center = pointBlock.getParameter<walberla::Vector3<real_t>>("center");

            if (!config.body.trajectoryPoints.empty() &&
                config.body.trajectoryPoints.back().step == point.step)
            {
                config.body.trajectoryPoints.back() = point;
            }
            else
            {
                config.body.trajectoryPoints.push_back(point);
            }
        }

        config.body.trajectoryPointCount = uint_t(config.body.trajectoryPoints.size());
    }

    return config;
}

} // namespace moving_body
} // namespace turbine_core
// Reason for edit end: extend the shared prm parser from the Phase 3 moving sphere to the Phase 4 moving-body setup so box dimensions and prescribed rotation can be configured without scattering prm reads through main.cu.
