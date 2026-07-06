// Reason for edit start: extend the shared prm parser from the Phase 2 stationary setup to a Phase 3 kinematic trajectory setup while keeping the legacy ABL path optional.
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

struct MovingBodyConfig
{
    bool enabled{ false };
    std::string mode{ "kinematic" };
    std::string representation{ "sphere" };
    real_t radius{ real_t(1) };
    walberla::Vector3< real_t > initialPosition{ real_t(0), real_t(0), real_t(0) };
    walberla::Vector3< real_t > initialVelocity{ real_t(0), real_t(0), real_t(0) };
    std::string trajectoryType{ "piecewise_linear" };
    uint_t trajectoryPointCount{ uint_t(0) };
    std::vector<TrajectoryPoint> trajectoryPoints{};
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
        config.body.enabled = movingBodyConfig.getParameter<bool>("enabled", false);
        config.body.mode = movingBodyConfig.getParameter<std::string>("mode", "kinematic");
        config.body.representation = movingBodyConfig.getParameter<std::string>("representation", "sphere");
        config.body.radius = movingBodyConfig.getParameter<real_t>("radius", real_t(1));
        config.body.initialPosition =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("initialPosition",
                                                                          walberla::Vector3<real_t>{ real_t(0), real_t(0), real_t(0) });
        config.body.initialVelocity =
                movingBodyConfig.getParameter<walberla::Vector3<real_t>>("initialVelocity",
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
// Reason for edit end: extend the shared prm parser from the Phase 2 stationary setup to a Phase 3 kinematic trajectory setup while keeping the legacy ABL path optional.
