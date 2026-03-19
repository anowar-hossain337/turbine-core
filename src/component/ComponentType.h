
#pragma once

#ifndef TURBINECORE_COMPONENTTYPE_H
#define TURBINECORE_COMPONENTTYPE_H

#include <type_traits>

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace component {

        enum class ComponentType : unsigned int {
            NONE       = 0u,
            BASE       = 1u << 0u,
            TOWER      = 1u << 1u,
            NACELLE    = 1u << 2u,
            HUB        = 1u << 3u,
            BLADE      = 1u << 4u,
            ROTOR_DISK = 1u << 5u,
            ANY        = ~NONE
        };

        using Base_T = std::underlying_type_t <ComponentType>;

        HOST_DEVICE_PREFIX inline constexpr ComponentType operator~ (ComponentType a) {
            return static_cast<ComponentType>(~static_cast<Base_T>(a));
        }

        HOST_DEVICE_PREFIX inline constexpr ComponentType operator| (ComponentType a, ComponentType b) {
            return static_cast<ComponentType>(static_cast<Base_T>(a) | static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr ComponentType operator& (ComponentType a, ComponentType b) {
            return static_cast<ComponentType>(static_cast<Base_T>(a) & static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr ComponentType operator^ (ComponentType a, ComponentType b) {
            return static_cast<ComponentType>(static_cast<Base_T>(a) ^ static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr bool operator!(ComponentType a) {
            return (static_cast<Base_T>(a) == static_cast<Base_T>(ComponentType::NONE));
        }

        HOST_DEVICE_PREFIX inline constexpr ComponentType& operator|=(ComponentType & a, ComponentType b) {
            return a = static_cast<ComponentType>(static_cast<Base_T>(a) | static_cast<Base_T>(b));
        }

        HOST_DEVICE_PREFIX inline constexpr ComponentType& operator&=(ComponentType & a, ComponentType b) {
            return a = static_cast<ComponentType>(static_cast<Base_T>(a) & static_cast<Base_T>(b));
        }

    }

    using component::ComponentType;

}

#endif //TURBINECORE_COMPONENTTYPE_H
