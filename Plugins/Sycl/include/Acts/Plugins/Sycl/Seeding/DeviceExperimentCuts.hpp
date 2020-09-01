// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Plugins/Sycl/Seeding/detail/Types.h"

namespace Acts::Sycl {
/// @class DeviceExperimentCuts can be used to increase or decrease seed weights
///  based on the space points used in a seed. Seed weights are also
/// influenced by the SeedFilter default implementation. This tool is also used
/// to decide if a seed passes a seed weight cut. As the weight is stored in
/// seeds, there are two distinct methods.
class DeviceExperimentCuts {
 public:
  /// Returns seed weight bonus/malus depending on detector considerations.
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return seed weight to be added to the seed's weight
  float seedWeight(const detail::deviceSpacePoint& bottom,
                   const detail::deviceSpacePoint& middle,
                   const detail::deviceSpacePoint& top) const {
    float weight = 0;
    if (bottom.r > 150) {
      weight = 400;
    }
    if (top.r < 150) {
      weight = 200;
    }
    return weight;
  };
  /// @param weight the current seed weight
  /// @param bottom bottom space point of the current seed
  /// @param middle middle space point of the current seed
  /// @param top top space point of the current seed
  /// @return true if the seed should be kept, false if the seed should be
  /// discarded
  bool singleSeedCut(float weight, const detail::deviceSpacePoint& bottom,
                     const detail::deviceSpacePoint& middle,
                     const detail::deviceSpacePoint& top) const {
    return !(bottom.r > 150. && weight < 380.);
  };
};
}  // namespace Acts::Sycl