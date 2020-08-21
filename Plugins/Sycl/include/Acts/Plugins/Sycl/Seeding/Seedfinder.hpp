// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Plugins/Sycl/Utilities/Helpers.h"

inline namespace cl{
  namespace sycl{
    class queue;
  }
};

namespace Acts::Sycl {

void offloadComputations( cl::sycl::queue* q,
                          const offloadSeedfinderConfig& configData,
                          const std::vector<offloadSpacePoint>& bottomSPs,
                          const std::vector<offloadSpacePoint>& middleSPs,
                          const std::vector<offloadSpacePoint>& topSPs,
                          std::vector<std::vector<SeedData>>& seeds);

cl::sycl::queue* createQueue();

template <typename external_spacepoint_t>
class Seedfinder {
  public:
  Seedfinder(Acts::SeedfinderConfig<external_spacepoint_t> config);

  ~Seedfinder() = default;
  Seedfinder() = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t>&) = delete;
  Seedfinder<external_spacepoint_t>& operator=(
    const Seedfinder<external_spacepoint_t>&) = delete;

  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t> > createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
  cl::sycl::queue* m_queue;
};

} // namespace Acts::Sycl

// Include the template implementation.
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.ipp"