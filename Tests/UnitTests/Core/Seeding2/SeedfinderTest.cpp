// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

// alternate SeedFilter
#include "Acts/Seeding2/SeedFilter.hpp"
#include "Acts/Seeding2/Seedfinder.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"

template <class config>
void initSeedfinderConfig(config& conf) {
  conf.rMax = 160.;
  conf.deltaRMin = 5.;
  conf.deltaRMax = 160.;
  conf.collisionRegionMin = -250.;
  conf.collisionRegionMax = 250.;
  conf.zMin = -2800.;
  conf.zMax = 2800.;
  conf.maxSeedsPerSpM = 5;
  
  conf.cotThetaMax = 7.40627;
  conf.sigmaScattering = 1.00000;
  conf.minPt = 500.;
  conf.bFieldInZ = 0.00199724;
  conf.beamPos = {-.5, -.5};
  conf.impactMax = 10.;
};

std::vector<const SpacePoint*> readFile(std::string filename) {
  std::string line;
  int layer;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, varianceR, varianceZ;
      if (linetype == "lxyz") {
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        r = std::sqrt(x * x + y * y);
        float f22 = varianceR;
        float wid = varianceZ;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }
        SpacePoint* sp =
            new SpacePoint{x, y, z, r, layer, varianceR, varianceZ};
        //     if(r < 200.){
        //       sp->setClusterList(1,0);
        //     }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int main(int argc, char** argv) {
  std::string file{"sp.txt"};
  bool help(false);
  int groups(500);

  int opt;
  while ((opt = getopt(argc, argv, "hf:g:")) != -1) {
    switch (opt) {
      case 'f':
        file = optarg;
        break;
      case 'g':
        groups = atoi(optarg);
        break;
      case 'h':
        help = true;
        [[fallthrough]];
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-h] [-f FILENAME] [-g groups]\n";
        if (help) {
          std::cout << "      -h : this help" << std::endl;
          std::cout
              << "      -f <FILE> : read spacepoints from FILE." << std::endl;
          std::cout << "      -g <NUM> : limit the number of groups to analyze" << std::endl;
        }
        exit(EXIT_FAILURE);
    }
  }

  std::ifstream f(file);
  if (!f.good()) {
    std::cerr << "input file \"" << file << "\" does not exist\n";
    exit(EXIT_FAILURE);
  }

  auto start_read = std::chrono::system_clock::now();
  std::vector<const SpacePoint*> spVec = readFile(file);
  auto end_read = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_read = end_read - start_read;

  std::cout << "read " << spVec.size() << " SP from file " << file << " in "
            << elapsed_read.count() << "s" << std::endl;

  Acts::SeedfinderConfig<SpacePoint> config;
  initSeedfinderConfig(config);

  Acts::detail::SeedfinderConfig<SpacePoint> config2;
  initSeedfinderConfig(config2);

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());

  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));

  config2.seedFilter = std::make_unique<Acts::detail::SeedFilter<SpacePoint>>(
      Acts::detail::SeedFilter<SpacePoint>(sfconf, &atlasCuts));

  Acts::Seedfinder<SpacePoint> a(config);
  Acts::detail::Seedfinder<SpacePoint> b(config2);  

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVectorA;
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVectorB;
  std::size_t nSeedsA = 0;
  std::size_t nSeedsB = 0;
  
  auto start = std::chrono::system_clock::now();
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  auto group_counter = 0;
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVectorA.push_back(a.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
    ++group_counter;
    if(group_counter == groups){
        break;
    }
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_secondsA = end - start;
  
  auto startB = std::chrono::system_clock::now();
  groupIt = spGroup.begin();
  endOfGroups = spGroup.end();
  group_counter = 0;
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVectorB.push_back(b.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
    ++group_counter;
    if(group_counter == groups){
        break;
    }
  }
  auto endB = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_secondsB = endB - startB;

  for (const auto& seeds : seedVectorA) {
    nSeedsA += seeds.size();
  }
  for (const auto& seeds : seedVectorB) {
    nSeedsB += seeds.size();
  }

  std::cout << "analyzed " << group_counter << " groups for both versions\n";


    // Count the total number of reconstructed seeds.
    std::size_t nMatch = 0;
    double matchPercentage = 0.0;

    assert(seedVectorA.size() == seedVectorB.size());
    for (size_t i = 0; i < seedVectorA.size(); i++) {
    // Access the seeds for this region.
    const auto& seeds_in_host_region = seedVectorA[i];
    const auto& seeds_in_device_region = seedVectorB[i];
    // Loop over all seeds found on the host.
    for (const auto& host_seed : seeds_in_host_region) {
        assert(host_seed.sp().size() == 3);
        // Try to find a matching seed that was found on the accelerator.
        for (const auto& device_seed : seeds_in_device_region) {
        assert(device_seed.sp().size() == 3);
        if ((*(host_seed.sp()[0]) == *(device_seed.sp()[0])) &&
            (*(host_seed.sp()[1]) == *(device_seed.sp()[1])) &&
            (*(host_seed.sp()[2]) == *(device_seed.sp()[2]))) {
            ++nMatch;
            break;
        }
        }
    }
    }
    matchPercentage = (100.0 * nMatch) / nSeedsA;  

  std::cout << std::endl;
  std::cout << "-------------------------- Results --------------------------------"
            << std::endl;
  std::cout << "|          |    Original    |    Alternate    | Speedup/agreement |"
            << std::endl;
  std::cout << "-------------------------------------------------------------------"
            << std::endl;
  std::cout << "| Time [s] |" << std::setw(12) << elapsed_secondsA.count() << std::setw(18);
  std::cout << elapsed_secondsB.count() << std::setw(23);
  std::cout << (elapsed_secondsA.count() / elapsed_secondsB.count()) << std::endl;
  std::cout << "|  Seeds   |" << std::setw(12) << nSeedsA << std::setw(18);
  std::cout << nSeedsB << std::setw(23);
  std::cout << matchPercentage << std::endl;

  return 0;
}
