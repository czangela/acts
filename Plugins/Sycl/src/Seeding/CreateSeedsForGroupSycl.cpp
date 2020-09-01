// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include <CL/sycl.hpp>

namespace Acts::Sycl {
/// Kernel classes in order of execution.
class duplet_search_bottom_kernel;
class duplet_search_top_kernel;
class ind_copy_bottom_kernel;
class ind_copy_top_kernel;
class transform_coord_bottom_kernel;
class transform_coord_top_kernel;
class triplet_search_kernel;
class filter_2sp_fixed_kernel;

uint32_t numWorkItems(uint32_t numThreads, uint32_t workGroupSize) {
  auto q = numThreads / workGroupSize + (numThreads % workGroupSize != 0);
  return q * workGroupSize;
}

void createSeedsForGroupSycl(
    const std::shared_ptr<cl::sycl::queue>& q,
    const detail::deviceSeedfinderConfig& configData,
    const DeviceExperimentCuts& deviceCuts,
    const std::vector<detail::deviceSpacePoint>& bottomSPs,
    const std::vector<detail::deviceSpacePoint>& middleSPs,
    const std::vector<detail::deviceSpacePoint>& topSPs,
    std::vector<std::vector<SeedData>>& seeds) {
  // Each vector stores data of space points in simplified
  // structures of float variables
  // M: number of middle space points
  // B: number of bottom space points
  // T: number of top space points
  const uint32_t M = middleSPs.size();
  const uint32_t B = bottomSPs.size();
  const uint32_t T = topSPs.size();

  // Up to the Nth space point, the sum of compatible bottom/top space points.
  // We need these for indexing other vectors later in the algorithm.
  std::vector<uint32_t> sumBotCompUptoMid(M + 1, 0);
  std::vector<uint32_t> sumTopCompUptoMid(M + 1, 0);
  std::vector<uint32_t> sumBotTopCombined(M + 1, 0);

  // Similarly to indBotCompMid and indTopCompMid, we store the indices of the
  // middle space points of the corresponding edges.
  std::vector<uint32_t> indMidBotComp;
  std::vector<uint32_t> indMidTopComp;

  // Number of edges for middle-bottom and middle-top duplet bipartite graphs.
  uint32_t edgesBottom = 0;
  uint32_t edgesTop = 0;
  uint32_t edgesComb = 0;

  try {
    using am = cl::sycl::access::mode;
    using at = cl::sycl::access::target;

    uint32_t globalBufferSize =
        q->get_device().get_info<cl::sycl::info::device::global_mem_size>();
    uint32_t maxWorkGroupSize =
        q->get_device().get_info<cl::sycl::info::device::max_work_group_size>();

    // Device allocations
    detail::deviceSpacePoint* deviceBottomSPs =
        cl::sycl::malloc_device<detail::deviceSpacePoint>(B, *q);
    detail::deviceSpacePoint* deviceMiddleSPs =
        cl::sycl::malloc_device<detail::deviceSpacePoint>(M, *q);
    detail::deviceSpacePoint* deviceTopSPs =
        cl::sycl::malloc_device<detail::deviceSpacePoint>(T, *q);

    // Store the number of compatible bottom/top space points per middle space
    // point.
    std::atomic_uint32_t* deviceCountBotDuplets =
        cl::sycl::malloc_device<std::atomic_uint32_t>(M, *q);
    std::atomic_uint32_t* deviceCountTopDuplets =
        cl::sycl::malloc_device<std::atomic_uint32_t>(M, *q);

    uint32_t* deviceNumTopDuplets = cl::sycl::malloc_device<uint32_t>(M, *q);

    // The limit of compatible bottom [top] space points per middle space point
    // is B [T]. Temporarily we reserve buffers of this size (M*B and M*T). We
    // store the indices of bottom [top] space points in bottomSPs [topSPs]. We
    // move the indices to optimal size vectors for easier indexing.
    uint32_t* deviceTmpIndBot = cl::sycl::malloc_device<uint32_t>(M * B, *q);
    uint32_t* deviceTmpIndTop = cl::sycl::malloc_device<uint32_t>(M * T, *q);

    q->memcpy(deviceBottomSPs, bottomSPs.data(),
              sizeof(detail::deviceSpacePoint) * (B));
    q->memcpy(deviceMiddleSPs, middleSPs.data(),
              sizeof(detail::deviceSpacePoint) * (M));
    q->memcpy(deviceTopSPs, topSPs.data(),
              sizeof(detail::deviceSpacePoint) * (T));
    q->memset(deviceCountBotDuplets, 0, M * sizeof(std::atomic_uint32_t));
    q->memset(deviceCountTopDuplets, 0, M * sizeof(std::atomic_uint32_t));
    q->wait();

    //*********************************************//
    // ********** DUPLET SEARCH - BEGIN ********** //
    //*********************************************//

    // After completing the duplet search, we'll have successfully contructed
    // two bipartite graphs for bottom-middle and top-middle space points. We
    // store the indices of the bottom/top space points of the edges of the
    // graphs. They index the bottomSPs and topSPs vectors.

    {
      q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<duplet_search_bottom_kernel>(
            cl::sycl::range<2>{M, B}, [=](cl::sycl::id<2> idx) {
              const auto mid = idx[0];
              const auto bot = idx[1];
              // We check whether this thread actually makes sense (within
              // bounds).
              if (mid < M && bot < B) {
                const auto midSP = deviceMiddleSPs[mid];
                const auto botSP = deviceBottomSPs[bot];

                const auto deltaR = midSP.r - botSP.r;
                const auto cotTheta = (midSP.z - botSP.z) / deltaR;
                const auto zOrigin = midSP.z - midSP.r * cotTheta;

                if (!(deltaR < configData.deltaRMin) &&
                    !(deltaR > configData.deltaRMax) &&
                    !(cl::sycl::abs(cotTheta) > configData.cotThetaMax) &&
                    !(zOrigin < configData.collisionRegionMin) &&
                    !(zOrigin > configData.collisionRegionMax)) {
                  // We keep counting duplets with atomic variables
                  const auto ind = deviceCountBotDuplets[mid].fetch_add(1);
                  deviceTmpIndBot[mid * B + ind] = bot;
                }
              }
            });
      });

      q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<duplet_search_top_kernel>(
            cl::sycl::range<2>{M, B}, [=](cl::sycl::id<2> idx) {
              const auto mid = idx[0];
              const auto top = idx[1];
              // We check whether this thread actually makes sense (within
              // bounds).
              if (mid < M && top < T) {
                const auto midSP = deviceMiddleSPs[mid];
                const auto topSP = deviceTopSPs[top];

                const auto deltaR = topSP.r - midSP.r;
                const auto cotTheta = (topSP.z - midSP.z) / deltaR;
                const auto zOrigin = midSP.z - midSP.r * cotTheta;

                if (!(deltaR < configData.deltaRMin) &&
                    !(deltaR > configData.deltaRMax) &&
                    !(cl::sycl::abs(cotTheta) > configData.cotThetaMax) &&
                    !(zOrigin < configData.collisionRegionMin) &&
                    !(zOrigin > configData.collisionRegionMax)) {
                  // We keep counting duplets with atomic access.
                  const auto ind = deviceCountTopDuplets[mid].fetch_add(1);
                  deviceTmpIndTop[mid * T + ind] = top;
                }
              }
            });
      });
    }  // sync

    q->memcpy(deviceNumTopDuplets, deviceCountTopDuplets, M * sizeof(uint32_t));

    //*********************************************//
    // *********** DUPLET SEARCH - END *********** //
    //*********************************************//

    // retrieve results from counting duplets
    {
      std::vector<uint32_t> countBotDuplets(M);
      std::vector<uint32_t> countTopDuplets(M);

      q->memcpy(countBotDuplets.data(), deviceCountBotDuplets,
                M * sizeof(std::atomic_uint32_t));
      q->memcpy(countTopDuplets.data(), deviceCountTopDuplets,
                M * sizeof(std::atomic_uint32_t));
      q->wait();

      for (uint32_t i = 1; i < M + 1; ++i) {
        sumBotCompUptoMid[i] +=
            sumBotCompUptoMid[i - 1] + countBotDuplets[i - 1];
        sumTopCompUptoMid[i] +=
            sumTopCompUptoMid[i - 1] + countTopDuplets[i - 1];
        sumBotTopCombined[i] += sumBotTopCombined[i - 1] +
                                countTopDuplets[i - 1] * countBotDuplets[i - 1];
      }

      edgesBottom = sumBotCompUptoMid[M];
      edgesTop = sumTopCompUptoMid[M];
      edgesComb = sumBotTopCombined[M];

      if (edgesBottom == 0 || edgesTop == 0)
        return;

      indMidBotComp.reserve(edgesBottom);
      indMidTopComp.reserve(edgesTop);

      for (uint32_t mid = 0; mid < M; ++mid) {
        std::fill_n(std::back_inserter(indMidBotComp), countBotDuplets[mid],
                    mid);
        std::fill_n(std::back_inserter(indMidTopComp), countTopDuplets[mid],
                    mid);
      }
    }  // sync

    cl::sycl::nd_range<1> edgesBotNdRange{
        cl::sycl::range<1>(numWorkItems(edgesBottom, maxWorkGroupSize)),
        cl::sycl::range<1>(maxWorkGroupSize)};

    cl::sycl::nd_range<1> edgesTopNdRange{
        cl::sycl::range<1>(numWorkItems(edgesTop, maxWorkGroupSize)),
        cl::sycl::range<1>(maxWorkGroupSize)};

    uint32_t* deviceIndBot = cl::sycl::malloc_device<uint32_t>(edgesBottom, *q);
    uint32_t* deviceIndTop = cl::sycl::malloc_device<uint32_t>(edgesTop, *q);
    uint32_t* deviceMidIndPerBot =
        cl::sycl::malloc_device<uint32_t>(edgesBottom, *q);
    uint32_t* deviceMidIndPerTop =
        cl::sycl::malloc_device<uint32_t>(edgesTop, *q);
    uint32_t* deviceSumBot = cl::sycl::malloc_device<uint32_t>(M + 1, *q);
    uint32_t* deviceSumTop = cl::sycl::malloc_device<uint32_t>(M + 1, *q);
    uint32_t* deviceSumComb = cl::sycl::malloc_device<uint32_t>(M + 1, *q);
    deviceLinEqCircle* deviceLinBot =
        cl::sycl::malloc_device<deviceLinEqCircle>(edgesBottom, *q);
    deviceLinEqCircle* deviceLinTop =
        cl::sycl::malloc_device<deviceLinEqCircle>(edgesTop, *q);

    q->memcpy(deviceMidIndPerBot, indMidBotComp.data(),
              sizeof(uint32_t) * edgesBottom);
    q->memcpy(deviceMidIndPerTop, indMidTopComp.data(),
              sizeof(uint32_t) * edgesTop);
    q->memcpy(deviceSumBot, sumBotCompUptoMid.data(),
              sizeof(uint32_t) * (M + 1));
    q->memcpy(deviceSumTop, sumTopCompUptoMid.data(),
              sizeof(uint32_t) * (M + 1));
    q->memcpy(deviceSumComb, sumBotTopCombined.data(),
              sizeof(uint32_t) * (M + 1));
    q->wait();

    // Copy indices from temporary matrices to final, optimal size vectors.
    {
      q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<ind_copy_bottom_kernel>(
            edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
              auto idx = item.get_global_linear_id();
              if (idx < edgesBottom) {
                uint32_t mid = deviceMidIndPerBot[idx];
                uint32_t ind =
                    deviceTmpIndBot[mid * B + idx - deviceSumBot[mid]];
                deviceIndBot[idx] = ind;
              }
            });
      });

      q->submit([&](cl::sycl::handler& h) {
        h.parallel_for<ind_copy_top_kernel>(
            edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
              auto idx = item.get_global_linear_id();
              if (idx < edgesTop) {
                uint32_t mid = deviceMidIndPerTop[idx];
                uint32_t ind =
                    deviceTmpIndTop[mid * T + idx - deviceSumTop[mid]];
                deviceIndTop[idx] = ind;
              }
            });
      });
    }  // sync

    //************************************************//
    // *** LINEAR EQUATION TRANSFORMATION - BEGIN *** //
    //************************************************//

    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0

    // coordinate transformation middle-bottom pairs
    auto linB = q->submit([&](cl::sycl::handler& h) {
      h.parallel_for<transform_coord_bottom_kernel>(
          edgesBotNdRange, [=](cl::sycl::nd_item<1> item) {
            auto idx = item.get_global_linear_id();
            if (idx < edgesBottom) {
              const auto midSP = deviceMiddleSPs[deviceMidIndPerBot[idx]];
              const auto botSP = deviceBottomSPs[deviceIndBot[idx]];

              const auto xM = midSP.x;
              const auto yM = midSP.y;
              const auto zM = midSP.z;
              const auto rM = midSP.r;
              const auto varianceZM = midSP.varZ;
              const auto varianceRM = midSP.varR;
              const auto cosPhiM = xM / rM;
              const auto sinPhiM = yM / rM;

              const auto deltaX = botSP.x - xM;
              const auto deltaY = botSP.y - yM;
              const auto deltaZ = botSP.z - zM;

              const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
              const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
              const auto iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
              const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
              const auto cot_theta = -(deltaZ * iDeltaR);

              deviceLinEqCircle L;
              L.cotTheta = cot_theta;
              L.zo = zM - rM * cot_theta;
              L.iDeltaR = iDeltaR;
              L.u = x * iDeltaR2;
              L.v = y * iDeltaR2;
              L.er = ((varianceZM + botSP.varZ) +
                      (cot_theta * cot_theta) * (varianceRM + botSP.varR)) *
                     iDeltaR2;

              deviceLinBot[idx] = L;
            }
          });
    });

    // coordinate transformation middle-top pairs
    auto linT = q->submit([&](cl::sycl::handler& h) {
      h.parallel_for<transform_coord_top_kernel>(
          edgesTopNdRange, [=](cl::sycl::nd_item<1> item) {
            auto idx = item.get_global_linear_id();
            if (idx < edgesTop) {
              const auto midSP = deviceMiddleSPs[deviceMidIndPerTop[idx]];
              const auto topSP = deviceTopSPs[deviceIndTop[idx]];

              const auto xM = midSP.x;
              const auto yM = midSP.y;
              const auto zM = midSP.z;
              const auto rM = midSP.r;
              const auto varianceZM = midSP.varZ;
              const auto varianceRM = midSP.varR;
              const auto cosPhiM = xM / rM;
              const auto sinPhiM = yM / rM;

              const auto deltaX = topSP.x - xM;
              const auto deltaY = topSP.y - yM;
              const auto deltaZ = topSP.z - zM;

              const auto x = deltaX * cosPhiM + deltaY * sinPhiM;
              const auto y = deltaY * cosPhiM - deltaX * sinPhiM;
              const auto iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
              const auto iDeltaR = cl::sycl::sqrt(iDeltaR2);
              const auto cot_theta = deltaZ * iDeltaR;

              deviceLinEqCircle L;
              L.cotTheta = cot_theta;
              L.zo = zM - rM * cot_theta;
              L.iDeltaR = iDeltaR;
              L.u = x * iDeltaR2;
              L.v = y * iDeltaR2;
              L.er = ((varianceZM + topSP.varZ) +
                      (cot_theta * cot_theta) * (varianceRM + topSP.varR)) *
                     iDeltaR2;

              deviceLinTop[idx] = L;
            }
          });
    });

    //************************************************//
    // **** LINEAR EQUATION TRANSFORMATION - END **** //
    //************************************************//

    //************************************************//
    // *********** TRIPLET SEARCH - BEGIN *********** //
    //************************************************//

    const auto maxMemoryAllocation = std::min(
        edgesComb, globalBufferSize /
                       uint32_t((sizeof(TripletData) + sizeof(SeedData)) * 2));

    TripletData* deviceCurvImpact =
        cl::sycl::malloc_device<TripletData>(maxMemoryAllocation, (*q));
    SeedData* deviceSeedArray =
        cl::sycl::malloc_device<SeedData>(maxMemoryAllocation, *q);

    seeds.resize(M);
    const float MIN = -100000.f;
    uint32_t lastMiddle = 0;
    for (uint32_t firstMiddle = 0; firstMiddle < M; firstMiddle = lastMiddle) {
      while (lastMiddle + 1 <= M && (sumBotTopCombined[lastMiddle + 1] -
                                         sumBotTopCombined[firstMiddle] <
                                     maxMemoryAllocation)) {
        ++lastMiddle;
      }

      const auto numCombinations =
          sumBotTopCombined[lastMiddle] - sumBotTopCombined[firstMiddle];
      if (numCombinations == 0)
        continue;

      uint32_t* offset = cl::sycl::malloc_device<uint32_t>(1, *q);
      std::atomic_uint32_t* countTriplets =
          cl::sycl::malloc_device<std::atomic_uint32_t>(1, *q);
      std::atomic_uint32_t* countSeeds =
          cl::sycl::malloc_device<std::atomic_uint32_t>(1, *q);
      q->memcpy(offset, &sumBotTopCombined[firstMiddle], sizeof(uint32_t));
      q->memset(countTriplets, 0, sizeof(std::atomic_uint32_t));
      q->memset(countSeeds, 0, sizeof(std::atomic_uint32_t));

      const auto tripletSearchWorkItems =
          numWorkItems(numCombinations, maxWorkGroupSize);
      cl::sycl::nd_range<1> tripletSearchNDRange{
          cl::sycl::range<1>(tripletSearchWorkItems),  // global work items
          cl::sycl::range<1>(maxWorkGroupSize)         // local work items
      };

      auto tripletEvent = q->submit([&](cl::sycl::handler& h) {
        h.depends_on({linB, linT});
        h.parallel_for<triplet_search_kernel>(
            tripletSearchNDRange, [=](cl::sycl::nd_item<1> item) {
              const uint32_t idx = item.get_global_linear_id();
              if (idx < numCombinations) {
                uint32_t L = firstMiddle;
                uint32_t R = lastMiddle;
                uint32_t mid = L;
                while (L < R - 1) {
                  mid = (L + R) / 2;
                  if (idx + offset[0] < deviceSumComb[mid])
                    R = mid;
                  else
                    L = mid;
                }
                mid = L;

                TripletData T = {MIN, MIN};
                deviceCurvImpact[idx] = T;
                const auto numT = deviceNumTopDuplets[mid];

                const auto ib = deviceSumBot[mid] +
                                ((idx - deviceSumComb[mid] + offset[0]) / numT);
                const auto it = deviceSumTop[mid] +
                                ((idx - deviceSumComb[mid] + offset[0]) % numT);

                const auto linBotEq = deviceLinBot[ib];
                const auto linTopEq = deviceLinTop[it];
                const auto midSP = deviceMiddleSPs[mid];

                const auto Vb = linBotEq.v;
                const auto Ub = linBotEq.u;
                const auto Erb = linBotEq.er;
                const auto cotThetab = linBotEq.cotTheta;
                const auto iDeltaRb = linBotEq.iDeltaR;

                const auto Vt = linTopEq.v;
                const auto Ut = linTopEq.u;
                const auto Ert = linTopEq.er;
                const auto cotThetat = linTopEq.cotTheta;
                const auto iDeltaRt = linTopEq.iDeltaR;

                const auto rM = midSP.r;
                const auto varianceRM = midSP.varR;
                const auto varianceZM = midSP.varZ;

                auto iSinTheta2 = (1. + cotThetab * cotThetab);
                auto scatteringInRegion2 =
                    configData.maxScatteringAngle2 * iSinTheta2;
                scatteringInRegion2 *=
                    configData.sigmaScattering * configData.sigmaScattering;
                auto error2 =
                    Ert + Erb +
                    2 * (cotThetab * cotThetat * varianceRM + varianceZM) *
                        iDeltaRb * iDeltaRt;
                auto deltaCotTheta = cotThetab - cotThetat;
                auto deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

                deltaCotTheta = cl::sycl::abs(deltaCotTheta);
                auto error = cl::sycl::sqrt(error2);
                auto dCotThetaMinusError2 =
                    deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
                auto dU = Ut - Ub;

                if ((!(deltaCotTheta2 - error2 > 0) ||
                     !(dCotThetaMinusError2 > scatteringInRegion2)) &&
                    !(dU == 0.)) {
                  auto A = (Vt - Vb) / dU;
                  auto S2 = 1. + A * A;
                  auto B = Vb - A * Ub;
                  auto B2 = B * B;

                  auto iHelixDiameter2 = B2 / S2;
                  auto pT2scatter =
                      4 * iHelixDiameter2 * configData.pT2perRadius;
                  auto p2scatter = pT2scatter * iSinTheta2;
                  auto Im = cl::sycl::abs((A - B * rM) * rM);

                  if (!(S2 < B2 * configData.minHelixDiameter2) &&
                      !((deltaCotTheta2 - error2 > 0) &&
                        (dCotThetaMinusError2 >
                         p2scatter * configData.sigmaScattering *
                             configData.sigmaScattering)) &&
                      !(Im > configData.impactMax)) {
                    countTriplets[0].fetch_add(1);
                    T.curvature = B / std::sqrt(S2);
                    T.impact = Im;
                    deviceCurvImpact[idx] = T;
                  }
                }
              }
            });
      });
      tripletEvent.wait();

      uint32_t sumTriplets;
      auto e0 =
          q->memcpy(&sumTriplets, countTriplets, sizeof(std::atomic_uint32_t));
      e0.wait();

      if (sumTriplets == 0)
        continue;

      q->submit([&](cl::sycl::handler& h) {
         h.depends_on(tripletEvent);
         h.parallel_for<filter_2sp_fixed_kernel>(
             tripletSearchNDRange, [=](cl::sycl::nd_item<1> item) {
               const uint32_t idx = item.get_global_linear_id();
               if (idx < numCombinations &&
                   deviceCurvImpact[idx].curvature != MIN) {
                 uint32_t L = firstMiddle;
                 uint32_t R = lastMiddle;
                 uint32_t mid = L;
                 while (L < R - 1) {
                   mid = (L + R) / 2;
                   if (idx + offset[0] < deviceSumComb[mid])
                     R = mid;
                   else
                     L = mid;
                 }
                 mid = L;

                 const auto sumMidComb = deviceSumComb[mid] - offset[0];
                 const auto idxOffset = idx - sumMidComb;
                 const auto numTopMid = deviceNumTopDuplets[mid];
                 // const auto numTopMid = 20;

                 const auto ib = deviceSumBot[mid] + ((idxOffset) / numTopMid);
                 const auto it = deviceSumTop[mid] + ((idxOffset) % numTopMid);

                 const auto bot = deviceIndBot[ib];
                 const auto top = deviceIndTop[it];

                 const auto current = deviceCurvImpact[idx];

                 // by default compatSeedLimit is 2 -> 2 is hard coded
                 // Variable length arrays are not supported in SYCL kernels.
                 float compatibleSeedR[2];

                 const auto invHelixDiameter = current.curvature;
                 const auto lowerLimitCurv =
                     invHelixDiameter - configData.deltaInvHelixDiameter;
                 const auto upperLimitCurv =
                     invHelixDiameter + configData.deltaInvHelixDiameter;
                 const auto currentTop_r = deviceTopSPs[top].r;
                 auto weight =
                     -(current.impact * configData.impactWeightFactor);

                 uint32_t compatCounter = 0;

                 const auto bottomOffset =
                     ((idxOffset) / numTopMid) * numTopMid;

                 for (uint32_t j = 0;
                      j < numTopMid &&
                      compatCounter < configData.compatSeedLimit;
                      ++j) {
                   uint32_t other_idx = sumMidComb + bottomOffset + j;
                   float otherCurv = deviceCurvImpact[other_idx].curvature;
                   if (otherCurv != MIN && other_idx != idx) {
                     uint32_t other_it = deviceSumTop[mid] + j;
                     float otherTop_r = deviceTopSPs[deviceIndTop[other_it]].r;
                     float deltaR = cl::sycl::abs(currentTop_r - otherTop_r);
                     if (deltaR >= configData.filterDeltaRMin &&
                         otherCurv >= lowerLimitCurv &&
                         otherCurv <= upperLimitCurv) {
                       uint32_t c = 0;
                       for (; c < compatCounter &&
                              cl::sycl::abs(compatibleSeedR[c] - otherTop_r) >=
                                  configData.filterDeltaRMin;
                            ++c) {
                       }
                       if (c == compatCounter) {
                         compatibleSeedR[c] = otherTop_r;
                         ++compatCounter;
                       }
                     }
                   }
                 }

                 weight += compatCounter * configData.compatSeedWeight;

                 const auto bottomSP = deviceBottomSPs[bot];
                 const auto middleSP = deviceMiddleSPs[mid];
                 const auto topSP = deviceTopSPs[top];

                 weight += deviceCuts.seedWeight(bottomSP, middleSP, topSP);

                 if (deviceCuts.singleSeedCut(weight, bottomSP, middleSP,
                                              topSP)) {
                   const auto i = countSeeds[0].fetch_add(1);
                   SeedData D;
                   D.bottom = bot;
                   D.top = top;
                   D.middle = mid;
                   D.weight = weight;
                   deviceSeedArray[i] = D;
                 }
               }
             });
       }).wait();

      uint32_t sumSeeds;
      auto e1 = q->memcpy(&sumSeeds, countSeeds, sizeof(std::atomic_uint32_t));
      e1.wait();

      std::vector<SeedData> hostSeedArray(sumSeeds);
      q->memcpy(&hostSeedArray[0], deviceSeedArray,
                sumSeeds * sizeof(SeedData));
      q->wait();

      for (uint32_t t = 0; t < sumSeeds; ++t) {
        auto m = hostSeedArray[t].middle;
        seeds[m].push_back(hostSeedArray[t]);
      }
    }

    //************************************************//
    // ************ TRIPLET SEARCH - END ************ //
    //************************************************//

    cl::sycl::free(deviceTmpIndBot, *q);
    cl::sycl::free(deviceTmpIndTop, *q);
    cl::sycl::free(deviceBottomSPs, *q);
    cl::sycl::free(deviceMiddleSPs, *q);
    cl::sycl::free(deviceTopSPs, *q);
    cl::sycl::free(deviceCountBotDuplets, *q);
    cl::sycl::free(deviceCountTopDuplets, *q);
    cl::sycl::free(deviceNumTopDuplets, *q);

    cl::sycl::free(deviceLinBot, *q);
    cl::sycl::free(deviceLinTop, *q);

    cl::sycl::free(deviceIndBot, *q);
    cl::sycl::free(deviceIndTop, *q);
    cl::sycl::free(deviceMidIndPerBot, *q);
    cl::sycl::free(deviceMidIndPerTop, *q);
    cl::sycl::free(deviceSumBot, *q);
    cl::sycl::free(deviceSumTop, *q);
    cl::sycl::free(deviceSumComb, *q);

    cl::sycl::free(deviceCurvImpact, *q);
    cl::sycl::free(deviceSeedArray, *q);
  } catch (cl::sycl::exception const& e) {
    std::cout << "Caught synchronous SYCL exception:\n"
              << e.what() << std::endl;
    exit(0);
  }
};
}  // namespace Acts::Sycl