// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <sstream>
#include <string>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace detail {

  /// @brief a struct that allows you to evaluate the intersection
  ///
  /// It actually modifies the position/direction of the intersection
  /// attempt to allow for a potential better estimate when being
  /// in a magnetic field setup
  struct IntersectionCorrector
  {

    Vector3D startPos   = Vector3D(0., 0., 0.);
    Vector3D startDir   = Vector3D(0., 0., 0.);
    double   pathLength = 0.;

    double straightLineStep = 100 * units::_um;

    IntersectionCorrector(const Vector3D& spos  = Vector3D(0., 0., 0.),
                          const Vector3D& sdir  = Vector3D(0., 0., 0.),
                          double          spath = 0.)
      : startPos(spos), startDir(sdir), pathLength(spath)
    {
    }

    /// Boolean() operator - returns false for void modifier
    explicit operator bool() const { return true; }

    /// empty correction interface
    bool
    operator()(Vector3D& pos, Vector3D& dir, double& path)
    {
      // approximation with straight line
      if (path * path < straightLineStep * straightLineStep) {
        return false;
      }
      // no correction at start positions
      if (pathLength == 0. || pos.isApprox(startPos)
          || dir.isApprox(startDir)) {
        // can not correct: this is the initial step, no reference point
        return false;
      }
      // change of direction per path length unit
      Vector3D deltaDir = (dir - startDir) * 1. / pathLength;
      dir               = dir + 0.5 * path * deltaDir;
      // return true
      return true;
    }
  };

}  // namespace detail
}  // namespace Acts