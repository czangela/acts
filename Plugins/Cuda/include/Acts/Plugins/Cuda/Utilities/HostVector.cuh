// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/Arrays.cuh"

// CUDA include(s).
#include <cuda.h>

namespace Acts {
namespace Cuda {

/// Vector holding data in host-pinned memory
template <typename T>
class HostVector {

public:
  /// The variable type being used
  typedef T Variable_t;

  /// Create a vector in host memory
  HostVector(size_t size);

  /// Get a specific element of the vector
  Variable_t& get(size_t offset = 0);
  /// Set a specific element of the vector
  void set(size_t offset, Variable_t val);

  /// Copy memory from a/the device.
  void copyFrom(const Variable_t* devPtr, size_t len, size_t offset);
  /// Copy memory from a/the device asynchronously.
  void copyFrom(const Variable_t* devPtr, size_t len, size_t offset,
                cudaStream_t stream);

  /// Reset the vector to all zeros
  void zeros();

private:
  /// The size of the vector
  size_t m_size;
  /// Smart pointer managing the vector's memory
  host_array< Variable_t > m_array;

}; // class HostVector

} // namespace Cuda
} // namespace Acts

// Include the template implementation.
#include "HostVector.ipp"
