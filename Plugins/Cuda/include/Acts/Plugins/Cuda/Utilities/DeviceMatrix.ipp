// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Utilities/CudaUtils.cu"

// CUDA include(s).
#include <cuda.h>

// System include(s).
#include <cassert>
#include <cstring>

namespace Acts {
namespace Cuda {

template <typename T>
DeviceMatrix<T>::DeviceMatrix(size_t nRows, size_t nCols)
: m_nRows(nRows), m_nCols(nCols),
  m_array(make_device_array<T>(nRows * nCols)) {

}

template <typename T>
typename DeviceMatrix<T>::Variable_t&
DeviceMatrix<T>::get(size_t row, size_t col) {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template <typename T>
const typename DeviceMatrix<T>::Variable_t&
DeviceMatrix<T>::get(size_t row, size_t col) const {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get()[row + col * m_nRows];
}

template <typename T>
typename DeviceMatrix<T>::Variable_t*
DeviceMatrix<T>::getPtr(size_t row, size_t col) {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template <typename T>
const typename DeviceMatrix<T>::Variable_t*
DeviceMatrix<T>::getPtr(size_t row, size_t col) const {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Return the requested element.
  return m_array.get() + row + col * m_nRows;
}

template <typename T>
void DeviceMatrix<T>::set(size_t row, size_t col, Variable_t val) {

  // Some security check(s).
  assert(row < m_nRows);
  assert(col < m_nCols);

  // Set the requested element.
  m_array.get()[row + col * m_nRows] = val;
  return;
}

template <typename T>
void DeviceMatrix<T>::copyFrom(const Variable_t* hostPtr, size_t len,
                               size_t offset) {

  // Some security check(s).
  assert(offset + len < m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpy(m_array.get() + offset, hostPtr,
                                   len * sizeof(Variable_t),
                                   cudaMemcpyHostToDevice));
  return;
}

template <typename T>
void DeviceMatrix<T>::copyFrom(const Variable_t* hostPtr, size_t len,
                               size_t offset, cudaStream_t stream) {

  // Some security check(s).
  assert(offset + len < m_nRows * m_nCols);

  // Do the copy.
  ACTS_CUDA_ERROR_CHECK(cudaMemcpyAsync(m_array.get() + offset, hostPtr,
                                        len * sizeof(Variable_t),
                                        cudaMemcpyHostToDevice, stream));
  return;
}

template <typename T>
void DeviceMatrix<T>::zeros() {

  memset(m_array.get(), 0, m_nRows * m_nCols * sizeof(Variable_t));
  return;
}

} // namespace Cuda
} // namespace Acts
