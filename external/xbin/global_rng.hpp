#pragma once
/** \file */

#include <random>

namespace xbin {

static std::mt19937& global_rng() {
  std::random_device rd;
  static std::mt19937 rng(rd());
  return rng;
}
}  // namespace rpxdock
