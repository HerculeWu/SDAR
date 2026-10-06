#pragma once
// Common includes of the unit tests: the library headers expect the user to define ASSERT
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#define ASSERT(x) do { if (!(x)) throw std::logic_error("ASSERT failed: " #x); } while (0)

#include "Common/Float.h"
