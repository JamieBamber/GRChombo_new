/* C++ implementation of some useful numpy functions
*/

#ifndef NUMPY_TOOLS_HPP_
#define NUMPY_TOOLS_HPP_

#include "math.h"
#include <vector>

namespace NumpyTools
{
// make a 1D array of equally spaced values
std::vector<double> linspace( double start, double stop, int N, bool endpoint);

// make a 1D array of equally logarithmically spaced values
std::vector<double> logspace(double start, double stop, int N, bool endpoint);
}

#include "numpy_tools.impl.hpp"

#endif /* NUMPY_TOOLS_HPP_ */
