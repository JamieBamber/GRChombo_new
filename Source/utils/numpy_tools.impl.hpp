/* Inline functions for c++ versions of some useful 
numpy functions */

#if !defined(NUMPY_TOOLS_HPP_)
#error "This file should only be included through numpy_tools.hpp"
#endif

#ifndef NUMPY_TOOLS_IMPL_HPP_
#define NUMPY_TOOLS_IMPL_HPP_

inline std::vector<double> NumpyTools::linspace( double start, double stop, int N, bool endpoint=true){
                std::vector<double> output;
                output.reserve(N);
                double x;
                for(int i=0; i<N; i++){
                        if (endpoint) {
                                x = start + i*(stop - start)/(N-1);
                        } else {
                                x = start + i*(stop - start)/N;
                        }
                        output.push_back(x);
                }
                return output;
}      

inline std::vector<double> NumpyTools::logspace(double start, double stop, int N, bool endpoint=true){
        std::vector<double> output;
        output.reserve(N);
        double x;
        for(int i=0; i<N; i++){
                if (endpoint) {
                        x = start*exp(i*(log(stop) - log(start))/(N-1));
                } else {
                        x = start*exp(i*(log(stop) - log(start))/(N));
                }
                output.push_back(x);
        }
        return output;
}


#endif /* NUMPY_TOOLS_IMPL_HPP_ */
