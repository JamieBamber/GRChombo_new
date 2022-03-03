 /* this calculates the associated legendre polynomials with normalisation 
coefficients appropriate for the spherical harmonics. I use the algorithm
described in https://arxiv.org/pdf/1410.1748.pdf */ 

#ifndef ASSOC_LEGENDRE_
#define ASSOC_LEGENDRE_
#include "math.h"

namespace AssocLegendre
{
template <class data_t> struct P_lm
{
    data_t Magnitude;
    data_t Derivative;
};

// calculate the associated legendre polynomials with coefficients suitable for the spherical harmonics
// c = cos(theta)
// s = sin(theta) 
template <class data_t> 
data_t assoc_legendre(const int l, const int m_, const data_t c, const data_t s) {
	// return P^l_m
	const int m = abs(m_);
	data_t P_00 = 1;
	if (l == 0){
		return P_00;
	};
	// Use the P^m-1_m-1 --> P^m_m reccurance relation
	data_t Pm_m = P_00;
	for(int n = 1; n <= m; ++n){
		Pm_m = -1 * Pm_m * sqrt(1 + 1.0/(2*n)) * s;
	};
	if (l == m){ 
		return Pm_m;
	}
	// Use P^m_m --> P^m_m+1 
	data_t Pm_m_plus_1 = sqrt(2*m + 3)*c*Pm_m;
	if (l == m+1){
		return Pm_m_plus_1;
	} else {
		// Use P^m_l-2, P^m_l-1 --> P^m_l	
		data_t Pm_l;
		data_t Pm_l_minus_2 = Pm_m;
		data_t Pm_l_minus_1 = Pm_m_plus_1;
		double a;
		double b;
		for(int L = m+2; L<= l; ++L){
			a = sqrt((double)(4*L*L - 1)/(L*L - m*m));
			b = -1 * sqrt((double)((L-1)*(L-1) - m*m)/(4*(L-1)*(L-1) - 1));
			Pm_l = a*(c*Pm_l_minus_1 + b*Pm_l_minus_2);
			Pm_l_minus_2 = Pm_l_minus_1;
			Pm_l_minus_1 = Pm_l;
		}
		return Pm_l;
	}
}

// define / set my version of the function to calculate the associated legendre polynomials
template <class data_t> 
P_lm<data_t> assoc_legendre_with_deriv (const int l, const int m_, const data_t c, const data_t s) {
	// return (P^l_m, d/d_theta(P^l_m) )
	int m = abs(m_);
	
	P_lm<data_t> P_lm;
	data_t P_00 = 1;

	if (l == 0) {
		P_lm.Magnitude = P_00;
		P_lm.Derivative = 0;
		return P_lm;	
	};
	// Use the P^m-1_m-1 --> P^m_m reccurance relation
	data_t Pm_m = P_00;
	data_t Pm_minus1_m_minus1 = P_00;
	data_t Pm_m_prime = 0;
	for(int n = 1; n <= m; ++n){
		Pm_m = -Pm_minus1_m_minus1 * sqrt(1 + 1.0/(2*n)) * s;
		Pm_m_prime = -sqrt(1 + 1.0/(2*n)) * (c * Pm_minus1_m_minus1 + s * Pm_m_prime);
		Pm_minus1_m_minus1 = Pm_m;
	}
	if (l == m){
		P_lm.Magnitude = Pm_m;
		P_lm.Derivative = Pm_m_prime; 
		return P_lm;
	};
	// Use P^m_m --> P^m_m+1 
	data_t Pm_m_plus_1 = sqrt(2*m + 3)*(c*Pm_m);
	data_t Pm_m_plus_1_prime = sqrt(2*m + 3)*(c*Pm_m_prime - s*Pm_m);
	if (l == m+1) {
		P_lm.Magnitude = Pm_m_plus_1;
		P_lm.Derivative = Pm_m_plus_1_prime;
		return P_lm;
	} else {
		// Use P^m_l-2, P^m_l-1 --> P^m_l	
		data_t Pm_l;
		data_t Pm_l_minus_2 = Pm_m;
		data_t Pm_l_minus_1 = Pm_m_plus_1;
		data_t Pm_l_prime;
                data_t Pm_l_minus_2_prime = Pm_m_prime;
                data_t Pm_l_minus_1_prime = Pm_m_plus_1_prime;
		double a;
		double b;
		for(int L = m+2; L<= l; ++L){
			a = sqrt((double)(4*L*L - 1)/(L*L - m*m));
                        b = -sqrt((double)((L-1)*(L-1) - m*m)/(4*(L-1)*(L-1) - 1));
                        Pm_l = a*(c*Pm_l_minus_1 + b*Pm_l_minus_2);
			Pm_l_prime = a*(c*Pm_l_minus_1_prime -s*Pm_l_minus_1 + b*Pm_l_minus_2_prime);
                        Pm_l_minus_2 = Pm_l_minus_1;
                        Pm_l_minus_1 = Pm_l;
			//
			Pm_l_minus_2_prime = Pm_l_minus_1_prime;
                        Pm_l_minus_1_prime = Pm_l_prime;
		}
		P_lm.Magnitude = Pm_l;
		P_lm.Derivative = Pm_l_prime;
		return P_lm;
	}
}
} // namespace AssocLegendre

#endif // ** ASSOC_LEGENDRE_ ** //	
