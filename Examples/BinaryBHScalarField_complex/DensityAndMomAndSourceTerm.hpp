/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DENSITYANDMOM_HPP_
#define DENSITYANDMOM_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "CCZ4Vars.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "DimensionDefinitions.hpp"

#include "CCZ4RHS.hpp"

//! Calculates the density rho with type matter_t and writes it to the grid
template <class matter_t> 
class DensityAndMomAndSourceTerm
{
  public:
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    /// CCZ4 variables
    template <class data_t> 
    using MetricVars = CCZ4Vars::VarsWithGauge<data_t>;

    // Inherit the variables from MatterVars and MetricVars
    template <class data_t>
    struct Vars : public MetricVars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            MetricVars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    DensityAndMomAndSourceTerm(matter_t a_matter, double a_dx, std::array<double, CH_SPACEDIM> a_center, double a_final_a,
	CCZ4_params_t<> a_ccz4_params)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx), m_center(a_center), m_final_a(a_final_a), m_ccz4_params(a_ccz4_params)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
	CH_TIME("DensityAndMom::compute");
        // copy data from chombo gridpoint into local variables, and derivs
	const auto vars = current_cell.template load_vars<Vars>();
        const auto d1 = m_deriv.template diff1<Vars>(current_cell);
        const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
	const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    	// Energy Momentum Tensor
    	const auto emtensor = m_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

	// spacial metric
	const data_t det_gamma = 1.0/(vars.chi*vars.chi*vars.chi);
	Tensor<2, data_t> gamma_UU, gamma, g_UU;
	FOR2(i, j){ gamma[i][j]= vars.h[i][j]/vars.chi; }
	FOR2(i, j){ gamma_UU[i][j]= h_UU[i][j]*(vars.chi); }
	FOR2(i, j){ g_UU[i][j] = gamma_UU[i][j] - vars.shift[i]*vars.shift[j]/(vars.lapse*vars.lapse); }	

	// cartesian coordinates
	const Coordinates<data_t> coords(current_cell, m_dx, m_center);
	data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
	data_t R = coords.get_radius();
	// calculate approximate Kerr Schild radius for the final black hole
	/*double a = m_final_a;
        data_t disc = simd_max((R*R - a*a), 0.001);
        data_t r2 = disc/2 + sqrt(disc*disc/4 + (a*z)*(a*z));
        data_t r = sqrt(r2);*/	
	
	// dx/daz vector
        Tensor<1, data_t> dxdaz;
        dxdaz[0] = - y;
        dxdaz[1] =   x;
        dxdaz[2] = 0;
	// radial direction normal vector
	Tensor <1, data_t> NR;
	NR[0] = x/R;
	NR[1] = y/R;
	NR[2] = z/R;
	// dx/dtheta vector
	data_t rxy = sqrt(x*x + y*y);
	Tensor <1, data_t> dxdtheta;
        dxdtheta[0] = x*z/rxy;
        dxdtheta[1] = y*z/rxy;
        dxdtheta[2] = -rxy;	

	// conserved rho = -sqrt(-g)T^0_0 = sqrt(det_gamma)*(alpha*rho_3+1 - beta^i * S_i)
        data_t rho = vars.lapse*emtensor.rho;
        FOR1(k){ rho += -vars.shift[k]*emtensor.Si[k];   }
        rho = sqrt(det_gamma)*rho;

	// conserved j_0^i = -sqrt(-g)T^i_0 
	//                 = sqrt(det_gamma))*alpha*g^ij( alpha * S_j - beta^k S_kj ) - beta^i rho 
	//                 = - sqrt(det_gamma)) [ alpha rho beta^i - S_j ( beta^i beta^j + alpha^2 gamma^{ij} ) + alpha gamma^{ij} beta^k S_jk ]
        Tensor<1, data_t> J;
	FOR1(i){ J[i] = vars.lapse * emtensor.rho * vars.shift[i]; }
	FOR2(i,j) { J[i] += - emtensor.Si[j] * ( vars.shift[i] * vars.shift[j] + vars.lapse*vars.lapse * gamma_UU[i][j] ); }
	FOR3(i,j,k) { J[i] += vars.lapse * gamma_UU[i][j] * vars.shift[k] * emtensor.Sij[j][k]; }
	FOR1(i){ J[i] = -sqrt(det_gamma) * J[i]; } 

	// conserved rho_azimuth = |gamma|(x * S_y - y * S_z)
        data_t rho_azimuth = sqrt(det_gamma)*(x * emtensor.Si[1] - y * emtensor.Si[0]);
	
	// azimuthal momentum covector j_i_phi 
	Tensor<1, data_t> J_azimuth_co;
	FOR2(i,j){ J_azimuth_co[i] += sqrt(det_gamma)*vars.lapse*(emtensor.Sij[i][j]*dxdaz[j]); } 	
	// azimuthal momentum vector j^i_phi
	Tensor<1, data_t> J_azimuth;
        FOR2(i,j){ J_azimuth[i] += g_UU[i][j] * J_azimuth_co[j]; }

	// projections of momentum vectors in the spherical coordinate directions
	data_t J_R = 0;
	data_t J_azimuth_R = 0;
	FOR1(i){ J_R += J[i]*NR[i]; }
	FOR1(i){ J_azimuth_R += J_azimuth[i]*NR[i]; }

	// Compute source term
	// compute time derivatives of alpha, beta^i, gamma_ij (could also potentially just get this from rhs calculation)

	// something to stop 1/chi blowing up
	data_t chi_reg = simd_max(vars.chi, 0.00000001);

	// dt_alpha determined by the CCZ4 gauge choice
	data_t dt_alpha = m_ccz4_params.lapse_coeff * pow(vars.lapse, m_ccz4_params.lapse_power) * (vars.K - 2 * vars.Theta);	
	FOR1(i){ dt_alpha +=  m_ccz4_params.lapse_advec_coeff * vars.shift[i]*d1.lapse[i]; }	

	// dt_beta^i
	Tensor<1, data_t> dt_shift;
	FOR1(i){ 
		dt_shift[i] = m_ccz4_params.shift_Gamma_coeff * vars.B[i];	
		FOR1(j){ dt_shift[i] += m_ccz4_params.shift_advec_coeff * vars.shift[j] * d1.shift[i][j]; }
	}
	// dt_gamma_ij
	Tensor<2, data_t> dt_gamma_ij, Di_betaj;
	FOR2(i,j){
		Di_betaj[i][j] = 0;
		FOR1(k){ Di_betaj[i][j] += gamma[j][k]*d1.shift[k][i] + chris.LLL[j][i][k]*vars.shift[k]; }
	}
	data_t Kij;
	FOR2(i,j){
		Kij = (vars.A[i][j] + (1.0/3.0)*vars.K*vars.h[i][j])/vars.chi;
		dt_gamma_ij[i][j] = - 2 * vars.lapse * Kij + Di_betaj[i][j] + Di_betaj[j][i];	
	}
	//
	data_t S_t;
	S_t = - emtensor.rho * dt_alpha;
	FOR1(i){ S_t += emtensor.Si[i] * dt_shift[i]; }
	Tensor<2, data_t> S_UUij;
	FOR2(i,j){
		S_UUij[i][j] = 0;
		FOR2(k,l){ S_UUij[i][j] += gamma_UU[i][k]*gamma_UU[j][l]*emtensor.Sij[k][l]; }
	}
	FOR2(i,j){ S_t += (vars.lapse / 2.0) * S_UUij[i][j] * dt_gamma_ij[i][j]; }
	data_t C_t = - sqrt(det_gamma) * S_t;
	//
	Tensor<1, data_t> S_i;
	FOR1(i){
		S_i[i] = - emtensor.rho * d1.lapse[i];
		FOR1(j){ S_i[i] += emtensor.Si[j] * d1.shift[j][i]; }
		FOR2(j,k) { S_i[i] += vars.lapse * S_UUij[j][k] * chris.LLL[j][i][k]; }
	}
	// Now compute C_phi
	data_t C_phi = 0;
	FOR1(i){
		C_phi += sqrt(det_gamma) * dxdaz[i] * S_i[i]; 
        }

    	// assign values of density in output box
    	current_cell.store_vars(rho, c_rho);
    	current_cell.store_vars(rho_azimuth, c_rho_azimuth);
	current_cell.store_vars(J_R, c_J_R);
	current_cell.store_vars(J_azimuth_R, c_J_azimuth_R);
        current_cell.store_vars(C_t, c_C_t);
	current_cell.store_vars(C_phi, c_C_phi); 
    }

  protected:
    const matter_t m_matter;              //!< The matter object
    const FourthOrderDerivatives m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx, m_final_a;
    const CCZ4_params_t<> m_ccz4_params;
};

#endif /* DENSITYANDMOM_HPP_ */
