/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DENSITYANDMOM_HPP_
#define DENSITYANDMOM_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "BSSNVars.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "DimensionDefinitions.hpp"

//! Calculates the density rho with type matter_t and writes it to the grid
template <class matter_t> 
class DensityAndMom
{
  public:
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    /// CCZ4 variables
    template <class data_t> 
    using MetricVars = BSSNVars::VarsWithGauge<data_t>;

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

    DensityAndMom(matter_t a_matter, double a_dx, std::array<double, CH_SPACEDIM> a_center, double a_final_a)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx), m_center(a_center), m_final_a(a_final_a)
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
	Tensor<2, data_t> gamma_UU;
	FOR2(i, j){ gamma_UU[i][j]= h_UU[i][j]*(1.0/vars.chi); }
	
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
        rho = rho;

	// conserved j^i = -sqrt(-g)T^i_0 = sqrt(det_gamma))*alpha*gamma^ij[ alpha * S_j - beta^k S_kj ] - beta^i rho for the cartesian coordinates
        Tensor<1, data_t> Sbeta; 
        FOR2(i, j){ Sbeta[i] += vars.shift[j]*emtensor.Sij[i][j]; }                      
        Tensor<1, data_t> J; 
        // conserved 3-current linear momentum vector in the cartesian coordinates 
	FOR2(i, j){ J[i] += sqrt(det_gamma)*vars.lapse*( gamma_UU[i][j]*(vars.lapse*emtensor.Si[j] - Sbeta[j]) ); } 
	FOR1(i){ J[i] += - vars.shift[i]*rho*sqrt(det_gamma); }

	// conserved rho_azimuth = |gamma|(x * S_y - y * S_z)
        data_t rho_azimuth = (x * emtensor.Si[1] - y * emtensor.Si[0]);
	
	// azimuthal momentum covector j_i_phi 
	Tensor<1, data_t> J_azimuth_co;
	FOR2(i,j){ J_azimuth_co[i] += sqrt(det_gamma)*vars.lapse*(emtensor.Sij[i][j]*dxdaz[j]); } 	
	// azimuthal momentum vector j^i_phi
	Tensor<1, data_t> J_azimuth;
        FOR2(i,j){ J_azimuth[i] += (gamma_UU[i][j] - vars.shift[i]*vars.shift[j]/(vars.lapse*vars.lapse))*J_azimuth_co[j]; }

	// **** projections of momentum vectors in the spherical coordinate directions
	data_t J_R = 0;
	data_t J_azimuth_R = 0;
	FOR1(i){ J_R += J[i]*NR[i]; }
	FOR1(i){ J_azimuth_R += J_azimuth[i]*NR[i]; }

	// Compute source term
	

    	// assign values of density in output box
    	current_cell.store_vars(rho, c_rho);
    	current_cell.store_vars(rho_azimuth, c_rho_azimuth);
	current_cell.store_vars(J_R, c_J_R);
	current_cell.store_vars(J_azimuth_R, c_J_azimuth_R);
     }

  protected:
    const matter_t m_matter;              //!< The matter object
    const FourthOrderDerivatives m_deriv; //!< An object for calculating derivatives of the variables
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx, m_final_a;
};

#endif /* DENSITYANDMOM_HPP_ */
