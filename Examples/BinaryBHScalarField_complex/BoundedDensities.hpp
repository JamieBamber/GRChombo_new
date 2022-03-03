/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BOUNDEDDENSITIES_HPP_
#define BOUNDEDDENSITIES_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "SimulationParameters.hpp"

//! Calculates the density rho with type matter_t and writes it to the grid
class BoundedDensities
{
  protected:

    //! Params for integration
    const integration_params_t m_params;
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
  public:
    BoundedDensities(integration_params_t a_params, double a_dx, std::array<double, CH_SPACEDIM> a_center)
        : m_dx(a_dx), m_center(a_center), m_params(a_params) 
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
	// get the metric vars from the background
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
	// coordinates
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;    
	data_t R = coords.get_radius();

	// get rho and rho_azimuth
	data_t rho, rho_azimuth, C_t, C_phi;
	rho = current_cell.load_vars(c_rho);
	rho_azimuth = current_cell.load_vars(c_rho_azimuth);
	C_t = current_cell.load_vars(c_C_t);
	C_phi = current_cell.load_vars(c_C_phi);

	//
	// data_t inside = simd_compare_lt(r,m_params.max_integration_radius)*simd_compare_gt(r,m_params.min_integration_radius);
	data_t inside = (R < m_params.max_integration_radius) && (R > m_params.min_integration_radius);
	
        // assign values of density in output box
        current_cell.store_vars(inside*rho, c_rho);
        current_cell.store_vars(inside*rho_azimuth, c_rho_azimuth);
        current_cell.store_vars(inside*C_t, c_C_t);
        current_cell.store_vars(inside*C_phi, c_C_phi);
    }
};

#endif /* BOUNDEDDENSITIES_HPP_ */
