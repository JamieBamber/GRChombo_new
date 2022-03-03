/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PHIEXTRACTION_HPP_
#define PHIEXTRACTION_HPP_

#include "SphericalExtraction.hpp"

// pout()
// #include "parstream.H"

//!  The class allows extraction of the scalar field phi in
//!  spherical shells at specified radii, and integration over those shells
/*!
    The values at each theta, phi point may then be
    written to an output file, or integrated across the surfaces.
*/
class PhiExtraction : public SphericalExtraction
{
  private:
	string m_data_subdir, m_suffix;
	std::vector<std::string> var_names = {"phi_Re", "phi_Im", "Pi_Re", "Pi_Im"};
	int m_num_vars = 4;
  public:
    //! The constructor
    PhiExtraction(SphericalExtraction::params_t &a_params, 
                        string a_data_subdir, string a_suffix, 
                    double a_dt,
                    double a_time, bool a_first_step,
                    double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time), 
                                 m_data_subdir(a_data_subdir), m_suffix(a_suffix)
    {
	add_var(c_phi_Re, VariableType::evolution);
	add_var(c_phi_Im, VariableType::evolution);
	add_var(c_Pi_Re, VariableType::evolution);
	add_var(c_Pi_Im, VariableType::evolution);
    }

//! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    PhiExtraction(SphericalExtraction::params_t a_params, string a_data_subdir, string a_suffix,
                double a_dt, double a_time, double a_restart_time = 0.0)
        : PhiExtraction(a_params, a_data_subdir, a_suffix, a_dt, a_time, (a_dt == a_time), 
                          a_restart_time)
    {
    }

//! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
     	// extract the values of the phi field on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
        {
            write_extraction("PhiExtractionOut_");
        }

	// now calculate and write the requested spherical harmonic modes
        std::vector<std::vector<double>>
            vars_mode_integrals(m_num_modes*m_num_vars);

	// const auto vars = current_cell.template load_vars<MatterVars>();

	// slightly pointless step
        for (int ivar =0; ivar < m_num_vars; ++ivar){
		const SphericalExtraction::real_function_t var_func = [ivar](std::vector<double> vars,
                                           	double r, double, double) {
            	// here the std::vector<double> passed will just have
            	// the real phi as its only component
            	return vars[ivar];
        	};	
	
		// add the modes that will be integrated
        	for (int imode = 0; imode < m_num_modes; ++imode)
        	{
            	const auto &mode = m_modes[imode];
            	constexpr int es = 0;
            	add_real_mode_integrand(es, mode.first, mode.second,
                               	var_func, vars_mode_integrals[ivar*m_num_modes+imode]);
        	}
        }
	// do the integration over the surface
        integrate();

	// write the integrals
        for (int ivar =0; ivar < m_num_vars; ++ivar){
		for (int imode = 0; imode < m_num_modes; ++imode)
        	{
            	const auto &mode = m_modes[imode];
            	std::string integrals_filename = "Vars_integral_" + m_data_subdir + "_" + var_names[ivar] + "_lm_" +
                                             	std::to_string(mode.first) +
                                             	std::to_string(mode.second) + m_suffix;
            	std::vector<std::vector<double>> integrals_for_writing = {
                	std::move(vars_mode_integrals[ivar*m_num_modes+imode])};
            	std::vector<std::string> labels = {var_names[ivar] + " integral"};
            	write_integrals(integrals_filename, integrals_for_writing, labels);
        	}
	}
    }
};

#endif /* PHIEXTRACTION_HPP_ */
