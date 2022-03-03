/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FORCEEXTRACTION_HPP_
#define FORCEEXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//!  The class allows extraction of the values of the force components on
//!  Spherical shells at specified radii, and integration over those shells
/*!
   The class allows the user to extract data from the grid for the force
   components over Spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class ForceExtraction : public SphericalExtraction

{
  private:
	std::string m_filename;

  public:
    //! The constructor
    ForceExtraction(SphericalExtraction::params_t &a_params, double a_dt,
                    double a_time, bool a_first_step, std::string a_filename,
                    double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                               a_restart_time), m_filename(a_filename)
    {
        add_var(c_J_R, VariableType::diagnostic);
        add_var(c_J_azimuth_R, VariableType::diagnostic);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    ForceExtraction(SphericalExtraction::params_t a_params, double a_dt,
                    double a_time, std::string a_filename, double a_restart_time = 0.0)
        : ForceExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                          a_filename, a_restart_time)
    {
    }

    // the references of the vars as used in the integrator
    enum M_VARS
    {
        m_J_R,
        m_J_azimuth_R
    };

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the Force scalars on the spheres
        extract(a_interpolator);

        // this would write out the values at every point on the sphere
        if (m_params.write_extraction)
        {
            write_extraction("Force4ExtractionOut_");
        }

        // Setup to integrate J_azimuth_R and J_R
        std::vector<std::vector<double>> force_integrals(2);
        add_var_integrand(m_J_R, force_integrals[m_J_R],
                          IntegrationMethod::simpson);
        add_var_integrand(m_J_azimuth_R, force_integrals[m_J_azimuth_R],
                          IntegrationMethod::simpson);

	pout() << "m_params.extraction_radii[0] = " << m_params.surface_param_values[0];
	pout() << "m_params.num_extraction_radii = " << m_params.num_surfaces;
	 
        // do the integration over the surface
        integrate();

        // write the integrals
        std::vector<std::string> labels(2);
        labels[m_J_R] = "J_R";
        labels[m_J_azimuth_R] = "J_azimuth_R";
        // std::string filename = "Force_integrals";
        write_integrals(m_filename, force_integrals, labels);
    }
};

#endif /* FORCEEXTRACTION_HPP_ */
