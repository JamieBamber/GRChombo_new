/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes

#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "BoostedBH.hpp"
#include "ComplexScalarPotential.hpp"
#include "ScalarRotatingCloud.hpp"

// numpy tools
#include "numpy_tools.hpp"

struct integration_params_t
{
    double min_integration_radius = 0;
    double max_integration_radius;
    // int &num_integration_radii;
  
    // std::vector<double> &integration_radii;
  
    std::array<double, CH_SPACEDIM> integration_center;
    std::vector<int> integration_levels;
    int min_integration_level;
};

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
	pout() << "@@@@@ start SimulationParameters::readParams(pp)" << endl;
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp)
    {
    	// Initial and SF data
        pp.load("G_Newton", G_Newton, 0.0);
        pp.load("scalar_mass", potential_params.scalar_mass);
        pp.load("field_amplitude", initial_params.field_amplitude);
	pp.load("scalar_omega", initial_params.omega);
	pp.load("scalar_kappa", initial_params.kappa);
	//
	pp.load("scalar_l", initial_params.l);
        pp.load("scalar_m", initial_params.m);
	pp.load("scalar_center", initial_params.center, center);
        pp.load("alignment", initial_params.alignment);
        pp.load("phase", initial_params.phase);
	
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);
	pp.load("input_time", input_time);
	pp.load("omega_BH", omega_BH);
	// pp.load("final_a", final_a); // final dimensionfull spin J/M of the merged black hole

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        std::array<double, CH_SPACEDIM> centerA, centerB;
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);

	// allow for rotation of input file
        bh1_params.center[0] = centerA[0] + offsetA[0] * cos(omega_BH * input_time) - offsetA[1] * sin(omega_BH * input_time);
        bh1_params.center[1] = centerA[1] + offsetA[1] * cos(omega_BH * input_time) + offsetA[0] * sin(omega_BH * input_time);
        bh1_params.center[2] = centerA[2] + offsetA[2];
        bh2_params.center[0] = centerB[0] + offsetB[0] * cos(omega_BH * input_time) - offsetB[1] * sin(omega_BH * input_time);
        bh2_params.center[1] = centerB[1] + offsetB[1] * cos(omega_BH * input_time) + offsetB[0] * sin(omega_BH * input_time);
        bh2_params.center[2] = centerB[2] + offsetB[2];

	// Time allowed for the field to evolve without the metric evolving
	pp.load("delay", delay, 0.0);

        // Do we want Weyl extraction and puncture tracking?
        pp.load("activate_Weyl_extraction", activate_Weyl_extraction, false);
        pp.load("activate_flux_extraction", activate_flux_extraction, false);
        pp.load("activate_integral", activate_integral, false);

        pp.load("track_punctures", track_punctures, false);
	pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
	pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);
	pp.load("calculate_constraints", calculate_constraints,
                false);

        // hard code num punctures to 2 for now
        int num_punctures = 2;
        initial_puncture_coords.resize(num_punctures);
        initial_puncture_coords[0] = bh1_params.center;
        initial_puncture_coords[1] = bh2_params.center;
	pout() << "@@@@@ Finished SimulationParameters: readParams(pp)" << endl;

	// ** Bounded Densities volume integration parameters **
	pp.load("inner_r", integration_params.min_integration_radius, 0.0);
        pp.load("outer_r", integration_params.max_integration_radius, L/2.0);

	// ** Params for second surface integration **
        pp.get("num_extraction_radii_2", extraction_params_2.num_extraction_radii);
        // -- make integration radius array
        if (pp.contains("min_integration_radius_2") && pp.contains("max_integration_radius_2")) {
                pp.load("linear_or_log", linear_or_log, true);
                pp.load("min_integration_radius_2", min_integration_radius_2);
                pp.load("max_integration_radius_2", max_integration_radius_2);
                if (linear_or_log) {
                        extraction_params_2.extraction_radii = 
                        NumpyTools::linspace(min_integration_radius_2, max_integration_radius_2, extraction_params_2.num_extraction_radii);
                } else {
                        extraction_params.extraction_radii = 
                        NumpyTools::logspace(min_integration_radius_2, max_integration_radius_2, extraction_params_2.num_extraction_radii);
                }       
                pout() << "extraction_params_2.extraction_radii = " << std::endl;
                        for(std::vector<double>::const_iterator i = extraction_params_2.extraction_radii.begin(); i != extraction_params_2.extraction_radii.end(); ++i){
                                pout() << std::to_string(*i) << std::endl;
                        }
                pout() << "end of radii list" << std::endl;
        }
	if (pp.contains("integration_radius_2"))
        {
            pp.load("integration_radius_2", extraction_params_2.extraction_radii, 1,
                    0.1);
        }
	if (pp.contains("extraction_radii_2"))
        {
             	pp.load("extraction_radii_2", extraction_params_2.extraction_radii,
                        extraction_params_2.num_extraction_radii);
        }
        
       	pp.load("num_points_phi_2", extraction_params_2.num_points_phi, 2);
        pp.load("num_points_theta_2", extraction_params_2.num_points_theta, 5);
        if (extraction_params_2.num_points_theta % 2 == 0)
        {
            extraction_params_2.num_points_theta += 1;
            pout() << "Parameter: num_points_theta incompatible with Simpson's "
                   << "rule so increased by 1.\n";
        }
	pp.load("extraction_center_2", extraction_params_2.center, center);
	pp.load("write_extraction_2", extraction_params_2.write_extraction, false);

	pp.load("integral_filename", integral_filename);
    }

    // Initial data
    bool activate_Weyl_extraction, activate_flux_extraction, activate_integral, track_punctures, calculate_constraint_norms, calculate_constraints;
    int puncture_tracking_level;
    std::vector<std::array<double, CH_SPACEDIM>> initial_puncture_coords;
    double G_Newton;
    double delay;
    double omega_BH, input_time;
    bool linear_or_log;
    double min_integration_radius_2, max_integration_radius_2;
    // double final_a;
    ScalarRotatingCloud::params_t initial_params;
    ComplexScalarPotential::params_t potential_params;
    std::string integral_filename;
    integration_params_t integration_params;
    SphericalExtraction::params_t extraction_params_2;
    
    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

