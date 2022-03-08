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
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"
#ifdef USE_TWOPUNCTURES
#include "TP_Parameters.hpp"
#endif

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
#ifdef USE_TWOPUNCTURES
        read_tp_params(pp);
#else
        read_bh_params(pp);
#endif
        read_shared_params(pp);
    }

#ifdef USE_TWOPUNCTURES
    void read_tp_params(GRParmParse &pp)
    {
        tp_params.verbose = (verbosity > 0);
        // check whether to calculate the target ADM masses or use provided bare
        // masses
        bool calculate_target_masses;
        pp.load("TP_calculate_target_masses", calculate_target_masses, false);
        tp_params.give_bare_mass = !calculate_target_masses;

        // masses
        if (calculate_target_masses)
        {
            pp.load("TP_target_mass_plus", tp_params.target_M_plus);
            pp.load("TP_target_mass_minus", tp_params.target_M_minus);
            pp.load("TP_adm_tol", tp_params.adm_tol, 1e-10);
            pout() << "The black holes have target ADM masses of "
                   << tp_params.target_M_plus << " and "
                   << tp_params.target_M_minus << "\n";
            bh1_params.mass = tp_params.target_M_minus;
            bh2_params.mass = tp_params.target_M_plus;
        }
        else
        {
            pp.load("TP_mass_plus", tp_params.par_m_plus);
            pp.load("TP_mass_minus", tp_params.par_m_minus);
            bh1_params.mass = tp_params.par_m_plus;
            bh2_params.mass = tp_params.par_m_minus;
            pout() << "The black holes have bare masses of "
                   << std::setprecision(16) << tp_params.par_m_plus << " and "
                   << tp_params.par_m_minus << "\n";
            // reset precision
            pout() << std::setprecision(6);
        }

        // BH spin and momenta
        std::array<double, CH_SPACEDIM> spin_minus, spin_plus;
        pp.load("TP_momentum_minus", bh1_params.momentum);
        pp.load("TP_momentum_plus", bh2_params.momentum);
        pp.load("TP_spin_plus", spin_plus);
        pp.load("TP_spin_minus", spin_minus);
        FOR(i)
        {
            tp_params.par_P_minus[i] = bh1_params.momentum[i];
            tp_params.par_P_plus[i] = bh2_params.momentum[i];
            tp_params.par_S_minus[i] = spin_minus[i];
            tp_params.par_S_plus[i] = spin_plus[i];
        }

        pout() << "The corresponding momenta are:";
        pout() << "\nP_plus = ";
        FOR(i) { pout() << tp_params.par_P_plus[i] << " "; }
        pout() << "\nP_minus = ";
        FOR(i) { pout() << tp_params.par_P_minus[i] << " "; }

        pout() << "\nThe corresponding spins are:";
        pout() << "\nS_plus = ";
        FOR(i) { pout() << tp_params.par_S_plus[i] << " "; }
        pout() << "\nS_minus = ";
        FOR(i) { pout() << tp_params.par_S_minus[i] << " "; }
        pout() << "\n";

        // interpolation type
        bool use_spectral_interpolation;
        pp.load("TP_use_spectral_interpolation", use_spectral_interpolation,
                false);
        tp_params.grid_setup_method =
            (use_spectral_interpolation) ? "evaluation" : "Taylor expansion";

        // initial_lapse (default to psi^n)
        pp.load("TP_initial_lapse", tp_params.initial_lapse,
                std::string("psi^n"));
        if (tp_params.initial_lapse != "twopunctures-antisymmetric" &&
            tp_params.initial_lapse != "twopunctures-averaged" &&
            tp_params.initial_lapse != "psi^n" &&
            tp_params.initial_lapse != "brownsville")
        {
            std::string message = "Parameter: TP_initial_lapse: ";
            message += tp_params.initial_lapse;
            message += " invalid";
            MayDay::Error(message.c_str());
        }
        if (tp_params.initial_lapse == "psi^n")
        {
            pp.load("TP_initial_lapse_psi_exponent",
                    tp_params.initial_lapse_psi_exponent, -2.0);
        }

        // Spectral grid parameters
        pp.load("TP_npoints_A", tp_params.npoints_A, 30);
        pp.load("TP_npoints_B", tp_params.npoints_B, 30);
        pp.load("TP_npoints_phi", tp_params.npoints_phi, 16);
        if (tp_params.npoints_phi % 4 != 0)
        {
            MayDay::Error("TP_npoints_phi must be a multiple of 4");
        }

        // Solver parameters and tolerances
        pp.load("TP_Newton_tol", tp_params.Newton_tol, 1e-10);
        pp.load("TP_Newton_maxit", tp_params.Newton_maxit, 5);
        pp.load("TP_epsilon", tp_params.TP_epsilon, 1e-6);
        pp.load("TP_Tiny", tp_params.TP_Tiny, 0.0);
        pp.load("TP_Extend_Radius", tp_params.TP_Extend_Radius, 0.0);

        // BH positions
        pp.load("TP_offset_minus", tp_offset_minus);
        pp.load("TP_offset_plus", tp_offset_plus);
        bh1_params.center = center;
        bh2_params.center = center;
        bh1_params.center[0] += tp_offset_minus;
        bh2_params.center[0] += tp_offset_plus;
        double center_offset_x = 0.5 * (tp_offset_plus + tp_offset_minus);
        tp_params.center_offset[0] = center_offset_x;
        // par_b is half the distance between BH_minus and BH_plus
        tp_params.par_b = 0.5 * (tp_offset_plus - tp_offset_minus);
        pp.load("TP_swap_xz", tp_params.swap_xz, false);

        // Debug output
        pp.load("TP_do_residuum_debug_output",
                tp_params.do_residuum_debug_output, false);
        pp.load("TP_do_initial_debug_output", tp_params.do_initial_debug_output,
                false);

        // Irrelevant parameters set to default value
        tp_params.keep_u_around = false;
        tp_params.use_sources = false;
        tp_params.rescale_sources = true;
        tp_params.use_external_initial_guess = false;
        tp_params.multiply_old_lapse = false;
        tp_params.schedule_in_ADMBase_InitialData = true;
        tp_params.solve_momentum_constraint = false;
        tp_params.metric_type = "something else";
        tp_params.conformal_storage = "not conformal at all";
        tp_params.conformal_state = 0;
        tp_params.mp = 0;
        tp_params.mm = 0;
        tp_params.mp_adm = 0;
        tp_params.mm_adm = 0;
    }
#else
    /// Read BH parameters if not using two punctures
    void read_bh_params(GRParmParse &pp)
    {
        // Initial data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> centerA, centerB;
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
  
        // find rotation angle
	    pp.load("input_time", input_time);
	    pp.load("omega_BH", omega_BH);

    	// allow for rotation of input binary
        bh1_params.center[0] = centerA[0] + offsetA[0] * cos(omega_BH * input_time) - offsetA[1] * sin(omega_BH * input_time);
        bh1_params.center[1] = centerA[1] + offsetA[1] * cos(omega_BH * input_time) + offsetA[0] * sin(omega_BH * input_time);
        bh1_params.center[2] = centerA[2] + offsetA[2];
        bh2_params.center[0] = centerB[0] + offsetB[0] * cos(omega_BH * input_time) - offsetB[1] * sin(omega_BH * input_time);
        bh2_params.center[1] = centerB[1] + offsetB[1] * cos(omega_BH * input_time) + offsetB[0] * sin(omega_BH * input_time);
        bh2_params.center[2] = centerB[2] + offsetB[2];

    }
#endif /* USE_TWOPUNCTURES */

    /// Read shared parameters
    void read_shared_params(GRParmParse &pp)
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

	// ** Bounded Densities volume integration parameters **
	pp.load("inner_r", integration_params.min_integration_radius, 0.0);
        pp.load("outer_r", integration_params.max_integration_radius, L/2.0);

	// ** Params for second surface integration **
        pp.get("num_extraction_radii_2", extraction_params_2.num_extraction_radii);
        // -- make integration radius array
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
    bool track_punctures, calculate_constraint_norms, calculate_constraints;
    bool activate_Weyl_extraction, activate_flux_extraction, activate_integral;
    int puncture_tracking_level;
    std::vector<std::array<double, CH_SPACEDIM>> initial_puncture_coords;
    double G_Newton;
    double omega_BH, input_time;
    ScalarRotatingCloud::params_t initial_params;
    ComplexScalarPotential::params_t potential_params;
    std::string integral_filename;
    integration_params_t integration_params;
    SphericalExtraction::params_t extraction_params_2;
    
    // Collection of parameters necessary for initial conditions
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;

#ifdef USE_TWOPUNCTURES
    double tp_offset_plus, tp_offset_minus;
    TP::Parameters tp_params;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */

