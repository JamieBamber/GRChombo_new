/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BoxLoops.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "NewMatterConstraints.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

#include "ComplexScalarField.hpp"
#include "ComplexScalarPotential.hpp"
#include "MatterOnly.hpp"
#include "DensityAndMomAndSourceTerm.hpp"
#include "ScalarRotatingCloud.hpp"

#include "BoundedDensities.hpp"
#include "ForceExtraction.hpp"
#include "Ylm_Extraction.hpp"

// Things to do during the advance step after RK4 steps
void BinaryBHLevel::specificAdvance()
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::specificAdvance " << m_level << endl;

    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// This initial data uses an approximation for the metric which
// is valid for small boosts
void BinaryBHLevel::initialData()
{
    CH_TIME("BinaryBHLevel::initialData");
    if (m_verbosity)
        pout() << "BinaryBHLevel::initialData " << m_level << endl;
#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    BoxLoops::loop(two_punctures_initial_data, m_state_new, m_state_new,
                   INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);
    
    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(make_compute_pack(SetValue(0.), binary), m_state_new,
                  m_state_new, INCLUDE_GHOST_CELLS);
#endif
    // scalar field compute class
    ScalarRotatingCloud initial_sf(m_p.initial_params, m_dx);
    BoxLoops::loop(initial_sf, m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
}

// Calculate RHS during RK4 substeps
void BinaryBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                    const double a_time)
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::specificEvalRHS " << m_level << endl;
 
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    ComplexScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);

        // Calculate CCZ4 right hand side
    if (m_p.max_spatial_derivative_order == 4)
    {
	MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge, FourthOrderDerivatives> my_ccz4_matter(
           scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
           m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge, SixthOrderDerivatives> my_ccz4_matter(
           scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
           m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void BinaryBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                      const GRLevelData &a_rhs, Real a_dt)
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::specificUpdateODE " << m_level << endl;

    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BinaryBHLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// specify the cells to tag
void BinaryBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                            const FArrayBox &current_state)
{
    if (m_p.track_punctures)
    {
        std::vector<double> puncture_masses;
#ifdef USE_TWOPUNCTURES
        // use calculated bare masses from TwoPunctures
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
#else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
#endif /* USE_TWOPUNCTURES */
        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
        BoxLoops::loop(ChiPunctureExtractionTaggingCriterion(
                           m_dx, m_level, m_p.max_level, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses),
                       current_state, tagging_criterion);
    }
    else
    {
        BoxLoops::loop(ChiExtractionTaggingCriterion(m_dx, m_level,
                                                     m_p.extraction_params,
                                                     m_p.activate_extraction),
                       current_state, tagging_criterion);
    }
}

void BinaryBHLevel::specificPostTimeStep()
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::specificPostTimeStep " << m_level << endl;
    CH_TIME("BinaryBHLevel::specificPostTimeStep");

        bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
    	if (m_verbosity)
		    pout() << "calculating Weyl" << endl;
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(Weyl4(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraints)
    {
        fillAllGhosts();
	    if (m_verbosity)
		    pout() << "calculating constraints" << endl;
	    // At any level, but after the coarsest timestep
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_constraints = at_level_timestep_multiple(min_level);
        if (calculate_constraints){
            ComplexScalarPotential potential(m_p.potential_params);
            ScalarFieldWithPotential scalar_field(potential);
	    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)), m_state_new, m_state_diagnostics,
        	               EXCLUDE_GHOST_CELLS);
		if (m_level == 0 && m_p.calculate_constraint_norms)
        	{
            	AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            	double L2_Ham = amr_reductions.norm(c_Ham);
            	double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            	SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                         	m_restart_time, SmallDataIO::APPEND,
                                         	first_step);
            	constraints_file.remove_duplicate_time_data();
            	if (first_step)
            	{
                	constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            	}
            	constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        	}
	}
    }
    
    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
	    if (m_verbosity)
		    pout() << "tracking punctures" << endl;
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time, m_dt,
                                                   write_punctures);
    }

    if (false) //(m_p.activate_flux_extraction == 1)
    {
	// At any level, but after the coarsest timestep
	if (m_verbosity)
                    pout() << "starting flux extraction section" << endl;
        int min_level = 0; //m_p.extraction_params_2.min_extraction_level();
	bool calculate_densities = at_level_timestep_multiple(min_level);
	if (calculate_densities)
	{
		// Calculate energy and angular momentum fluxes and densities
		ComplexScalarPotential potential(m_p.potential_params);
		ScalarFieldWithPotential scalar_field(potential);
		if (m_verbosity)
        	pout() << "Now making densities and momenta" << endl;
		BoxLoops::loop(DensityAndMomAndSourceTerm<ScalarFieldWithPotential>(
	                	scalar_field, m_dx, m_p.center, 0, m_p.ccz4_params),
	                	m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
	}
	
	if (m_level == 0)
	{
    	// fill only specific ghosts
    	bool fill_ghosts = false;
    	m_gr_amr.m_interpolator->refresh(fill_ghosts);
    	m_gr_amr.fill_multilevel_ghosts(
        	VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
        	min_level);
		//
	 	if (m_verbosity)
		   	pout() << "now doing flux integrals" << endl;
		ForceExtraction my_flux_extraction(m_p.extraction_params_2, m_dt, m_time,
	                              	"Force_integrals", m_restart_time);
		my_flux_extraction.execute_query(m_gr_amr.m_interpolator);
		if (m_verbosity)
               	pout() << "now doing phi integrals" << endl;
		PhiExtraction my_phi_extraction(m_p.extraction_params_2, "", "", m_dt, m_time, m_restart_time);
		my_phi_extraction.execute_query(m_gr_amr.m_interpolator);		
	}
    }

    if (false) //(m_p.activate_integral)
    {
        if (m_verbosity)
                    pout() << "starting integral section" << endl;
	int min_level = 0; 
	bool calculate_densities = at_level_timestep_multiple(min_level);
        if (calculate_densities)
	    {
	        // Now, over all levels, we instead exclude an inner sphere around the black holes
	        // with radius set by r_inner
	        BoundedDensities density_bound(m_p.integration_params, m_dx, m_p.center);
	        BoxLoops::loop(
	                       density_bound,
	                       m_state_diagnostics, m_state_diagnostics, INCLUDE_GHOST_CELLS, disable_simd());
	    }

	    if (m_level == 0)
	    {
	        // Now do the integral
		    if (m_verbosity)
			    pout() << "now integrating over the volume" << endl;
	        AMRReductions<VariableType::diagnostic> amr_reductions_2(m_gr_amr);
	        double rho_sum = amr_reductions_2.sum(c_rho);
	        double rhoJ_sum = amr_reductions_2.sum(c_rho_azimuth);
	        double C_t_sum = amr_reductions_2.sum(c_C_t);
	        double C_phi_sum = amr_reductions_2.sum(c_C_phi);
	        
	       	SmallDataIO integral_file_bh(m_p.integral_filename, m_dt, m_time,
	                                  m_restart_time, SmallDataIO::APPEND,
	                                  first_step);
	        // remove any duplicate data if this is post restart
	        integral_file_bh.remove_duplicate_time_data();
	        std::vector<double> data_for_writing_bh = {rho_sum, rhoJ_sum, C_t_sum, C_phi_sum};
	        // write data
	        if (first_step)
	        {
	          integral_file_bh.write_header_line({"# r_min = "+std::to_string(m_p.integration_params.min_integration_radius), "r_max = "+std::to_string(m_p.integration_params.max_integration_radius)});
	          integral_file_bh.write_header_line({"rho", "rhoJ", "C_t", "C_phi"});
	        }
		    integral_file_bh.write_time_data_line(data_for_writing_bh);    
	    }
     }
}

#ifdef CH_USE_HDF5
// Things to do before a plot level
void BinaryBHLevel::prePlotLevel()
{
    if (m_verbosity)
        pout() << "BinaryBHLevel::prePlotLevel " << m_level << endl;   

    ComplexScalarPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    fillAllGhosts();
    if (m_p.calculate_constraints) {
	    BoxLoops::loop(MatterConstraints<ScalarFieldWithPotential>(scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom1, c_Mom3)), m_state_new, m_state_diagnostics,
        	               EXCLUDE_GHOST_CELLS);
    }
    pout() << "Now making densities and momenta" << endl;
    BoxLoops::loop(DensityAndMomAndSourceTerm<ScalarFieldWithPotential>(
                       scalar_field, m_dx, m_p.center, 0, m_p.ccz4_params),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    // Populate the Weyl Scalar values on the grid                                                                                                                              
    if (m_p.activate_extraction == 1)
        BoxLoops::loop(Weyl4(m_p.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
}
#endif /* CH_USE_HDF5 */
