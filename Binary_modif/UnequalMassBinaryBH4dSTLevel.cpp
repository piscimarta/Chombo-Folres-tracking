/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "UnequalMassBinaryBH4dSTLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BinaryPunctureTaggingCriterion.hpp"
#include "BoxLoops.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "ComputePack.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "ModifiedGravityWeyl4.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "PunctureTracker.hpp"
#include "RHSDiagnostics.hpp"
#include "RhoDiagnostics.hpp"
#include "ScalarExtraction.hpp"
#include "SetValue.hpp"
#include "BoundaryConditions4dST.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "WeylExtraction.hpp"

// Things to do during the advance step after RK4 steps
void UnequalMassBinaryBH4dSTLevel::specificAdvance()
{
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
void UnequalMassBinaryBH4dSTLevel::initialData()
{
    CH_TIME("UnequalMassBinaryBH4dSTLevel::initialData");
    if (m_verbosity)
        pout() << "UnequalMassBinaryBH4dSTLevel::initialData " << m_level
               << endl;
#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(
        m_dx, m_p.center, m_tp_amr.m_two_punctures);
    // Can't use simd with this initial data
    // BoxLoops::loop(two_punctures_initial_data, m_state_new, m_state_new,
    //               INCLUDE_GHOST_CELLS, disable_simd());
    BoxLoops::loop(
        make_compute_pack(two_punctures_initial_data,
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
#else
    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), binary,
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);
#endif
}

// Calculate RHS during RK4 substeps
void UnequalMassBinaryBH4dSTLevel::specificEvalRHS(GRLevelData &a_soln,
                                                   GRLevelData &a_rhs,
                                                   const double a_time)
{
    // Enforce positive chi and alpha and trace free A
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate ModifiedCCZ4 right hand side with theory_t =
    // FourDerivScalarTensor
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    if (m_p.max_spatial_derivative_order == 4)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, FourthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma, m_p.rescale_sigma,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, SixthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma, m_p.rescale_sigma,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// enforce trace removal during RK4 substeps
void UnequalMassBinaryBH4dSTLevel::specificUpdateODE(GRLevelData &a_soln,
                                                     const GRLevelData &a_rhs,
                                                     Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void UnequalMassBinaryBH4dSTLevel::preTagCells()
{
    // We only use chi in the tagging criterion so only fill the ghosts for chi
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
}

// specify the cells to tag
void UnequalMassBinaryBH4dSTLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state)
{
    if (m_p.track_punctures)
    {
        std::vector<double> puncture_masses;
        // std::array<double,2> puncture_masses;
#ifdef USE_TWOPUNCTURES
        // use calculated bare masses from TwoPunctures
        puncture_masses = {m_tp_amr.m_two_punctures.mm,
                           m_tp_amr.m_two_punctures.mp};
#else
        puncture_masses = {m_p.bh1_params.mass, m_p.bh2_params.mass};
#endif /* USE_TWOPUNCTURES */
        auto puncture_coords =
            m_bh_amr.m_puncture_tracker.get_puncture_coords();
        // different tagging criterion for unequal mass binaries
        BoxLoops::loop(BinaryPunctureTaggingCriterion<FourthOrderDerivatives>(
                           m_dx, m_level, m_p.tag_horizons_max_levels,
                           m_p.tag_punctures_max_levels, m_p.extraction_params,
                           puncture_coords, m_p.activate_extraction,
                           m_p.track_punctures, puncture_masses,
                           m_p.bh_tagging_buffers,
                           m_p.puncture_tag_min_separation),
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

void UnequalMassBinaryBH4dSTLevel::specificPostTimeStep()
{
    CH_TIME("UnequalMassBinaryBH4dSTLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1 || m_p.activate_scalar_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            CouplingAndPotential coupling_and_potential(
                m_p.coupling_and_potential_params);
            FourDerivScalarTensorWithCouplingAndPotential fdst(
                coupling_and_potential, m_p.G_Newton);
            ModifiedPunctureGauge modified_puncture_gauge(
                m_p.modified_ccz4_params);
            ModifiedGravityWeyl4<FourDerivScalarTensorWithCouplingAndPotential,
                                 ModifiedPunctureGauge, FourthOrderDerivatives>
                weyl4(fdst, m_p.modified_ccz4_params, modified_puncture_gauge,
                      m_p.extraction_params.extraction_center, m_dx, m_p.sigma,
                      CCZ4RHS<>::USE_CCZ4, m_p.rescale_sigma);
            // CCZ4 is required since this code only works in this
            // formulation
            BoxLoops::loop(weyl4, m_state_new, m_state_diagnostics,
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
                if (m_p.activate_extraction)
                {
                    WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                                 m_time, first_step,
                                                 m_restart_time);
                    my_extraction.execute_query(m_gr_amr.m_interpolator);
                }

                if (m_p.activate_scalar_extraction)
                {
                    ScalarExtraction phi_extraction(
                        m_p.scalar_extraction_params, m_dt, m_time, first_step,
                        m_restart_time);
                    phi_extraction.execute_query(m_gr_amr.m_interpolator);
                }
            }
        }
    }

    fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);
    RHSDiagnostics<FourDerivScalarTensorWithCouplingAndPotential,
                   ModifiedPunctureGauge, FourthOrderDerivatives>
        rhs_diagnostics(fdst, m_p.modified_ccz4_params, modified_puncture_gauge,
                        m_dx, m_p.sigma, m_p.rescale_sigma, m_p.center, m_p.G_Newton);
    auto compute_pack = make_compute_pack(rho_diagnostics, rhs_diagnostics);
    BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);

    if (m_p.calculate_constraint_norms)
    {
        CouplingAndPotential coupling_and_potential(
            m_p.coupling_and_potential_params);
        FourDerivScalarTensorWithCouplingAndPotential fdst(
            coupling_and_potential, m_p.G_Newton);
        fillAllGhosts();
        BoxLoops::loop(ModifiedGravityConstraints<
                           FourDerivScalarTensorWithCouplingAndPotential>(
                           fdst, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                           Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            double L2_Ham = amr_reductions.norm(c_Ham);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
            SmallDataIO constraints_file(m_p.data_path + "constraint_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            constraints_file.remove_duplicate_time_data();
            if (first_step)
            {
                constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
            }
            constraints_file.write_time_data_line({L2_Ham, L2_Mom});
        }
    }

    // do puncture tracking on requested level
    if (m_p.track_punctures && m_level == m_p.puncture_tracking_level)
    {
        CH_TIME("PunctureTracking");
        // only do the write out for every coarsest level timestep
        int coarsest_level = 0;
        bool write_punctures = at_level_timestep_multiple(coarsest_level);
        m_bh_amr.m_puncture_tracker.execute_tracking(m_time, m_restart_time,
                                                     m_dt, write_punctures);
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        if (m_p.AH_set_origins_to_punctures && m_p.track_punctures)
        {
            m_bh_amr.m_ah_finder.set_origins(
                m_bh_amr.m_puncture_tracker.get_puncture_coords());
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif
}

#ifdef CH_USE_HDF5
// Things to do before a plot level - need to calculate the Weyl scalars
void UnequalMassBinaryBH4dSTLevel::prePlotLevel()
{
    fillAllGhosts();
    if (m_p.activate_extraction == 1)
    {
        CouplingAndPotential coupling_and_potential(
            m_p.coupling_and_potential_params);
        FourDerivScalarTensorWithCouplingAndPotential fdst(
            coupling_and_potential, m_p.G_Newton);
        ModifiedGravityConstraints<
            FourDerivScalarTensorWithCouplingAndPotential>
            constraints(fdst, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                        Interval(c_Mom1, c_Mom3));
        ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
        ModifiedGravityWeyl4<FourDerivScalarTensorWithCouplingAndPotential,
                             ModifiedPunctureGauge, FourthOrderDerivatives>
            weyl4(fdst, m_p.modified_ccz4_params, modified_puncture_gauge,
                  m_p.extraction_params.extraction_center, m_dx, m_p.sigma,
                  CCZ4RHS<>::USE_CCZ4, m_p.rescale_sigma);
        // CCZ4 is required since this code only works in this formulation
        auto compute_pack = make_compute_pack(weyl4, constraints);
        BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
    }
}
#endif /* CH_USE_HDF5 */
