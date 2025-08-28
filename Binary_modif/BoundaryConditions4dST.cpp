#include "BoundaryConditions4dST.hpp"

// You may need these if used inside your implementation
#include "BoxIterator.H"
#include "Coordinates.hpp"
#include "UserVariables.hpp"

void BoundaryConditions4dST::fill_rhs_boundaries(const Side::LoHiSide a_side,
                                                 const GRLevelData &a_soln,
                                                 GRLevelData &a_rhs)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions4dST::fill_rhs_boundaries");

    FOR(idir)
    {
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            constexpr bool filling_rhs = true;
            custom_fill_boundary_cells_dir(a_side, a_soln, a_rhs, idir,
                                           boundary_condition,
                                           Interval(0, NUM_VARS - 1),
                                           VariableType::evolution, filling_rhs);
        }
    }
}

/// fill solution boundary conditions, e.g. after interpolation
void BoundaryConditions4dST::fill_solution_boundaries(const Side::LoHiSide a_side,
                                                  GRLevelData &a_state,
                                                  const Interval &a_comps)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_solution_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic and solution
        // boundary enforced in this direction
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            // same copying of cells which we require for the rhs solution
            // but tell it we are not filling the rhs for mixed condition
            if ((boundary_condition == REFLECTIVE_BC) ||
                (boundary_condition == EXTRAPOLATING_BC) ||
                (boundary_condition == MIXED_BC))
            {
                const bool filling_rhs = false;
                custom_fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                        boundary_condition, a_comps,
                                        VariableType::evolution, filling_rhs);
            }
        }
    }
}

/// fill diagnostic boundaries
void BoundaryConditions4dST::fill_diagnostic_boundaries(const Side::LoHiSide a_side,
                                                    GRLevelData &a_state,
                                                    const Interval &a_comps)
{
    CH_assert(is_defined);
    CH_TIME("BoundaryConditions::fill_diagnostic_boundaries");

    // cycle through the directions
    FOR(idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);
            // for any non reflective BC, we just want to fill the ghosts with
            // something non nan so set the boundary condition to be
            // EXTRAPOLATING
            if (boundary_condition != REFLECTIVE_BC)
            {
                boundary_condition = EXTRAPOLATING_BC;
            }
            const bool filling_rhs = false;
            custom_fill_boundary_cells_dir(a_side, a_state, a_state, idir,
                                    boundary_condition, a_comps,
                                    VariableType::diagnostic, filling_rhs);
        }
    }
}

/// Fill the boundary values appropriately based on the params set
/// in the direction dir
void BoundaryConditions4dST::custom_fill_boundary_cells_dir(
    const Side::LoHiSide a_side, const GRLevelData &a_soln, GRLevelData &a_out,
    const int dir, const int boundary_condition, const Interval &a_comps,
    const VariableType var_type, const bool filling_rhs)
{
    std::vector<int> comps_vector, sommerfeld_comps_vector,
        extrapolating_comps_vector;
    if (boundary_condition != MIXED_BC)
    {
        comps_vector.resize(a_comps.size());
        std::iota(comps_vector.begin(), comps_vector.end(), a_comps.begin());
    }
    else
    {
        for (int icomp = a_comps.begin(); icomp <= a_comps.end(); ++icomp)
        {
            if (m_params.mixed_bc_vars_map[icomp] == SOMMERFELD_BC)
                sommerfeld_comps_vector.push_back(icomp);
            else if (m_params.mixed_bc_vars_map[icomp] == EXTRAPOLATING_BC)
                extrapolating_comps_vector.push_back(icomp);
        }
    }

    // iterate through the boxes, shared amongst threads
    DataIterator dit = a_out.dataIterator();
    int nbox = dit.size();
#pragma omp parallel for default(shared)
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex dind = dit[ibox];
        FArrayBox &out_box = a_out[dind];
        const FArrayBox &soln_box = a_soln[dind];
        Box this_box = out_box.box();
        IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
        IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        this_box &= m_domain_box;
        // get the boundary box (may be Empty)
        Box boundary_box =
            get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

        // now we have the appropriate box, fill it!
        BoxIterator bit(boundary_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            switch (boundary_condition)
            {
            // simplest case - boundary values are set to zero
            case STATIC_BC:
            {
                for (int icomp = a_comps.begin(); icomp <= a_comps.end();
                     ++icomp)
                {
                    out_box(iv, icomp) = 0.0;
                }
                break;
            }
            // Sommerfeld is outgoing radiation - only applies to rhs
            case SOMMERFELD_BC:
            {
                custom_fill_sommerfeld_cell(out_box, soln_box, iv, comps_vector);
                break;
            }
            // Enforce a reflective symmetry in some direction
            case REFLECTIVE_BC:
            {
                fill_reflective_cell(out_box, iv, a_side, dir, comps_vector,
                                     var_type);
                break;
            }
            case EXTRAPOLATING_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir, comps_vector,
                                        m_params.extrapolation_order);
                break;
            }
            case MIXED_BC:
            {
                fill_extrapolating_cell(out_box, iv, a_side, dir,
                                        extrapolating_comps_vector,
                                        m_params.extrapolation_order);
                if (filling_rhs)
                {
                    custom_fill_sommerfeld_cell(out_box, soln_box, iv,
                                         sommerfeld_comps_vector);
                }
                break;
            }
            default:
                MayDay::Error(
                    "BoundaryCondition::Supplied boundary not supported.");
            } // end switch
        }     // end iterate over box
    }         // end iterate over boxes
}

void BoundaryConditions4dST::custom_fill_sommerfeld_cell(
    FArrayBox &rhs_box, const FArrayBox &soln_box, const IntVect iv,
    const std::vector<int> &sommerfeld_comps) const
{
    // assumes an asymptotic value + radial waves and permits them
    // to exit grid with minimal reflections
    // get real position on the grid
    RealVect loc(iv + 0.5 * RealVect::Unit);
    loc *= m_dx;
    loc -= m_center;
    double radius_squared = 0.0;
    FOR(i) { radius_squared += loc[i] * loc[i]; }
    double radius = sqrt(radius_squared);
    IntVect lo_local_offset = iv - soln_box.smallEnd();
    IntVect hi_local_offset = soln_box.bigEnd() - iv;

    // Apply Sommerfeld BCs to each variable in sommerfeld_comps
    for (int icomp : sommerfeld_comps)
    {
        rhs_box(iv, icomp) = 0.0;
        FOR(idir2)
        {
            IntVect iv_offset1 = iv;
            IntVect iv_offset2 = iv;
            double d1;
            // bit of work to get the right stencils for near
            // the edges of the domain, only using second order
            // stencils for now
            if (lo_local_offset[idir2] < 1)
            {
                // near lo end
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += +2;
                d1 = 1.0 / m_dx *
                     (-1.5 * soln_box(iv, icomp) +
                      2.0 * soln_box(iv_offset1, icomp) -
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else if (hi_local_offset[idir2] < 1)
            {
                // near hi end
                iv_offset1[idir2] += -1;
                iv_offset2[idir2] += -2;
                d1 = 1.0 / m_dx *
                     (+1.5 * soln_box(iv, icomp) -
                      2.0 * soln_box(iv_offset1, icomp) +
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else
            {
                // normal case
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += -1;
                d1 =
                    0.5 / m_dx *
                    (soln_box(iv_offset1, icomp) - soln_box(iv_offset2, icomp));
            }

            // for each direction add dphidx * x^i / r
            rhs_box(iv, icomp) += -d1 * loc[idir2] / radius;
        }

        // asymptotic values - these need to have been set in
        // the params file
        rhs_box(iv, icomp) +=
            (m_params.vars_asymptotic_values[icomp] - soln_box(iv, icomp)) /
            radius;
	if (icomp == c_phi || icomp == c_Pi)
        {
            rhs_box(iv, icomp) += -m_phys_params.scalar_mass * soln_box(iv, icomp);
        }
    }
}
