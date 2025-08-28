#ifndef BOUNDARYCONDITIONS4DST_HPP_
#define BOUNDARYCONDITIONS4DST_HPP_

#include "BoundaryConditions.hpp"
#include "CouplingAndPotential.hpp"

class BoundaryConditions4dST : public BoundaryConditions
{
  protected:
    CouplingAndPotential::params_t m_phys_params;

  public:
    BoundaryConditions4dST() = default; // ✅ Only default constructor

    void define(const BoundaryConditions::params_t &bc_params,
                const CouplingAndPotential::params_t &phys_params,
                double dx,
                std::array<double, CH_SPACEDIM> center,
                ProblemDomain domain,
                int num_ghosts)
    {
        m_phys_params = phys_params;
        BoundaryConditions::define(dx, center, bc_params, domain, num_ghosts); // ✅ base call
    }

    virtual void fill_solution_boundaries(const Side::LoHiSide side,
                                      GRLevelData &state,
                                      const Interval &comps);

    virtual void fill_rhs_boundaries(const Side::LoHiSide side,
                                 const GRLevelData &soln,
                                 GRLevelData &rhs);

    virtual void fill_diagnostic_boundaries(const Side::LoHiSide side,
                                        GRLevelData &state,
                                        const Interval &comps);

    // 🔽 New declarations — required if defined in the .cpp
    void custom_fill_boundary_cells_dir(const Side::LoHiSide a_side,
                                        const GRLevelData &a_soln,
                                        GRLevelData &a_out,
                                        const int dir,
                                        const int boundary_condition,
                                        const Interval &a_comps,
                                        const VariableType var_type,
                                        const bool filling_rhs);

    void custom_fill_sommerfeld_cell(FArrayBox &rhs_box,
                                     const FArrayBox &soln_box,
                                     const IntVect iv,
                                     const std::vector<int> &sommerfeld_comps) const;
};

#endif /* BOUNDARYCONDITIONS4DST_HPP_ */

