#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

// How rho is filled at t=0. Kept apart from the bias below: the bias is what the
// dissipator holds the leads at for all time, this is only the starting fill.
struct InitialCondition {
    bool from_bias{true};   // Fill each lead at its own bias, leave the molecule empty
    double mu{0.0};         // Common starting level (eV) when from_bias is false
    double T{0.0};          // Starting temperature (K) when from_bias is false
    bool connected{false};  // Fill the eigenstates of the whole H instead of block by block
};

struct TBModelParameters {
    // Molecule grid dimensions
    int N_Mx{6};  // Molecule grid width  (N_My=1 -> 1D chain)
    int N_My{1};  // Molecule grid height
    double alpha_M{0.0};
    double beta_M{-0.2};

    // Left Lead (L) & Buffer (ML)
    bool has_L{true};
    int N_L{300};
    int N_ML{50};
    double alpha_L{0.0};
    double beta_L{-0.2};
    double beta_LM{-0.2};

    // Right Lead (R) & Buffer (MR)
    bool has_R{true};
    int N_R{300};
    int N_MR{50};
    double alpha_R{0.0};
    double beta_R{-0.2};
    double beta_MR{-0.2};

    // Upper Lead (U) & Buffer (MU)
    bool has_U{false};
    int N_U{300};
    int N_MU{50};
    double alpha_U{0.0};
    double beta_U{-0.2};
    double beta_UM{-0.2};

    // Down Lead (D) & Buffer (MD)
    bool has_D{false};
    int N_D{300};
    int N_MD{50};
    double alpha_D{0.0};
    double beta_D{-0.2};
    double beta_DM{-0.2};

    // Chemical potentials / Bias (eV), held by the dissipator for all time
    double bias_L{0.15};
    double bias_R{-0.15};
    double bias_U{0.0};
    double bias_D{0.0};

    // Fill at t=0, independent of the bias above
    InitialCondition init;
};

// Builds the Hamiltonian matrix in the site representation.
Eigen::MatrixXd build_hamiltonian(const TBModelParameters& params);

struct LeadInfo {
    std::string name; // "L", "R", "U", "D"
    int start_idx;
    int size;
    double mu; // Chemical potential for this lead
};

struct SystemPartition {
    int EM_start;
    int EM_size;
    std::vector<LeadInfo> leads;
};

struct StateRepresentation {
    Eigen::MatrixXd U; // Global transformation matrix
    Eigen::VectorXd eigenvalues; // Eigenvalues of isolated blocks
    Eigen::MatrixXd H_tilde; // Hamiltonian in state representation
    SystemPartition partition;
};

// Computes the state representation (Eq 5, 6, 7 of the paper)
StateRepresentation compute_state_representation(const Eigen::MatrixXd& H, const TBModelParameters& params);

