#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

struct TBModelParameters {
    // Left / Right leads (default: enabled)
    bool has_L{true};
    bool has_R{true};

    int N_L;
    int N_R;
    int N_ML;
    int N_MR;
    int N_Mx;  // Molecule grid width  (N_My=1 → 1D chain = paper default)
    int N_My;  // Molecule grid height

    // Optional Up / Down branches
    bool has_U{false};
    int N_U;
    int N_MU;
    double alpha_U;
    double beta_U;
    double beta_UM;

    bool has_D{false};
    int N_D;
    int N_MD;
    double alpha_D;
    double beta_D;
    double beta_DM;

    double alpha_L;
    double alpha_M;
    double alpha_R;

    double beta_L;
    double beta_M;
    double beta_R;
    double beta_LM;
    double beta_MR;

    // Chemical potentials
    double bias_L{0.15};
    double bias_R{-0.15};
    double bias_U{0.0};
    double bias_D{0.0};
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

