#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "hamiltonian.h"

struct LeadWorkspace {
    std::string name;
    int start;
    int size;
    Eigen::MatrixXcd V_lead_EM;
    Eigen::MatrixXcd V_EM_lead;
    Eigen::MatrixXcd rho0;
};

// Pre-allocated workspace for RK4 to avoid dynamic memory allocation in the inner loop
struct RK4Workspace {
    Eigen::MatrixXcd k1;
    Eigen::MatrixXcd k2;
    Eigen::MatrixXcd k3;
    Eigen::MatrixXcd k4;
    Eigen::MatrixXcd tmp_rho;
    
    Eigen::VectorXcd H_diag;
    Eigen::MatrixXcd omega; // Precomputed omega(i,j) = H_diag(i) - H_diag(j)
    int EM_start;
    int EM_size;
    std::vector<LeadWorkspace> leads;
    
    void setup(const StateRepresentation& rep, double T);
};

Eigen::MatrixXd compute_target_density_matrix_block(
    const Eigen::VectorXd& eigenvalues, 
    int start_idx, int size, 
    double mu, double T);

// In-place evaluation of the RHS to avoid allocations
void dlvn_rhs_state_rep_inplace(
    Eigen::MatrixXcd& drho_dt,
    const Eigen::MatrixXcd& rho_tilde,
    double gamma,
    RK4Workspace& ws);

void rk4_step_state_rep(
    Eigen::MatrixXcd& rho_tilde, 
    double gamma, double dt,
    RK4Workspace& ws);


