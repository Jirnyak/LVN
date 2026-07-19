#include "lvn_dynamics.h"
#include <cmath>
#include <iostream>

Eigen::MatrixXd compute_target_density_matrix_block(
    const Eigen::VectorXd& eigenvalues, 
    int start_idx, int size, 
    double mu, double T) 
{
    Eigen::MatrixXd rho0 = Eigen::MatrixXd::Zero(size, size);
    const double k_B = 8.617333262145e-5; // Boltzmann constant in eV/K
    
    for (int i = 0; i < size; ++i) {
        double E = eigenvalues(start_idx + i);
        double f = 0.0;
        if (T < 1e-10) {
            f = (E <= mu) ? 1.0 : 0.0;
        } else {
            f = 1.0 / (std::exp((E - mu) / (k_B * T)) + 1.0);
        }
        rho0(i, i) = f;
    }
    return rho0;
}

void RK4Workspace::setup(const StateRepresentation& rep, double T) {
    int N = rep.H_tilde.rows();
    k1.resize(N, N);
    k2.resize(N, N);
    k3.resize(N, N);
    k4.resize(N, N);
    tmp_rho.resize(N, N);

    Eigen::MatrixXcd H_c = rep.H_tilde.cast<std::complex<double>>();
    H_diag = H_c.diagonal();
    // Precompute energy difference matrix: omega(i,j) = E_i - E_j
    omega.resize(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            omega(i, j) = H_diag(i) - H_diag(j);
    EM_start = rep.partition.EM_start;
    EM_size = rep.partition.EM_size;

    leads.clear();
    for (size_t i = 0; i < rep.partition.leads.size(); ++i) {
        const auto& li = rep.partition.leads[i];
        if (li.size > 0) {
            LeadWorkspace lw;
            lw.name = li.name;
            lw.start = li.start_idx;
            lw.size = li.size;
            lw.V_lead_EM = H_c.block(li.start_idx, EM_start, li.size, EM_size);
            lw.V_EM_lead = H_c.block(EM_start, li.start_idx, EM_size, li.size);
            lw.rho0 = compute_target_density_matrix_block(rep.eigenvalues, li.start_idx, li.size, li.mu, T).cast<std::complex<double>>();
            leads.push_back(lw);
        }
    }
}

void dlvn_rhs_state_rep_inplace(
    Eigen::MatrixXcd& drho_dt,
    const Eigen::MatrixXcd& rho_tilde,
    double gamma,
    RK4Workspace& ws) 
{
    const double hbar_eV_fs = 0.6582119569;
    const int N = rho_tilde.rows();

    // 1. Commutator [H, rho] = H * rho - rho * H
    // 1a. Diagonal part of H (precomputed energy differences)
    drho_dt.noalias() = ws.omega.cwiseProduct(rho_tilde);

    // 1b. Off-diagonal couplings between each Lead and EM
    for (size_t l = 0; l < ws.leads.size(); ++l) {
        const auto& lead = ws.leads[l];
        // H_offdiag * rho
        drho_dt.block(lead.start, 0, lead.size, N).noalias() += lead.V_lead_EM * rho_tilde.block(ws.EM_start, 0, ws.EM_size, N);
        drho_dt.block(ws.EM_start, 0, ws.EM_size, N).noalias() += lead.V_EM_lead * rho_tilde.block(lead.start, 0, lead.size, N);

        // -rho * H_offdiag
        drho_dt.block(0, lead.start, N, lead.size).noalias() -= rho_tilde.block(0, ws.EM_start, N, ws.EM_size) * lead.V_EM_lead;
        drho_dt.block(0, ws.EM_start, N, ws.EM_size).noalias() -= rho_tilde.block(0, lead.start, N, lead.size) * lead.V_lead_EM;
    }

    drho_dt *= std::complex<double>(0.0, -1.0 / hbar_eV_fs);

    // 2. Dissipator / Driving term D[rho]
    for (size_t l = 0; l < ws.leads.size(); ++l) {
        const auto& lead = ws.leads[l];
        // Self-relaxation on diagonal block (alpha, alpha)
        drho_dt.block(lead.start, lead.start, lead.size, lead.size) -= gamma * (rho_tilde.block(lead.start, lead.start, lead.size, lead.size) - lead.rho0);

        // Dephasing with Extended Molecule (alpha, EM) and (EM, alpha)
        drho_dt.block(lead.start, ws.EM_start, lead.size, ws.EM_size) -= 0.5 * gamma * rho_tilde.block(lead.start, ws.EM_start, lead.size, ws.EM_size);
        drho_dt.block(ws.EM_start, lead.start, ws.EM_size, lead.size) -= 0.5 * gamma * rho_tilde.block(ws.EM_start, lead.start, ws.EM_size, lead.size);

        // Dephasing with other leads (alpha, beta) for beta != alpha
        for (size_t m = 0; m < ws.leads.size(); ++m) {
            const auto& other = ws.leads[m];
            if (other.start != lead.start) {
                drho_dt.block(lead.start, other.start, lead.size, other.size) -= gamma * rho_tilde.block(lead.start, other.start, lead.size, other.size);
            }
        }
    }
}

void rk4_step_state_rep(
    Eigen::MatrixXcd& rho_tilde, 
    double gamma, double dt,
    RK4Workspace& ws) 
{
    dlvn_rhs_state_rep_inplace(ws.k1, rho_tilde, gamma, ws);
    ws.k1 *= dt;

    ws.tmp_rho.noalias() = rho_tilde + 0.5 * ws.k1;
    dlvn_rhs_state_rep_inplace(ws.k2, ws.tmp_rho, gamma, ws);
    ws.k2 *= dt;

    ws.tmp_rho.noalias() = rho_tilde + 0.5 * ws.k2;
    dlvn_rhs_state_rep_inplace(ws.k3, ws.tmp_rho, gamma, ws);
    ws.k3 *= dt;

    ws.tmp_rho.noalias() = rho_tilde + ws.k3;
    dlvn_rhs_state_rep_inplace(ws.k4, ws.tmp_rho, gamma, ws);
    ws.k4 *= dt;
    
    rho_tilde += (ws.k1 + 2.0 * ws.k2 + 2.0 * ws.k3 + ws.k4) / 6.0;
}

