#include "hamiltonian.h"
#include <cmath>
#include <iostream>

Eigen::MatrixXd build_hamiltonian(const TBModelParameters& p) {
    int N_L  = p.has_L ? p.N_L * p.N_My  : 0;
    int N_R  = p.has_R ? p.N_R * p.N_My  : 0;
    int N_ML = p.has_L ? p.N_ML * p.N_My : 0;
    int N_MR = p.has_R ? p.N_MR * p.N_My : 0;
    
    int N_MU = (p.has_U && p.N_U > 0) ? p.N_MU * p.N_Mx : 0;
    int N_MD = (p.has_D && p.N_D > 0) ? p.N_MD * p.N_Mx : 0;
    int N_U = (p.has_U && p.N_U > 0) ? p.N_U * p.N_Mx : 0;
    int N_D = (p.has_D && p.N_D > 0) ? p.N_D * p.N_Mx : 0;
    
    int N_M_total = p.N_Mx * p.N_My;
    int N_EM = N_ML + N_M_total + N_MR + N_MU + N_MD;
    int total_sites = N_L + N_EM + N_R + N_U + N_D;

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_sites, total_sites);

    int idx_L = 0;
    int idx_ML = idx_L + N_L;
    int idx_M = idx_ML + N_ML;
    int idx_MR = idx_M + N_M_total;
    int idx_MU = idx_MR + N_MR;
    int idx_MD = idx_MU + N_MU;
    int idx_R = idx_MD + N_MD;
    int idx_U = idx_R + N_R;
    int idx_D = idx_U + N_U;

    // Helper to build a 2D block with nearest-neighbor hopping
    // nx is the length of the first coordinate, ny is the length of the second coordinate.
    // Indexing: s = start_idx + i * ny + j
    auto build_2d_block = [&](int start_idx, int nx, int ny, double alpha, double beta) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int s = start_idx + i * ny + j;
                H(s, s) = alpha;
                if (i < nx - 1) {
                    int sr = start_idx + (i + 1) * ny + j;
                    H(s, sr) = beta;
                    H(sr, s) = beta;
                }
                if (j < ny - 1) {
                    int sd = start_idx + i * ny + (j + 1);
                    H(s, sd) = beta;
                    H(sd, s) = beta;
                }
            }
        }
    };

    auto l_idx = [&](int ix, int iy) { return idx_L + ix * p.N_My + iy; };
    auto ml_idx = [&](int ix, int iy) { return idx_ML + ix * p.N_My + iy; };
    auto m_idx = [&](int ix, int iy) { return idx_M + ix * p.N_My + iy; };
    auto mr_idx = [&](int ix, int iy) { return idx_MR + ix * p.N_My + iy; };
    auto r_idx = [&](int ix, int iy) { return idx_R + ix * p.N_My + iy; };
    
    // For U and D leads, they expand along X axis. 
    // Their first index is length (Y axis), second index is width (X axis).
    auto mu_idx = [&](int iy, int ix) { return idx_MU + iy * p.N_Mx + ix; };
    auto u_idx = [&](int iy, int ix) { return idx_U + iy * p.N_Mx + ix; };
    auto md_idx = [&](int iy, int ix) { return idx_MD + iy * p.N_Mx + ix; };
    auto d_idx = [&](int iy, int ix) { return idx_D + iy * p.N_Mx + ix; };

    // 1. Left Lead (if enabled)
    if (p.has_L && p.N_L > 0) {
        build_2d_block(idx_L, p.N_L, p.N_My, p.alpha_L, p.beta_L);
        // Coupling L to ML (1-to-1)
        if (p.N_ML > 0) {
            for (int iy = 0; iy < p.N_My; ++iy) {
                int sl = l_idx(p.N_L - 1, iy);
                int sml = ml_idx(0, iy);
                H(sl, sml) = p.beta_L;
                H(sml, sl) = p.beta_L;
            }
        }
    }

    // 2. ML (Extended Molecule Left part)
    if (p.has_L && p.N_ML > 0) {
        build_2d_block(idx_ML, p.N_ML, p.N_My, p.alpha_L, p.beta_L);
        // Coupling ML[last] to M (1-to-1)
        if (p.N_Mx > 0) {
            for (int iy = 0; iy < p.N_My; ++iy) {
                int sml = ml_idx(p.N_ML - 1, iy);
                int sm = m_idx(0, iy);
                H(sml, sm) = p.beta_LM;
                H(sm, sml) = p.beta_LM;
            }
        }
    }

    // 3. Molecule
    if (p.N_Mx > 0 && p.N_My > 0) {
        build_2d_block(idx_M, p.N_Mx, p.N_My, p.alpha_M, p.beta_M);
    }

    // 4. MR (Extended Molecule Right part)
    if (p.has_R && p.N_MR > 0) {
        build_2d_block(idx_MR, p.N_MR, p.N_My, p.alpha_R, p.beta_R);
        // Coupling M to MR[0] (1-to-1)
        if (p.N_Mx > 0) {
            for (int iy = 0; iy < p.N_My; ++iy) {
                int sm = m_idx(p.N_Mx - 1, iy);
                int smr = mr_idx(0, iy);
                H(sm, smr) = p.beta_MR;
                H(smr, sm) = p.beta_MR;
            }
        }
        // Coupling MR to R (1-to-1)
        if (p.N_R > 0) {
            for (int iy = 0; iy < p.N_My; ++iy) {
                int smr = mr_idx(p.N_MR - 1, iy);
                int sr = r_idx(0, iy);
                H(smr, sr) = p.beta_R;
                H(sr, smr) = p.beta_R;
            }
        }
    }

    // 5. Right Lead
    if (p.has_R && p.N_R > 0) {
        build_2d_block(idx_R, p.N_R, p.N_My, p.alpha_R, p.beta_R);
    }

    // 6. Upper Branch (MU + U)
    if (p.has_U && p.N_U > 0) {
        if (p.N_MU > 0) {
            build_2d_block(idx_MU, p.N_MU, p.N_Mx, p.alpha_U, p.beta_U);
            // Coupling MU[0] to M top edge (iy = 0) (1-to-1)
            if (p.N_My > 0) {
                for (int ix = 0; ix < p.N_Mx; ++ix) {
                    int smu = mu_idx(0, ix);
                    int sm = m_idx(ix, 0);
                    H(smu, sm) = p.beta_UM;
                    H(sm, smu) = p.beta_UM;
                }
            }
            // Coupling MU[last] to U[0]
            for (int ix = 0; ix < p.N_Mx; ++ix) {
                int smu = mu_idx(p.N_MU - 1, ix);
                int su = u_idx(0, ix);
                H(smu, su) = p.beta_U;
                H(su, smu) = p.beta_U;
            }
        }
        build_2d_block(idx_U, p.N_U, p.N_Mx, p.alpha_U, p.beta_U);
    }

    // 7. Down Branch (MD + D)
    if (p.has_D && p.N_D > 0) {
        if (p.N_MD > 0) {
            build_2d_block(idx_MD, p.N_MD, p.N_Mx, p.alpha_D, p.beta_D);
            // Coupling MD[0] to M bottom edge (iy = N_My - 1) (1-to-1)
            if (p.N_My > 0) {
                for (int ix = 0; ix < p.N_Mx; ++ix) {
                    int smd = md_idx(0, ix);
                    int sm = m_idx(ix, p.N_My - 1);
                    H(smd, sm) = p.beta_DM;
                    H(sm, smd) = p.beta_DM;
                }
            }
            // Coupling MD[last] to D[0]
            for (int ix = 0; ix < p.N_Mx; ++ix) {
                int smd = md_idx(p.N_MD - 1, ix);
                int sd = d_idx(0, ix);
                H(smd, sd) = p.beta_D;
                H(sd, smd) = p.beta_D;
            }
        }
        build_2d_block(idx_D, p.N_D, p.N_Mx, p.alpha_D, p.beta_D);
    }
    
    return H;
}

StateRepresentation compute_state_representation(const Eigen::MatrixXd& H, const TBModelParameters& p) {
    int N_L  = p.has_L ? p.N_L * p.N_My  : 0;
    int N_R  = p.has_R ? p.N_R * p.N_My  : 0;
    int N_ML = p.has_L ? p.N_ML * p.N_My : 0;
    int N_MR = p.has_R ? p.N_MR * p.N_My : 0;
    
    int N_MU = (p.has_U && p.N_U > 0) ? p.N_MU * p.N_Mx : 0;
    int N_MD = (p.has_D && p.N_D > 0) ? p.N_MD * p.N_Mx : 0;
    int N_U = (p.has_U && p.N_U > 0) ? p.N_U * p.N_Mx : 0;
    int N_D = (p.has_D && p.N_D > 0) ? p.N_D * p.N_Mx : 0;
    
    int N_M_total = p.N_Mx * p.N_My;
    int N_EM = N_ML + N_M_total + N_MR + N_MU + N_MD;
    int N = H.rows();

    int idx_L = 0;
    int idx_EM = idx_L + N_L;
    int idx_R = idx_EM + N_EM;
    int idx_U = idx_R + N_R;
    int idx_D = idx_U + N_U;

    SystemPartition partition;
    partition.EM_start = idx_EM;
    partition.EM_size = N_EM;
    
    // Independent bias for each lead
    if (N_L > 0) partition.leads.push_back({"L", idx_L, N_L, p.bias_L});
    if (N_R > 0) partition.leads.push_back({"R", idx_R, N_R, p.bias_R});
    if (N_U > 0) partition.leads.push_back({"U", idx_U, N_U, p.bias_U});
    if (N_D > 0) partition.leads.push_back({"D", idx_D, N_D, p.bias_D});

    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd evals = Eigen::VectorXd::Zero(N);

    // Diagonalize blocks
    if (N_L > 0) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H.block(idx_L, idx_L, N_L, N_L));
        U.block(idx_L, idx_L, N_L, N_L) = solver.eigenvectors();
        evals.segment(idx_L, N_L) = solver.eigenvalues();
    }
    if (N_EM > 0) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H.block(idx_EM, idx_EM, N_EM, N_EM));
        U.block(idx_EM, idx_EM, N_EM, N_EM) = solver.eigenvectors();
        evals.segment(idx_EM, N_EM) = solver.eigenvalues();
    }
    if (N_R > 0) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H.block(idx_R, idx_R, N_R, N_R));
        U.block(idx_R, idx_R, N_R, N_R) = solver.eigenvectors();
        evals.segment(idx_R, N_R) = solver.eigenvalues();
    }
    if (N_U > 0) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H.block(idx_U, idx_U, N_U, N_U));
        U.block(idx_U, idx_U, N_U, N_U) = solver.eigenvectors();
        evals.segment(idx_U, N_U) = solver.eigenvalues();
    }
    if (N_D > 0) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H.block(idx_D, idx_D, N_D, N_D));
        U.block(idx_D, idx_D, N_D, N_D) = solver.eigenvectors();
        evals.segment(idx_D, N_D) = solver.eigenvalues();
    }

    // H_tilde = U^\dagger * H * U
    Eigen::MatrixXd H_tilde = U.transpose() * H * U;
    
    return {U, evals, H_tilde, partition};
}
