#include <SDL.h>
#if defined(__APPLE__)
#include <OpenGL/gl3.h>
#else
#include <GLES3/gl3.h>
#endif
#include "imgui.h"
#include "imgui_impl_sdl2.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <string>
#include <cmath>
#include "hamiltonian.h"
#include "lvn_dynamics.h"

struct SimulationRun {
    std::vector<float> time;
    std::vector<float> current_L;
    std::vector<float> current_R;
    std::vector<float> current_U;
    std::vector<float> current_D;
    bool has_L{false};
    bool has_R{false};
    bool has_U{false};
    bool has_D{false};
    double gamma;
    double temperature;
    std::string label;
};

struct SimulationState {
    std::vector<SimulationRun> history;
    SimulationRun current_run;
    std::atomic<bool> is_running{false};
    std::atomic<float> progress{0.0f};
    std::mutex data_mutex;
    
    void start_new_run(const TBModelParameters& p, double g, double t, bool hl, bool hr, bool hu, bool hd) {
        std::lock_guard<std::mutex> lock(data_mutex);
        if (current_run.time.size() > 0) {
            history.push_back(current_run);
        }
        current_run = SimulationRun();
        current_run.gamma = g;
        current_run.temperature = t;
        current_run.has_L = hl;
        current_run.has_R = hr;
        current_run.has_U = hu;
        current_run.has_D = hd;
        
        std::string lbl = "G=" + std::to_string(g).substr(0, 5) + ", T=" + std::to_string((int)t) + "K";
        if (hl) lbl += " L:" + std::to_string(p.bias_L).substr(0, 5);
        if (hr) lbl += " R:" + std::to_string(p.bias_R).substr(0, 5);
        if (hu) lbl += " U:" + std::to_string(p.bias_U).substr(0, 5);
        if (hd) lbl += " D:" + std::to_string(p.bias_D).substr(0, 5);
        current_run.label = lbl;
        progress = 0.0f;
    }
    
    void clear_history() {
        std::lock_guard<std::mutex> lock(data_mutex);
        history.clear();
        current_run = SimulationRun();
        progress = 0.0f;
    }
};

static double calculate_current_between(const Eigen::MatrixXcd& rho_site, int i, int j, double beta) {
    const double factor = 0.97358;
    return factor * beta * rho_site(i, j).imag();
}

struct PlotChannels {
    bool L{true};
    bool R{false};
    bool U{false};
    bool D{false};
};

void run_simulation(TBModelParameters params, double gamma, double temperature, double dt, double max_time_fs, PlotChannels channels, SimulationState* state) {
    state->is_running = true;
    state->start_new_run(params, gamma, temperature, channels.L, channels.R, channels.U, channels.D);

    Eigen::MatrixXd H_site = build_hamiltonian(params);
    StateRepresentation rep = compute_state_representation(H_site, params);
    
    int N = H_site.rows();

    RK4Workspace ws;
    ws.setup(rep, temperature);
    
    Eigen::MatrixXcd rho_tilde = Eigen::MatrixXcd::Zero(N, N);
    for (size_t i = 0; i < ws.leads.size(); ++i) {
        const auto& lw = ws.leads[i];
        rho_tilde.block(lw.start, lw.start, lw.size, lw.size) = lw.rho0;
    }

    // dt is now passed as an argument 
    int steps = (int)(max_time_fs / dt);
    
    int idx_EM = rep.partition.EM_start;
    int N_ML = params.has_L ? params.N_ML * params.N_My : 0;
    int N_M_total = params.N_Mx * params.N_My;
    int N_MR = params.has_R ? params.N_MR * params.N_My : 0;
    int N_MU = (params.has_U && params.N_U > 0) ? params.N_MU * params.N_Mx : 0;
    int N_MD = (params.has_D && params.N_D > 0) ? params.N_MD * params.N_Mx : 0;

    // Only need back-transform if any current channel is requested
    bool need_current = channels.L || channels.R || channels.U || channels.D;

    for (int step = 0; step <= steps; ++step) {
        if (!state->is_running) break; 
        
        rk4_step_state_rep(rho_tilde, gamma, dt, ws);
        
        if (need_current && step % 5 == 0) {
            // Block-diagonal U → exact EM back-transform
            Eigen::MatrixXcd rho_tilde_EM = rho_tilde.block(idx_EM, idx_EM, ws.EM_size, ws.EM_size);
            Eigen::MatrixXd U_EM = rep.U.block(idx_EM, idx_EM, ws.EM_size, ws.EM_size);
            Eigen::MatrixXcd rho_site_EM = U_EM * rho_tilde_EM * U_EM.transpose();
            
            std::lock_guard<std::mutex> lock(state->data_mutex);
            state->current_run.time.push_back(step * dt);

            // All currents summed over full lead-molecule edge bonds (1-to-1)
            // No normalization needed for 2D leads
            if (channels.L && params.has_L && N_ML > 0) {
                double c_L = 0.0;
                for (int iy = 0; iy < params.N_My; ++iy) {
                    int idx_ml = (params.N_ML - 1) * params.N_My + iy;
                    int idx_m = N_ML + iy;
                    c_L += calculate_current_between(rho_site_EM, idx_ml, idx_m, params.beta_LM);
                }
                state->current_run.current_L.push_back((float)c_L);
            }
            if (channels.R && params.has_R && N_MR > 0) {
                double c_R = 0.0;
                for (int iy = 0; iy < params.N_My; ++iy) {
                    int idx_m = N_ML + (params.N_Mx - 1) * params.N_My + iy;
                    int idx_mr = N_ML + N_M_total + iy;
                    c_R += calculate_current_between(rho_site_EM, idx_m, idx_mr, params.beta_MR);
                }
                state->current_run.current_R.push_back((float)c_R);
            }
            if (channels.U && params.has_U && N_MU > 0) {
                double c_U = 0.0;
                int first_MU = N_ML + N_M_total + N_MR;
                for (int ix = 0; ix < params.N_Mx; ++ix) {
                    int idx_mu = first_MU + ix;
                    int idx_m = N_ML + ix * params.N_My; // iy = 0
                    c_U += calculate_current_between(rho_site_EM, idx_mu, idx_m, params.beta_UM);
                }
                state->current_run.current_U.push_back((float)c_U);
            }
            if (channels.D && params.has_D && N_MD > 0) {
                double c_D = 0.0;
                int first_MD = N_ML + N_M_total + N_MR + N_MU;
                for (int ix = 0; ix < params.N_Mx; ++ix) {
                    int idx_md = first_MD + ix;
                    int idx_m = N_ML + ix * params.N_My + (params.N_My - 1); // iy = N_My - 1
                    c_D += calculate_current_between(rho_site_EM, idx_md, idx_m, params.beta_DM);
                }
                state->current_run.current_D.push_back((float)c_D);
            }
            state->progress = (float)step / steps;
        }
    }
    
    if (state->is_running) {
        state->progress = 1.0f;
    }
    state->is_running = false;
}

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;
    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0) {
        std::cerr << "Error: " << SDL_GetError() << std::endl;
        return -1;
    }

#if defined(__APPLE__)
    const char* glsl_version = "#version 150";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG); 
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
#else
    const char* glsl_version = "#version 300 es";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_ES);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
#endif

    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_Window* window = SDL_CreateWindow("DLVN Dynamics Simulator", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720, window_flags);
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, gl_context);
    SDL_GL_SetSwapInterval(1); 

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init(glsl_version);

    TBModelParameters params;
    params.has_L = true;
    params.has_R = true;
    params.N_L = 300; 
    params.N_R = 300;
    params.N_ML = 50;
    params.N_MR = 50;
    params.N_Mx = 6;   // Molecule grid width  (1D chain for paper default)
    params.N_My = 1;   // Molecule grid height (1 = 1D chain)
    params.alpha_L = 0.0;
    params.alpha_M = 0.0;
    params.alpha_R = 0.0;
    params.beta_L = -0.2;
    params.beta_M = -0.2;
    params.beta_R = -0.2;
    params.beta_LM = -0.2;
    params.beta_MR = -0.2;

    params.has_U = false;
    params.N_U = 300;
    params.N_MU = 50;
    params.alpha_U = 0.0;
    params.beta_U = -0.2;
    params.beta_UM = -0.2;

    params.has_D = false;
    params.N_D = 300;
    params.N_MD = 50;
    params.alpha_D = 0.0;
    params.beta_D = -0.2;
    params.beta_DM = -0.2;

    double gamma = 0.000;
    double max_time = 3000.0;
    double temperature = 0.0;
    double dt = 1.0;

    // Which current channels to compute and plot
    PlotChannels channels;

    SimulationState sim_state;
    std::thread sim_thread;

    bool done = false;
    while (!done) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT)
                done = true;
            if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window))
                done = true;
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 420), ImGuiCond_FirstUseEver);
        ImGui::Begin("DLVN Model Control");
        
        ImGui::Text("Hamiltonian Leads");
        ImGui::Checkbox("LEFT Lead (L)", &params.has_L); ImGui::SameLine();
        ImGui::Checkbox("RIGHT Lead (R)", &params.has_R);
        ImGui::Checkbox("UP Lead (U)", &params.has_U); ImGui::SameLine();
        ImGui::Checkbox("DOWN Lead (D)", &params.has_D);
        
        ImGui::Separator();
        ImGui::Text("Hamiltonian Parameters");
        ImGui::Text("Molecule Grid");
        ImGui::InputInt("N_Mx", &params.N_Mx);
        ImGui::InputInt("N_My", &params.N_My);
        if (params.N_Mx < 1) params.N_Mx = 1;
        if (params.N_My < 1) params.N_My = 1;
        
        if (ImGui::TreeNode("Tight-Binding Parameters (eV)")) {
            ImGui::Text("Molecule");
            ImGui::InputDouble("alpha_M", &params.alpha_M);
            ImGui::InputDouble("beta_M", &params.beta_M);
            ImGui::Separator();
            ImGui::Text("Couplings (Lead-Molecule)");
            ImGui::InputDouble("beta_LM", &params.beta_LM);
            ImGui::InputDouble("beta_MR", &params.beta_MR);
            ImGui::InputDouble("beta_UM", &params.beta_UM);
            ImGui::InputDouble("beta_DM", &params.beta_DM);
            ImGui::Separator();
            ImGui::Text("Lead L");
            ImGui::InputDouble("alpha_L", &params.alpha_L);
            ImGui::InputDouble("beta_L", &params.beta_L);
            ImGui::Separator();
            ImGui::Text("Lead R");
            ImGui::InputDouble("alpha_R", &params.alpha_R);
            ImGui::InputDouble("beta_R", &params.beta_R);
            ImGui::Separator();
            ImGui::Text("Lead U");
            ImGui::InputDouble("alpha_U", &params.alpha_U);
            ImGui::InputDouble("beta_U", &params.beta_U);
            ImGui::Separator();
            ImGui::Text("Lead D");
            ImGui::InputDouble("alpha_D", &params.alpha_D);
            ImGui::InputDouble("beta_D", &params.beta_D);
            ImGui::TreePop();
        }
        
        ImGui::Separator();
        ImGui::Text("Lead Configuration");
        if (params.has_L) {
            ImGui::Text("LEFT Lead");
            ImGui::InputInt("N_L", &params.N_L);
            ImGui::InputInt("N_ML", &params.N_ML);
            ImGui::InputDouble("Bias L (eV)", &params.bias_L);
        }
        if (params.has_R) {
            ImGui::Separator();
            ImGui::Text("RIGHT Lead");
            ImGui::InputInt("N_R", &params.N_R);
            ImGui::InputInt("N_MR", &params.N_MR);
            ImGui::InputDouble("Bias R (eV)", &params.bias_R);
        }
        if (params.has_U) {
            ImGui::Separator();
            ImGui::Text("UP Lead");
            ImGui::InputInt("N_U", &params.N_U);
            ImGui::InputInt("N_MU", &params.N_MU);
            ImGui::InputDouble("Bias U (eV)", &params.bias_U);
        }
        if (params.has_D) {
            ImGui::Separator();
            ImGui::Text("DOWN Lead");
            ImGui::InputInt("N_D", &params.N_D);
            ImGui::InputInt("N_MD", &params.N_MD);
            ImGui::InputDouble("Bias D (eV)", &params.bias_D);
        }
        
        ImGui::Separator();
        ImGui::Text("Simulation Parameters");
        ImGui::InputDouble("Gamma (fs^-1)", &gamma);
        ImGui::InputDouble("Max Time (fs)", &max_time);
        ImGui::InputDouble("Temperature (K)", &temperature);
        if (temperature < 0.0) temperature = 0.0;
        ImGui::InputDouble("Time Step (dt) (fs)", &dt);
        if (dt <= 0.0) dt = 0.1;
        
        ImGui::Separator();
        
        if (sim_state.is_running) {
            ImGui::Text("Simulation running...");
            ImGui::ProgressBar(sim_state.progress, ImVec2(0.0f, 0.0f));
            if (ImGui::Button("Cancel Simulation")) {
                sim_state.is_running = false;
                if (sim_thread.joinable()) {
                    sim_thread.join();
                }
            }
        } else {
            if (ImGui::Button("Run Simulation")) {
                if (sim_thread.joinable()) {
                    sim_thread.join();
                }
                sim_thread = std::thread(run_simulation, params, gamma, temperature, dt, max_time, channels, &sim_state);
            }
        }
        
        if (ImGui::Button("Clear Graph")) {
            sim_state.clear_history();
        }
        
        ImGui::Separator();
        ImGui::Text("Current Channels");
        ImGui::Checkbox("I_L", &channels.L); ImGui::SameLine();
        ImGui::Checkbox("I_R", &channels.R); ImGui::SameLine();
        ImGui::Checkbox("I_U", &channels.U); ImGui::SameLine();
        ImGui::Checkbox("I_D", &channels.D);
        
        ImGui::End();

        ImGui::SetNextWindowPos(ImVec2(320, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(900, 600), ImGuiCond_FirstUseEver);
        ImGui::Begin("Transient Current");
        
        if (ImPlot::BeginPlot("##CurrentPlot", ImVec2(-1, -1))) {
            ImPlot::SetupAxes("Time (fs)", "Current (mA)");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, max_time, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.03, 0.04, ImPlotCond_Once);
            
            ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

            std::lock_guard<std::mutex> lock(sim_state.data_mutex);
            
            // Helper lambda to draw a single run
            auto draw_run = [](const SimulationRun& run) {
                if (run.time.empty()) return;
                if (run.has_L && !run.current_L.empty()) {
                    std::string lbl = run.label + " [L]";
                    int n = (int)std::min(run.time.size(), run.current_L.size());
                    ImPlot::PlotLine(lbl.c_str(), run.time.data(), run.current_L.data(), n);
                }
                if (run.has_R && !run.current_R.empty()) {
                    std::string lbl = run.label + " [R]";
                    int n = (int)std::min(run.time.size(), run.current_R.size());
                    ImPlot::PlotLine(lbl.c_str(), run.time.data(), run.current_R.data(), n);
                }
                if (run.has_U && !run.current_U.empty()) {
                    std::string lbl = run.label + " [U]";
                    int n = (int)std::min(run.time.size(), run.current_U.size());
                    ImPlot::PlotLine(lbl.c_str(), run.time.data(), run.current_U.data(), n);
                }
                if (run.has_D && !run.current_D.empty()) {
                    std::string lbl = run.label + " [D]";
                    int n = (int)std::min(run.time.size(), run.current_D.size());
                    ImPlot::PlotLine(lbl.c_str(), run.time.data(), run.current_D.data(), n);
                }
            };

            // Draw past runs
            for (size_t i = 0; i < sim_state.history.size(); ++i) {
                draw_run(sim_state.history[i]);
            }
            
            // Draw current run
            draw_run(sim_state.current_run);

            ImPlot::EndPlot();
        }
        
        ImGui::End();

        ImGui::Render();
        glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);
    }

    if (sim_state.is_running) {
        sim_state.is_running = false;
    }
    if (sim_thread.joinable()) {
        sim_thread.join();
    }

    ImPlot::DestroyContext();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
