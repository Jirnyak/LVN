# DLVN Quantum Transport Simulation Engine

A C++23 / Eigen / OpenGL / ImGui simulation framework for nanoscale quantum electron transport based on the Driven Liouville-von Neumann (DLVN) methodology in the state representation (*Zelovich et al., J. Chem. Theory Comput. 2014, 10, 2927–2941*).

---

## Theoretical Framework & Methodology

Simulating open quantum systems where a central molecular junction is coupled to macroscopic leads requires specialized boundary treatment. Under the standard Liouville-von Neumann equation ($\frac{d\rho}{dt} = -\frac{i}{\hbar}[H, \rho]$), finite spatial representations of leads inevitably produce non-physical boundary reflections as propagating electronic wavepackets reach the finite edges.

The Driven Liouville-von Neumann (DLVN) framework eliminates artificial reflections by enforcing relaxation toward thermodynamic equilibrium across the extended lead spaces. This engine implements the complete formulation from first principles for arbitrary multi-terminal geometries, including standard two-terminal ($L-R$) and four-terminal ($L-R-U-D$) configurations.

### Mathematical Formulation

1. **Tight-Binding Hamiltonian (Site Representation)**  
   The total system comprising Left Lead ($L$), Right Lead ($R$), optional Up/Down Leads ($U, D$), and the central Extended Molecule ($EM$) is discretized on an atomic lattice. The Hamiltonian operator is defined as:
   $$\hat{H} = \sum_n \alpha_n \hat{c}_n^\dagger \hat{c}_n + \sum_{n \neq m} \beta_{nm} \hat{c}_n^\dagger \hat{c}_m$$
   where $\alpha_n$ is the on-site energy and $\beta_{nm}$ is the nearest-neighbor hopping integral. In matrix form:
   $$H = \begin{pmatrix} H_{EM} & V_{EM,L} & V_{EM,R} & V_{EM,U} & V_{EM,D} \\ V_{L,EM} & H_L & 0 & 0 & 0 \\ V_{R,EM} & 0 & H_R & 0 & 0 \\ V_{U,EM} & 0 & 0 & H_U & 0 \\ V_{D,EM} & 0 & 0 & 0 & H_D \end{pmatrix}$$
   Direct lead-to-lead coupling matrices are strictly zero ($V_{\alpha,\beta} = 0$ for $\alpha \neq \beta$).

2. **Transformation to the State Representation**  
   Thermodynamic properties such as temperature $T$ and chemical potential $\mu$ (bias voltage) are rigorously defined exclusively for stationary energy eigenstates via the Fermi-Dirac distribution. Applying damping directly in the spatial atomic site basis ($\{n\}$) is non-physical and violates the Pauli exclusion principle.
   
   To establish proper thermodynamic reservoirs, each isolated lead block is diagonalized: $H_\alpha U_\alpha = U_\alpha \tilde{H}_\alpha$, yielding diagonal eigenvalue matrices $\tilde{H}_\alpha$. Applying the global unitary transformation $U = U_{EM} \oplus \left(\bigoplus_\alpha U_\alpha\right)$ converts the Hamiltonian and density matrix into the state representation:
   $$\tilde{H} = U^\dagger H U, \quad \tilde{\rho} = U^\dagger \rho U$$

3. **Equation of Motion**  
   In the state representation, the target equilibrium density matrix $\tilde{\rho}^0_\alpha$ for each active lead $\alpha \in \mathcal{P} = \{L, R, U, D\}$ is diagonal with entries given by the Fermi-Dirac statistics:
   $$\tilde{\rho}^0_\alpha(E_k) = \frac{1}{1 + e^{(E_k - \mu_\alpha)/k_B T_\alpha}}$$
   The time evolution of the open quantum system is governed by:
   $$\frac{d\tilde{\rho}}{dt} = -\frac{i}{\hbar}[\tilde{H}, \tilde{\rho}] + \mathcal{D}[\tilde{\rho}]$$
   where the multi-lead dissipator $\mathcal{D}[\tilde{\rho}]$ acts block-wise across pairs of sub-spaces $(\alpha, \beta)$:
   $$\mathcal{D}[\tilde{\rho}]_{\alpha,\beta} = \begin{cases} -\Gamma (\tilde{\rho}_{\alpha,\alpha} - \tilde{\rho}^0_\alpha) & \text{if } \alpha = \beta \in \mathcal{P} \\ -\frac{1}{2}\Gamma \tilde{\rho}_{\alpha,EM} & \text{if } \alpha \in \mathcal{P}, \beta = EM \\ -\Gamma \tilde{\rho}_{\alpha,\beta} & \text{if } \alpha \neq \beta \in \mathcal{P} \\ 0 & \text{if } \alpha = \beta = EM \end{cases}$$

4. **Numerical Integration & Observable Current Calculation**  
   The differential equation is integrated forward in time using a fourth-order Runge-Kutta (RK4) scheme. At each time step $t$, the density matrix is transformed back to the spatial site representation ($\rho(t) = U \tilde{\rho}(t) U^\dagger$). The time-dependent electrical current $I_{n,n+1}(t)$ between adjacent lattice sites $n$ and $n+1$ inside the molecular junction is calculated via the probability continuity equation:
   $$I_{n,n+1}(t) = \frac{4 e^2}{\hbar} \beta_{n,n+1}^{\text{eV}} \text{Im}[\rho_{n,n+1}(t)]$$
   where $\beta_{n,n+1}^{\text{eV}}$ is expressed in electron-volts, yielding exact currents scaling directly in milliamperes ($\text{mA}$).

5. **Multi-Terminal Extension & Von Neumann Neighborhood Topology**  
   Because the DLVN formalism operates on abstract block sub-spaces, extending the central molecule coupling to four orthogonal leads ($\mathcal{P} = \{L, R, U, D\}$) realizes a four-point von Neumann neighborhood topology. 
   
   This four-terminal micro-solver provides a rigorous local computational unit for multiscale 2D quantum transport modeling. Instead of performing full $O(N^3)$ diagonalizations of macroscopic 2D lattice networks, the dynamic directional transmission properties obtained from this local four-point DLVN junction can be embedded into coarse-grained cellular automata or grid-based macroscopic simulations.

---

## Software Architecture & Implementation

The simulation engine is structured into modular layers focusing on zero-allocation runtime performance and clean mathematical abstraction:

```
src/
 ├── hamiltonian.h / .cpp    # Subsystem partitioning (LeadInfo, SystemPartition, TBModelParameters),
 │                           # tight-binding assembly, and unitary block diagonalization.
 ├── lvn_dynamics.h / .cpp   # Pre-allocated RK4Workspace containers and exact in-place
 │                           # evaluation of right-hand-side dynamics (dlvn_rhs_state_rep_inplace).
 └── main.cpp                # SDL2 / OpenGL / ImGui render loop, thread synchronization
                             # via std::mutex, and interactive parameter adjustment.
```

### Computational Optimization & Scaling

* **Zero-Allocation Hot Path**: The integration loop (`dlvn_rhs_state_rep_inplace`) executes without heap allocations. All Runge-Kutta slope buffers ($k_1, k_2, k_3, k_4$), scratch space, and pre-computed interaction blocks ($V_{\alpha, EM}$) are allocated once within `RK4Workspace::setup()`. Matrix evaluations use strictly in-place Eigen block operations (`.noalias()`).
* **Lead Discretization Requirements in 4-Terminal Mode**: While traditional tight-binding methods without volumetric damping require large lead sizes ($N_{\text{lead}} \ge 300$) to delay boundary reflections, the DLVN state representation relaxes every lead eigenstate with rate $\Gamma$. Consequently, the energy levels are broadened by $\hbar\Gamma$, creating a continuous effective density of states. Setting $N_{\text{lead}} = 80$ inside the configuration panel for four-terminal runs preserves the exact physical steady-state plateau and transient dynamics while accelerating the numerical integration by $10\times$ to $38\times$.

---

## Build and Execution (macOS / Linux)

### Prerequisites
* CMake 3.16 or higher
* C++23 compatible compiler (Clang or GCC)
* SDL2 (`brew install sdl2` on macOS or `sudo apt install libsdl2-dev` on Linux)

> **Note on Dependencies**: Manual installation of `Eigen`, `ImGui`, or `ImPlot` is not required. The build system uses CMake `FetchContent` to download and link exact pinned releases (Eigen 3.4.0, ImGui v1.91.5, ImPlot v0.16) directly inside `build/_deps/` during configuration.

### Build Commands

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

### Execution

```bash
./build/LVN
```

---

## Configuration & Simulation Parameters

### Molecule Grid
* **`N_Mx`**: Number of atomic sites along the horizontal axis of the molecular scattering region. Default: `6` (reproduces paper reference).
* **`N_My`**: Number of atomic sites along the vertical axis. Default: `1` (1D chain, degenerates to original linear topology). Set `N_My > 1` for 2D grid transport with quantum path interference. Each 2D lead couples site-by-site (1-to-1) to the corresponding edge of the scattering region.

### Lead Sizes (independent per lead)
* **`N_L`, `N_R`, `N_U`, `N_D`**: Number of atomic sites in each macroscopic lead reservoir. Recommended: `300` for 2-terminal runs; `80–100` for 4-terminal runs.
* **`N_ML`, `N_MR`, `N_MU`, `N_MD`**: Number of sites in each extended molecule buffer arm. These screen the lead boundary from the molecular junction.

### Simulation Parameters
* **`Bias Voltage (eV)`**: Chemical potential (Fermi level shift) configured independently for each lead (`bias_L`, `bias_R`, `bias_U`, `bias_D`), allowing arbitrary non-equilibrium multi-terminal configurations.
* **`Gamma (fs^-1)`**: Volumetric relaxation rate toward the Fermi-Dirac equilibrium distribution ($\Gamma$).
* **`Max Time (fs)`**: Total time duration of the Runge-Kutta integration window.

### Lead Toggles & Current Channels
* **`LEFT/RIGHT/UP/DOWN Lead`**: Enable or disable individual leads in the Hamiltonian.
* **`I_L`, `I_R`, `I_U`, `I_D`**: Select which transient current channels to compute and plot. Currents are measured at the lead-molecule junction bonds.

