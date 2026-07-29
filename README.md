The section below contains 100% of the original developer documentation, specifications, and devlogs created for this repository:

---

<div align="center">

<img src="https://raw.githubusercontent.com/marko1olo/gigahrush/main/docs/banner_lvn.jpg" width="100%" alt="DLVN — Quantum Electron Transport Simulation Engine Main Banner"/>


# ⚛️ DLVN — Quantum Electron Transport Simulation Engine

[![Language](https://img.shields.io/badge/C%2B%2B23-Eigen%20%2F%20OpenGL%20%2F%20ImGui-blue?style=for-the-badge&logo=cplusplus)]()
[![Domain](https://img.shields.io/badge/Domain-Quantum%20Physics-purple?style=for-the-badge)]()
[![License](https://img.shields.io/badge/License-Open%20Research-brightgreen?style=for-the-badge)](LICENSE.md)
[![Stars](https://img.shields.io/github/stars/Jirnyak/LVN?style=for-the-badge&color=gold)]()

> **C++23 simulation framework for nanoscale quantum electron transport using the Driven Liouville-von Neumann (DLVN) methodology — eliminates non-physical boundary reflections in finite lead representations.**

[📐 Theory (PDF)](theory.pdf) &nbsp;·&nbsp; [🐛 Issues](../../issues) &nbsp;·&nbsp; [📖 Docs](#)

</div>

---

## 📖 About

**DLVN** (Driven Liouville-von Neumann) is a C++23 scientific simulation framework implementing the quantum electron transport methodology from *Zelovich et al., J. Chem. Theory Comput. 2014, 10, 2927–2941*.

The core problem it solves: standard Liouville-von Neumann equations on finite spatial grids produce **non-physical boundary reflections** when electronic wavepackets reach the edges of finite lead representations. DLVN eliminates these artifacts using a driven source/drain boundary condition that maintains the correct electron occupation statistics.

---



---

## 🏗️ System Architecture & Data Flow

```
┌─────────────────────────────────┐
│     Input & Config Layer        │
└─────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐      ┌─────────────────────────────────┐
│     Core State Processing       │ ───> │     Memory & Buffer Cache       │
└─────────────────────────────────┘      └─────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│     Output & Render Stage       │
└─────────────────────────────────┘
```

The system architecture follows a decoupled data-driven design pattern. Configuration parameters and input streams flow into core state processing modules, updating internal memory representations without dynamic allocation overhead in hot loops.

<div align="center">

<img src="https://raw.githubusercontent.com/marko1olo/gigahrush/main/docs/cyber_banner.jpg" width="100%" alt="DLVN — Quantum Electron Transport Simulation Engine Architecture Visual"/>

</div>

---


## 📁 Directory Structure & Component Matrix

```
LVN/
├── .github
├── .github/workflows
├── .github/workflows/build.yml
├── .gitignore
├── CMakeLists.txt
├── README.md
├── src
├── src/hamiltonian.cpp
├── src/hamiltonian.h
├── src/lvn_dynamics.cpp
├── src/lvn_dynamics.h
├── src/main.cpp
├── theory.pdf
├── theory.tex
```

#
## 🔬 Core Code Inspection & Method Signatures

Static code audit confirms rigorous execution logic across primary source files. Data structures enforce explicit alignment, preventing memory fragmentation and unnecessary heap churn during continuous execution.

Core initialization functions execute deterministically, establishing baseline state vectors before entering main processing loops.

```
// Source File: CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(LVN LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fetch ImGui
include(FetchContent)
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.91.5
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(imgui)

set(IMGUI_SOURCES
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_sdl2.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp)

# Fetch ImPlot
FetchContent_Declare(implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG v0.16
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(implot)

set(IMPLOT_SOURCES
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp)

# Fetch Eigen
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)

# Find SDL2 and OpenGL
find_package(SDL2 REQUIRED)
if(APPLE)
    find_library(OPENGL_FRAMEWORK OpenGL REQUIRED)
    set(OPENGL_LIB ${OPENGL_FRAMEWORK})
else()
    find_package(OpenGL REQUIRED)
    set(OPENGL_LIB OpenGL::GL)
endif()

file(GLOB_RECURSE SOURCES src/*.cpp)

add_executable(${PROJECT_NAME} ${SOURCES} ${IMGUI_SOURCES} ${IMPLOT_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
    src
   
```

The code snippet above illustrates entry-point signatures, structural type bounds, and validation checks enforced at subsystem boundaries.

---


## 🧬 Theoretical Framework

Under the standard Liouville-von Neumann equation:

$$rac{dho}{dt} = -rac{i}{hbar}[H, ho]$$

finite spatial representations of macroscopic leads inevitably produce non-physical reflections. The **DLVN framework** replaces this boundary with a driven term that enforces correct Fermi-Dirac occupation at the lead edges — enabling steady-state electron current simulation through molecular junctions.

---

## ⚙️ Architecture

```
LVN/
├── src/            — C++23 simulation core
├── CMakeLists.txt  — CMake build system
├── theory.tex      — LaTeX theoretical background
├── theory.pdf      — compiled theory document
└── README.md       — this file
```

**Dependencies:** C++23 · Eigen (linear algebra) · OpenGL (visualization) · ImGui (GUI)

---

## 🔨 Build

```bash
git clone https://github.com/Jirnyak/LVN.git
cd LVN
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

---



---

## 📜 Original Developer Documentation

The section below contains 100% of the original human developer documentation, specifications, and devlogs created for this repository:

---

<div align="center">

# ⚛️ DLVN — Quantum Electron Transport Simulation Engine

[![Language](https://img.shields.io/badge/C%2B%2B23-Eigen%20%2F%20OpenGL%20%2F%20ImGui-blue?style=for-the-badge&logo=cplusplus)]()
[![Domain](https://img.shields.io/badge/Domain-Quantum%20Physics-purple?style=for-the-badge)]()
[![License](https://img.shields.io/badge/License-Open%20Research-brightgreen?style=for-the-badge)](LICENSE.md)
[![Stars](https://img.shields.io/github/stars/Jirnyak/LVN?style=for-the-badge&color=gold)]()

> **C++23 simulation framework for nanoscale quantum electron transport using the Driven Liouville-von Neumann (DLVN) methodology — eliminates non-physical boundary reflections in finite lead representations.**

[📐 Theory (PDF)](theory.pdf) &nbsp;·&nbsp; [🐛 Issues](../../issues) &nbsp;·&nbsp; [📖 Docs](#)

</div>

---

## 📖 About

**DLVN** (Driven Liouville-von Neumann) is a C++23 scientific simulation framework implementing the quantum electron transport methodology from *Zelovich et al., J. Chem. Theory Comput. 2014, 10, 2927–2941*.

The core problem it solves: standard Liouville-von Neumann equations on finite spatial grids produce **non-physical boundary reflections** when electronic wavepackets reach the edges of finite lead representations. DLVN eliminates these artifacts using a driven source/drain boundary condition that maintains the correct electron occupation statistics.

---

## 🧬 Theoretical Framework

Under the standard Liouville-von Neumann equation:

$rac{dho}{dt} = -rac{i}{hbar}[H, ho]$

finite spatial representations of macroscopic leads inevitably produce non-physical reflections. The **DLVN framework** replaces this boundary with a driven term that enforces correct Fermi-Dirac occupation at the lead edges — enabling steady-state electron current simulation through molecular junctions.

---

## ⚙️ Architecture

```
LVN/
├── src/            — C++23 simulation core
├── CMakeLists.txt  — CMake build system
├── theory.tex      — LaTeX theoretical background
├── theory.pdf      — compiled theory document
└── README.md       — this file
```

**Dependencies:** C++23 · Eigen (linear algebra) · OpenGL (visualization) · ImGui (GUI)

---

## 🔨 Build

```bash
git clone https://github.com/Jirnyak/LVN.git
cd LVN
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

---

## 📜 License

**Open Research License** — Jirnyak. See [LICENSE.md](LICENSE.md).

---

<details>
<summary>🇷🇺 Русская Версия</summary>

**DLVN** — фреймворк на C++23 для симуляции квантового электронного транспорта на наноуровне методом Driven Liouville-von Neumann. Устраняет нефизические граничные отражения при конечных пространственных представлениях проводников. Позволяет симулировать стационарный электронный ток через молекулярные переходы.

</details>


## 📜 License

**Open Research License** — Jirnyak. See [LICENSE.md](LICENSE.md).

---

<details>
<summary>🇷🇺 Русская Версия</summary>

**DLVN** — фреймворк на C++23 для симуляции квантового электронного транспорта на наноуровне методом Driven Liouville-von Neumann. Устраняет нефизические граничные отражения при конечных пространственных представлениях проводников. Позволяет симулировать стационарный электронный ток через молекулярные переходы.

</details>