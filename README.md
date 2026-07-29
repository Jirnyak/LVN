<div align="center">

<img src="https://raw.githubusercontent.com/marko1olo/gigahrush/main/docs/banner_lvn.jpg" width="100%" alt="DLVN Banner"/>

# ⚛️ DLVN — Quantum Electron Transport Simulation Engine

[![Language](https://img.shields.io/badge/Language-C%2B%2B23%20%2F%20Eigen-blue?style=for-the-badge&logo=cplusplus)]()
[![Domain](https://img.shields.io/badge/Domain-Quantum%20Physics-purple?style=for-the-badge)]()
[![License](https://img.shields.io/badge/License-Open%20Research-brightgreen?style=for-the-badge)](LICENSE.md)

> **C++23 quantum electron transport simulation framework using Driven Liouville-von Neumann methodology — eliminates boundary reflections in nanoscale molecular junctions.**

</div>

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

## 📜 License

**Open Research License** — Jirnyak. See [LICENSE.md](LICENSE.md).

---

<details>
<summary>🇷🇺 Русская Версия</summary>

**DLVN** — фреймворк на C++23 для симуляции квантового электронного транспорта на наноуровне методом Driven Liouville-von Neumann. Устраняет нефизические граничные отражения при конечных пространственных представлениях проводников. Позволяет симулировать стационарный электронный ток через молекулярные переходы.

</details>
