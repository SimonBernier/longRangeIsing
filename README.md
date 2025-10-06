# Long-Range Transverse Field Ising Model - Moving Front Quenches

Tensor network simulations of the transverse field Ising model with long-range algebraic interactions using ITensor. This project investigates spatiotemporal quenches and quantum dynamics in systems with power-law decaying interactions J(r) ~ 1/r<sup>α</sup>, relevant to experimental quantum simulators based on Rydberg atoms and trapped ions.

## 🔬 Overview

This repository contains the numerical implementations used to generate all figures in our Physical Review B publication on spatiotemporal quenches in long-range Hamiltonians. The work extends previous studies on short-range models to realistic experimental systems where interactions decay algebraically with distance.

### Physical Model

The long-range transverse field Ising Hamiltonian:

```
H = -Σᵢⱼ J(rᵢⱼ) σˣᵢ σˣⱼ - h Σᵢ σᶻᵢ
```

where `J(r) ~ 1/r^α` is the long-range coupling strength, h is the transverse field, and α controls the range of interactions.

### Key Physics

**Long-range interactions fundamentally change:**
- Critical behavior and universality classes
- Dynamical critical exponents (z ≠ 1 for α < 3)
- Speed of information propagation
- Efficiency of spatiotemporal quenches for ground state preparation

---

## 📂 Project Structure

### Core Modules

#### 1. **gapLR** - Energy Gap Calculations ⭐
Calculates the energy gap between ground and first excited states in the long-range TFI model.

**Key Applications:**
- Determining critical points via gap scaling
- Extracting critical exponents
- **Validating exponential sum approximation** of algebraic interactions
- Testing the accuracy of the sum-of-exponentials representation

**Technical Note:** The gap calculations verify that representing `1/r^α` as a sum of exponentials (using optimized parameters from `input_alpha_...` files) accurately reproduces the physics of true algebraic interactions.

---

#### 2. **vCrit** - Speed of Excitations
Determines the effective speed of light in the long-range quantum system.

**Method:**
- Introduces local perturbation to ground state
- Evolves using 4th order TDVP
- Tracks von Neumann entropy spreading
- Extracts propagation velocity from entanglement light cone

**Physical Significance:** In long-range models, the effective speed of excitations depends on α and can differ from the nearest-neighbor case, affecting the optimal quench protocol.

---

#### 3. **superluminal** - Moving Front Quenches (Main Focus) 🌟
Simulates inhomogeneous quenches with moving quench fronts in the long-range TFI model.

**Protocol:**
- Quench front propagates at constant velocity v
- Initial state: Ground state in gapped phase
- Final state: Critical ground state
- Evolution: 4th order Time-Dependent Variational Principle (TDVP)

**Observables:**
- Local and global energy density evolution
- Spin correlation functions
- Von Neumann entanglement entropy dynamics

**Key Results:**
- For α ≳ 3 (where z = 1), optimal cooling occurs when v ≈ c
- For α < 3, the long-range nature modifies optimal quench protocols
- Spatiotemporal quenches remain efficient for ground state preparation

---

#### 4. **critical** - Critical Ground State Properties
Calculates equilibrium properties at the critical point, including:
- Ground state correlations
- Critical exponents
- Scaling behavior

---

#### 5. **staticFront** - Static Front Analysis
Investigates ground state properties with a static (non-moving) quench front.

**Purpose:** Provides baseline comparison for moving front dynamics.

---

## 🛠️ Technical Implementation

### Framework
**ITensor C++ Library** - Tensor network calculations optimized for 1D systems

### Key Innovation: Sum-of-Exponentials Representation

**Challenge:** True algebraic interactions J(r) ~ 1/r<sup>α</sup> are computationally expensive in tensor network methods.

**Solution:** Approximate power-law decay as sum of exponentials:

```
J(r) ~ 1/r^α  ≈  Σᵢ Aᵢ exp(-r/ξᵢ)
```

**Implementation:**
- Parameters {Aᵢ, ξᵢ} obtained via **nonlinear optimization**
- Stored in `input_alpha_...` files for different α values
- **TFIsingLR.h** header file implements the long-range Hamiltonian using these optimized parameters
- Validation: `gapLR` module confirms approximation accuracy

**Why This Matters:**
- Enables efficient Matrix Product State (MPS) representation
- Preserves physical properties (critical points, exponents, dynamics)
- Makes long-range simulations tractable while maintaining accuracy

### Algorithms
- 4th order Time-Dependent Variational Principle (TDVP)
- Matrix Product State (MPS) evolution
- Finite-size scaling
- Entanglement entropy calculations

**Language:** C++

---

## 📊 Physical Insights

This work demonstrates:
- **Spatiotemporal quenches** work efficiently even with long-range interactions
- **Optimal quench velocity** depends on the dynamical critical exponent z
- For α ≳ 3: Lorentz-like behavior (z = 1), optimal v ≈ c
- For α < 3: Modified dynamics, non-Lorentzian behavior
- **Quantum simulators** (Rydberg atoms, trapped ions) naturally realize these long-range models

---

## 🔬 Experimental Relevance

This work directly applies to:
- **Rydberg atom arrays:** Van der Waals interactions (α = 6) or dipolar (α = 3)
- **Trapped ion systems:** Coulomb interactions with tunable range
- **Polar molecules:** Electric dipole-dipole interactions
- **Quantum simulation** of many-body dynamics

---

## 📄 Publication

This code was used to generate all figures in:

**[Spatiotemporal Quenches in Long-Range Hamiltonians](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.108.024310)**  
Simon Bernier and Kartiek Agarwal  
*Physical Review B* **108**, 024310 (2023)  
[arXiv:2212.07499](https://arxiv.org/abs/2212.07499)

---

## 🎓 Academic Context

This project extends our previous work on 2D nearest-neighbor models to experimentally relevant systems with long-range interactions, bridging:
- Tensor network methods
- Quantum simulation platforms
- Non-equilibrium quantum dynamics
- Critical phenomena beyond nearest neighbors

---

## 🔗 Dependencies

- **ITensor library** (C++)
- C++11 compatible compiler
- LAPACK/BLAS libraries
- Optimized `input_alpha_...` parameter files (included in repository)

---

## 💡 Key Files

### TFIsingLR.h
Core header file implementing the long-range transverse field Ising Hamiltonian using the sum-of-exponentials approximation.

**Usage:** This file must be properly configured with the appropriate `input_alpha_...` parameter file for your chosen interaction range α.

### input_alpha_... Files
Pre-computed optimization results for different values of α, containing:
- Exponential amplitudes {Aᵢ}
- Decay lengths {ξᵢ}
- Obtained via nonlinear least-squares fitting to 1/r<sup>α</sup>

---

## 📈 Computational Strategy

1. **Choose interaction range** α (determines physics regime)
2. **Load optimized parameters** from corresponding `input_alpha_...` file
3. **Configure TFIsingLR.h** with these parameters
4. **Run simulations** (equilibrium or time evolution)
5. **Validate approximation** using gapLR module if needed

---

## 🔗 Related Work

This repository complements our 2D work:
- [ising2d](https://github.com/SimonBernier/ising2d) - Short-range 2D TFI model
- Both projects study spatiotemporal quenches but in different geometries and interaction ranges

---

## 📚 Physics Background

**Dynamical Critical Exponent (z):**
- Controls how time and space scale at criticality: τ ~ ξ<sup>z</sup>
- Short-range models: z = 1 (Lorentz invariance)
- Long-range models: z depends on α
  - α > 3: z = 1 (short-range-like)
  - α < 3: z < 1 (true long-range behavior)

**Kibble-Zurek Mechanism:**
Predicts defect scaling during quenches. Spatiotemporal quenches can outperform uniform quenches by respecting causal structure.

---

## 📧 Contact

**Simon Bernier**
- Email: simon.bernier@mail.mcgill.ca
- LinkedIn: [simon-bernier-6701a9285](https://www.linkedin.com/in/simon-bernier-6701a9285)

---

## 📝 Citation

If you use this code or build upon this work, please cite:

```bibtex
@article{bernier2023spatiotemporal,
  title={Spatiotemporal Quenches in Long-Range Hamiltonians},
  author={Bernier, Simon and Agarwal, Kartiek},
  journal={Physical Review B},
  volume={108},
  pages={024310},
  year={2023},
  publisher={American Physical Society},
  doi={10.1103/PhysRevB.108.024310}
}
```

---

## 🎯 Future Directions

- Extension to 2D long-range models
- Floquet engineering with periodic driving
- Machine learning for parameter optimization
- Connections to near-term quantum devices

---

*This project demonstrates expertise in: computational quantum many-body physics, tensor networks, long-range interactions, nonlinear optimization, C++ programming, and quantum simulation relevant to experimental platforms.*
