# 🔥 Thermodynamic-PCA Simulation

<p align="center">
  <img src="docs/banner_heatmap.svg" width="80%">
</p>

---

## 📘 Overview

We simulate *n* artificial chemical species, each with *p* discrete energy
levels *E<sub>ij</sub>* drawn from a chosen distribution.  
Using textbook statistical-mechanics we compute, per species *i*:

| Symbol | Meaning | Units |
| ------ | ------- | ----- |
| Z      | Canonical partition function              | dimensionless |
| P<sub>ij</sub> | Boltzmann probability of level *j* | – |
| ⟨E⟩    | Average energy                           | eV |
| S      | Entropy (Shannon / Gibbs)                | eV K⁻¹ |
| F      | Helmholtz free energy                    | eV |
| %NA    | % “inaccessible” micro-states (`P < 10⁻⁶`)| % |

The resulting feature matrix (species × thermo-features) is then inspected with
Principal Component Analysis.  A companion 3-D visualiser shows both *scores*
(points) and *loadings* (vectors) for intuitive interpretation.

---

## 🧠 Conceptual Mapping

|  ML / Data Science Concept |  In This Project                                     |
| -------------------------- | ---------------------------------------------------- |
| Samples                    | Chemical species                                     |
| Features                   | Thermodynamic quantities (⟨E⟩, S, Z, F, %NA)         |
| Labels (future)            | Phases / clusters revealed by PCA or clustering      |
| Controlled variations      | Temperature, chemical potential, energy-level shape |

---

## 🔧 Current Tech Stack

* **Language** MATLAB (R2023a +)  
  *Python port on the roadmap—see “Future work”.*
* **Math** Statistical mechanics · Thermodynamics · Linear algebra  
* **External libs** None beyond base MATLAB (no Toolboxes required)

---
```
1. Partition Function and Temperature Dependence
Expectation: The canonical partition function ( Z ) and derived properties like average energy ( \langle E \rangle ) and entropy ( S ) should exhibit strong temperature dependence.
Observation: As temperature increases, the partition function should increase, leading to higher average energy and entropy. This reflects the increased population of higher energy states.
PCA Insight: PCA should reveal that temperature is a dominant factor in the first principal component, as it significantly influences multiple thermodynamic properties.
2. Boltzmann Distribution
Expectation: The Boltzmann probabilities ( P_i ) should follow the exponential distribution ( P_i \propto e^{-\beta E_i} ), where ( \beta = 1/(k_B T) ).
Observation: At low temperatures, only the lowest energy states should be significantly populated, while at high temperatures, higher energy states become more accessible.
PCA Insight: The distribution of Boltzmann probabilities across energy levels should cluster in a way that PCA can distinguish based on temperature.
3. Entropy and Disorder
Expectation: Entropy ( S ) should increase with temperature, reflecting greater disorder in the system.
Observation: Systems with more energy levels or broader distributions of energy levels should exhibit higher entropy.
PCA Insight: PCA should identify entropy as a key factor in distinguishing different species or conditions, particularly in the context of phase transitions or clustering.
4. Free Energy and Stability
Expectation: The Helmholtz free energy ( F ) should decrease with temperature for systems where entropy increases more rapidly than energy.
Observation: Systems with lower free energy are more stable. Changes in chemical potential or energy-level shape should affect the free energy landscape.
PCA Insight: Free energy should be a significant factor in the principal components, especially when comparing systems under different conditions.
5. Phase Transitions
Expectation: For systems with multiple energy levels, phase transitions (e.g., solid to liquid to gas) can be observed as temperature changes.
Observation: Sharp changes in thermodynamic properties (e.g., entropy, energy) should indicate phase transitions.
PCA Insight: PCA should reveal distinct clusters corresponding to different phases, with principal components capturing the transition points.
6. Energy-Level Shape and Chemical Potential
Expectation: The shape of the energy-level distribution and chemical potential should significantly affect the thermodynamic properties.
Observation: Systems with exponentially spaced energy levels might behave differently from those with linearly spaced levels.
PCA Insight: PCA should distinguish between different energy-level shapes and chemical potentials, highlighting their impact on the overall system behavior.
7. Inaccessible Micro-States
Expectation: The percentage of "inaccessible" micro-states (where ( P < 10^{-6} )) should decrease with increasing temperature.
Observation: At low temperatures, many high-energy states should be effectively inaccessible, while at high temperatures, more states become accessible.
PCA Insight: The percentage of inaccessible states should be a significant factor in the principal components, particularly for systems with a wide range of energy levels.
8. Linear Algebra and PCA Interpretation
Expectation: PCA should reduce the dimensionality of the feature matrix (species × thermodynamic properties) while retaining the most significant variations.
Observation: The first few principal components should capture the majority of the variance in the data, with each component corresponding to a physically interpretable combination of thermodynamic properties.
PCA Insight: The loading vectors (eigenvectors) should reveal which thermodynamic properties are most influential in distinguishing different species or conditions.
Educational Implications
These expected observations can be used to create educational modules that guide users through the underlying physics and mathematics. For example:

Temperature Effects: Show how increasing temperature affects the partition function, energy, and entropy, and how PCA captures these changes.
Phase Transitions: Demonstrate how PCA can identify phase transitions by clustering thermodynamic properties.
Energy-Level Shapes: Explore how different energy-level distributions affect the system's behavior and how PCA distinguishes these cases.
By aligning your app's outputs with these well-known scientific principles, you can provide a robust educational tool that not only simulates thermodynamic systems but also validates fundamental concepts in statistical mechanics and PCA. Would you like to explore how to incorporate these expectations into your app's user interface or documentation?
```
🛠️ Future Work (Roadmap)

Structured energy spectra
- harmonic oscillator, double-well, clustered levels
- compare random vs structured PCA signatures
Temperature / $\mu$-sweeps to build richer, higher-dimensional datasets


Advanced analytics
– kernel-PCA, UMAP, clustering
– information-theoretic measures (KL divergence between species)

Feel free to open issues or PRs for any of the above!

👤 Author & Licence

    L Shriver 
    Licensed under the MIT Licence. Contributions welcome.
