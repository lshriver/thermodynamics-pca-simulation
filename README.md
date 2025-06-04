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

## 🚀 Quick Start

```bash
git clone https://github.com/your-handle/thermo-pca.git
cd thermo-pca
sd
```

# Launch MATLAB and run the main demo

```
>> thermo_pca_demo        % or `run("thermo_pca_demo.m")`
```
The script will:
1. Print two summary tables to the Command Window;
2. Open six publication-ready figures (set the CI env-var to suppress GUIs).


📁 Repository Layout
```
.
├── code/            % MATLAB scripts & functions
│   └── thermo_pca_demo.m
├── docs/            % exported figures, SVGs, and extra reading
├── notes/           % derivations, scratch, future-work ideas
└── README.md
```
Python assets will live in py/ once the port begins.

🛠️ Future Work (Roadmap)

Structured energy spectra
- harmonic oscillator, double-well, clustered levels
- compare random vs structured PCA signatures
Temperature / $\mu$-sweeps to build richer, higher-dimensional datasets

Full Python rewrite
– NumPy/SciPy + pandas workflow
– interactive Plotly / Streamlit dashboard for portfolio demos
– package on PyPI (pip install thermo-pca)

Advanced analytics
– kernel-PCA, UMAP, clustering
– information-theoretic measures (KL divergence between species)

Feel free to open issues or PRs for any of the above!

👤 Author & Licence

    L Shriver 
    Licensed under the MIT Licence. Contributions welcome.
