import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go 

from analysis_utils import (
    clr_transform,
    calculate_boltzmann_statistics,
    perform_pca_analysis,
    create_3d_scatter,
    create_correlation_heatmap,
    create_variance_plot
)

# Streamlit page configuration
st.set_page_config(
    page_title="Thermodynamics & PCA Analysis",
    page_icon=":fire:",
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    st.title("Thermodynamics & PCA Analysis")
    st.markdown("*Statistical mechanics simulation with Boltzmann statistics and compositional data analysis")

    # Sidebar
    seed = st.sidebar.number_input("Random Seed", 0, 10000, 27)
    np.random.seed(seed)
    n_species = st.sidebar.slider("Number of species", 5, 50, 10)
    n_levels = st.sidebar.slider("Energy Levels per Species", 3, 5, 10)
    T = st.sidebar.selectbox("Temperature (K)", [100,200,300,400,500], index=2)
    thresh_p = st.sidebar.number_input("Inaccessible Probability Threshold", 1e-8, 1e-6, 1e-3, forma="%.2e")
    color_choice = st.sidebar.selectbox(
        "3D Plot Color Variable",
        ['AvgE', 'Entropy', 'Z', 'F', 'PctInaccess'], index=1
    )

    # Data generation & Boltzmann
    E = np.random.rand(n_species, n_levels)
    P, Z, mean_E, S, F, pct_inacc = calculate_boltzmann_statistics(E, T, thres_p)

    # Tabs
    tabs = st.tabs([
        "Data Overview", "Boltzmann Analysis", 
        "PCA - Probabilities", "PCA - State Variables", "Correlations"
    ])

    # -- Tab 1: Overview --
    with tabs[0]:
        st.header("System Overview")
        c1, c2 = st.columns(2)
        # Energy table & heatmap
        with c1:
            st.subheader("Energy Landscape (eV)")

        # Thermo features & metrics
        with c2:
            st.subheader("Thermodynamic Properties")

# -- Tab 2: Boltzmann Analysis --
with tabs[1]:
    st.header("Boltzmann Probability Analysis")

# -- Tab 3: PCA on Boltzmann Probabilities -- 
    st.header("PCA Analysis - Probability Distributions")

# -- Tab 4: PCA on Thermodynamic Features --
with tabs[3]:
    st.healder("PCA Analysis - Thermodynamic Features")

# -- Tab 5: COrrelations --
with tabs[4]:
    st.header("Correlation Analysis")


    st.subheader("Cross-correlation: Probablity PCA vs. Thermo PCA")

st.markdown("---")
st.markdown("""
            **About this analysis:**
            - Random E landscapes (0-1 eV)
            - Boltzmann stats at selected T (k_B=8.617e-5 eV/K)
            - CLR transform $\rightarrow$ PCA on both distributions and thermo features
            - Interactive 3D plots and correlation heatmaps
            """)

if __name__ == "__main__":
    main()