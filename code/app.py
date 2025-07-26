import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from colormaps import *
import style
import os

from analysis_utils import (
    clr_transform,
    calculate_boltzmann_statistics,
    perform_pca_analysis,
    create_3d_scatter,
    create_correlation_heatmap,
    create_variance_plot
)

env = os.getenv('APP_ENV', 'remote')     # default to local
if env == 'remote':
    background_path = 'code/static/images/wisp.jpg'
elif env == 'local':
    background_path = 'pca_projects/thermo_pca/code/static/images/wisp.jpg'
else:
    raise ValueError(f"Unknown environment: {env}")

# load your CSS & background
style.load_custom_css()
style.apply_background(background_path)

st.set_page_config(
    page_title="Thermodynamics & PCA Analysis",
    page_icon="code/static/images/ember.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    st.title("Thermodynamics & PCA Analysis")
    st.markdown("Statistical mechanics simulation with Boltzmann statistics and compositional data analysis")

    # ─── Sidebar inputs ────────────────────────────────────────
    seed         = st.sidebar.number_input("Random Seed", 0, 10000, 27)
    np.random.seed(seed)
    n_species    = st.sidebar.slider("Number of species", 5, 50, 10)
    n_levels     = st.sidebar.slider("Energy Levels per Species", 3, 10, 5)
    T            = st.sidebar.selectbox("Temperature (K)", [100,200,300,400,500], index=2)
    thresh_p     = st.sidebar.number_input(
        "Inaccessible Probability Threshold",
        min_value=1e-8, max_value=1e-3, value=1e-6,
        format="%.2e"
    )
    color_choice = st.sidebar.selectbox(
        "3D Plot Color Variable",
        ['AvgE','Entropy','Z','F','PctInaccess'], index=1
    )

    # ─── Generate data & Boltzmann stats ───────────────────────
    E = np.random.rand(n_species, n_levels)
    P, Z, mean_E, S, F, pct_inacc = calculate_boltzmann_statistics(E, T, thresh_p)

    # ─── PCA on probabilities ──────────────────────────────────
    P_clr      = clr_transform(P)
    scores_P, load_P, ev_P, _ = perform_pca_analysis(
        P_clr, [f'E{i+1}' for i in range(n_levels)]
    )

    # ─── PCA on thermodynamic features ────────────────────────
    X_thermo      = np.column_stack([mean_E, S, Z, F, pct_inacc])
    feature_names = ['<E>','S','Z','F','%Inacc']
    scores_X, load_X, ev_X, _ = perform_pca_analysis(X_thermo, feature_names)

    # ─── Correlation heatmaps ─────────────────────────────────
    fig_corr_P = create_correlation_heatmap(
        scores_P, P_clr, [f'E{i+1}' for i in range(n_levels)],
        "Correlations: CLR-Probabilities vs PCA Scores"
    )
    fig_corr_X = create_correlation_heatmap(
        scores_X, X_thermo, feature_names,
        "Correlations: Thermo Features vs PCA Scores"
    )

    # ─── Cross-correlation matrix ─────────────────────────────
    n_pcs = min(scores_P.shape[1], scores_X.shape[1], 5)
    cross_corr = np.zeros((n_pcs, n_pcs))
    for i in range(n_pcs):
        for j in range(n_pcs):
            cross_corr[i,j] = np.corrcoef(scores_P[:,i], scores_X[:,j])[0,1]
    fig_cross = go.Figure(go.Heatmap(
        z=cross_corr,
        x=[f"Thermo PC{i+1}" for i in range(n_pcs)],
        y=[f"Prob PC{i+1}"   for i in range(n_pcs)],
        colorscale=GREEN, zmid=0,
        colorbar=dict(title="Correlation")
    ))
    fig_cross.update_layout(height=300)

    # ─── Three columns: 1∶2∶1 ───────────────────────────────────
    col1, col2, col3 = st.columns([1, 2, 1])

    # ─ COLUMN 1: Energy + Thermo Metrics ─────────────────────
    with col1:
        st.subheader("Energy Landscape (eV)")
        figE = go.Figure(go.Heatmap(
            z=E,
            x=[f"Level {i+1}"   for i in range(n_levels)],
            y=[f"Species {i+1}" for i in range(n_species)],
            colorscale=BLUE,
            colorbar=dict(title="E (eV)")
        ))
        figE.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=40, r=40, t=40, b=40),
            height=300
        )
        st.plotly_chart(figE, use_container_width=True)

        st.subheader("Thermodynamic Properties")
        st.metric("Temperature",       f"{T} K")
        st.metric("Avg Entropy",       f"{np.mean(S):.6f} eV/K")
        st.metric("Avg Free Energy",   f"{np.mean(F):.3f} eV")
        st.metric("Pct Inaccessible",  f"{np.mean(pct_inacc)*100:.2f}%")

    # ─ COLUMN 2: PCA View + Variance ──────────────────────────
    with col2:
        view = st.radio(
            "Select PCA View",
            ("Boltzmann Probabilities", "Thermodynamic Features"),
            horizontal=True
        )

        if view == "Boltzmann Probabilities":
            st.subheader("PCA: Boltzmann Distributions")
            fig3d = create_3d_scatter(
                scores_P,
                {'AvgE': mean_E,'Entropy': S,'Z': Z,'F': F,'PctInaccess': pct_inacc}[color_choice],
                color_choice,
                "CLR-transformed",
                vector_loadings=load_P[:5],
                vector_labels=[f"P{i+1}" for i in range(5)]
            )
            st.plotly_chart(fig3d, use_container_width=True)

            st.subheader("Explained Variance – Probabilities")
            figVar = create_variance_plot(ev_P, "Probabilities")
            st.plotly_chart(figVar, use_container_width=True)

        else:
            st.subheader("PCA: Thermodynamic Features")
            fig3d = create_3d_scatter(
                scores_X,
                {'AvgE': mean_E,'Entropy': S,'Z': Z,'F': F,'PctInaccess': pct_inacc}[color_choice],
                color_choice,
                "Thermo features",
                vector_loadings=load_X,
                vector_labels=feature_names
            )
            st.plotly_chart(fig3d, use_container_width=True)

            st.subheader("Explained Variance – Thermo Features")
            figVar = create_variance_plot(ev_X, "Thermodynamic Features")
            st.plotly_chart(figVar, use_container_width=True)

    # ─ COLUMN 3: Correlation Heatmaps ─────────────────────────
    with col3:
        st.subheader("Correlations: Probabilities")
        st.plotly_chart(fig_corr_P, use_container_width=True)

        st.subheader("Correlations: Thermo Features")
        st.plotly_chart(fig_corr_X, use_container_width=True)

        st.subheader("Cross-Correlation of PCA Spaces")
        st.plotly_chart(fig_cross, use_container_width=True)

if __name__ == "__main__":
    main()