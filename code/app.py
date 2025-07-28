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
    raise ValueError(f"Unknown enviornment: {env}")

style.load_custom_css()
style.apply_background(background_path)

# Streamlit page configuration
st.set_page_config(
    page_title="Thermodynamics & PCA Analysis",
    page_icon="code/static/images/ember.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    st.title("Thermodynamics & PCA Analysis")
    st.markdown("Statistical mechanics simulation with Boltzmann statistics and compositional data analysis")

    # Sidebar
    seed = st.sidebar.number_input("Random Seed", 0, 10000, 27)
    np.random.seed(seed)
    n_species = st.sidebar.slider("Number of species", 5, 50, 10)
    n_levels = st.sidebar.slider("Energy Levels per Species", 3, 10, 5)
    T = st.sidebar.selectbox("Temperature (K)", [100,200,300,400,500], index=2)
    thresh_p = st.sidebar.number_input("Inaccessible Probability Threshold", 
                                       min_value=1e-8, max_value=1e-3, value=1e-6,
                                       format="%.2e")
    color_choice = st.sidebar.selectbox(
        "3D Plot Color Variable",
        ['AvgE', 'Entropy', 'Z', 'F', 'PctInaccess'], index=1
    )

    # Data generation & Boltzmann
    E = np.random.rand(n_species, n_levels)
    P, Z, mean_E, S, F, pct_inacc = calculate_boltzmann_statistics(E, T, thresh_p)

    # Tabs
    tabs = st.tabs([
        "Data Overview", "Boltzmann Analysis", 
        "PCA - Probabilities", "PCA - State Variables", "Correlations"
    ])

    # -- Tab 1: Overview --
    with tabs[0]:
        st.header("System Overview")
        c1, c2 = st.columns([3, 1])
        # Energy table & heatmap
        with c1:
            st.subheader("Energy Landscape (eV)")
            dfE = pd.DataFrame(E, columns=[f'E{i+1}' for i in range(n_levels)])
            dfE.index = [f'Species {i+1}' for i in range (n_species)]
            #st.dataframe(dfE, use_container_width=True)
            figE = go.Figure(data=go.Heatmap(
                z=E, x=[f'Level {i+1}' for i in range(n_levels)],
                y=[f'Species {i+1}' for i in range(n_species)], colorscale=BLUE,
                colorbar=dict(title='Energy (eV)')
            ))
            figE.update_layout(height=400)
            st.plotly_chart(figE, use_container_width=True)
        # Thermo features & metrics
        with c2:
            st.subheader("Thermodynamic Properties")
            thermo_df = pd.DataFrame({
                'Species': [f'S{i+1}' for i in range(n_species)],
                'AvgE (eV)': mean_E,
                'Entropy (eV/K)': S,
                'Z': Z,
                'F (eV)': F,
                'PctInacces (%)': pct_inacc
            })
            #st.dataframe(thermo_df, use_container_width=True)

            # System statistics
            st.metric("Temperature", f"{T} K")
            st.metric("Average Entropy", f"{np.mean(S):.6f} eV/K")
            st.metric("Average Free Energy", f"{np.mean(F):.3f} eV")

    # -- Tab 2: Boltzmann Analysis --
    with tabs[1]:
        st.markdown("<h1 class='gradient_text1'> Boltzmann Probability Analysis </h1>", unsafe_allow_html=True)
        fig_prob = go.Figure(data=go.Heatmap(
            z=P,
            x=[f'Level {i+1}' for i in range(n_levels)],
            y=[f'Species {i+1}' for i in range(n_species)],
            colorscale=GREEN,
            colorbar=dict(title='Probability')
        ))
        fig_prob.update_layout(
            title="Boltzmann Probabillity Matrix",
            xaxis_title="Energy Level Index",
            yaxis_title="Species Index",
            height=500
        )
        st.plotly_chart(fig_prob, use_container_width=True)

        # Statistics summary
        col1, col2, col3 = st.columns(3)
        with col1: 
            st.metric("Max Probability", f"{np.max(P):.4f}")
        with col2:
            st.metric("Min Probability", f"{np.min(P):.2e}")
        with col3:
            st.metric("States < Threshold", f"{np.sum(P < thresh_p)} / {P.size}")

    # -- Tab 3: PCA on Boltzmann Probabilities -- 
    with tabs[2]:
        st.header("PCA Analysis - Probability Distributions")

        # Apply CLR transformation to probability matrix
        P_clr = clr_transform(P)

        # Perform PCA on CLR-transformed probabilities
        scores_P, loadings_P, explained_var_P, pca_P = perform_pca_analysis(
            P_clr, [f'E{i+1}' for i in range(n_levels)]
        )

        col1, col2 = st.columns(2)

        with col1: 
            # 3D scatter plot
            color_map = {
                'AvgE': mean_E,
                'Entropy': S,
                'Z': Z,
                'F': F,
                'PctInaccess': pct_inacc
            }
            color_values = color_map[color_choice]

            fig_3d_P = create_3d_scatter(scores_P, color_values, color_choice,
                                         "PCA of Boltzmann Distributions (CLR-transformed)",
                                         vector_loadings=loadings_P[:5],
                                         vector_labels=[f'P{i+1}' for i in range(5)])
            
            st.plotly_chart(fig_3d_P, use_container_width=True)

            # Loadings table
            st.subheader("Loadings (PC1-PC3)")
            loadings_df = pd.DataFrame(
                data=loadings_P[:, :3],
                columns=['PC1', 'PC2', 'PC3'],
                index=[f'E{i+1}' for i in range(n_levels)]
            )
            st.dataframe(loadings_df, use_container_width=True)

        with col2:
            # Explained variance plot
            fig_var_P = create_variance_plot(explained_var_P, "Explained Variance - Probabilities")
            st.plotly_chart(fig_var_P, use_container_width=True)

            # Cumulative variance
            cumulative_var = np.cumsum(explained_var_P)
            st.subheader("Cumulative Variance Explained")
            for i in range(min(5, len(cumulative_var))):
                st.write(f"PC1-PC{i+1}: {cumulative_var[i]:.1f}%")

    # -- Tab 4: PCA on Thermodynamic Features --
    with tabs[3]:
        st.header("PCA Analysis - Thermodynamic Features")

        # Combine thermodynamic features
        X_thermo = np.column_stack([mean_E, S, Z, F, pct_inacc])
        feature_names = ['<E>', 'S', 'Z', 'F', '%Inaccess']

        # Perform PCA on thermodynamic features
        scores_X, loadings_X, explained_var_X, pca_X = perform_pca_analysis(X_thermo, feature_names)

        col1, col2 = st.columns(2)

        with col1:
            # 3D scatter plot
            color_map = {
                'AvgE': mean_E,
                'Entropy': S,
                'Z': Z,
                'F': F,
                'PctInaccess': pct_inacc
            }
            color_values = color_map[color_choice]

            fig_3d_X = create_3d_scatter(scores_X, color_values, color_choice,
                                         "PCA of Thermodynamic Features",
                                         vector_loadings=loadings_X,
                                         vector_labels=['E_avg', 'S', 'Z', 'F', 'Pct_Inacc'])
            st.plotly_chart(fig_3d_X, use_container_width=True)

            # Loadings table
            st.subheader("Loadings (PC1-PC3)")
            loadings_df = pd.DataFrame(
                data=loadings_X[:, :3],
                columns=['PC1', 'PC2', 'PC3'],
                index=feature_names
            )
            st.dataframe(loadings_df, use_container_width=True)

        with col2:
            # Explained variance plot
            fig_var_X = create_variance_plot(explained_var_X,
                                             "Explained Variance - Thermodynamic Feautres")
            st.plotly_chart(fig_var_X, use_container_width=True)

            # Cumulative Variance
            cumulative_var = np.cumsum(explained_var_X)
            st.subheader("Cumulative Variance Explained")
            for i in range(min(5, len(cumulative_var))):
                st.write(f"PC1-PC{i+1}: {cumulative_var[i]:.1f}%")

    # -- Tab 5: Correlations --
    with tabs[4]:
        st.header("Correlation Analysis")
        
        # Correlation heatmaps
        fig_corr_P = create_correlation_heatmap(
            scores_P, P_clr, [f'E{i+1}' for i in range(n_levels)],
            "Correlations: CLR-transformed Probabilities vs. PCA Scores"
        )
        st.plotly_chart(fig_corr_P, use_container_width=True)

        fig_corr_X = create_correlation_heatmap(
            scores_X, X_thermo, feature_names,
            "Correlations: Thermodynamics Features vs. PCA Scores"
        )
        st.plotly_chart(fig_corr_X, use_container_width=True)

        # Cross-correlation between the two PCA spaces
        st.subheader("Cross-correlation: Probablity PCA vs. Thermo PCA")
        n_pcs = min(5, scores_P.shape[1], scores_X.shape[1])
        cross_corr = np.zeros((n_pcs, n_pcs))
        for i in range(n_pcs):
            for j in range(n_pcs):
                cross_corr[i, j] = np.corrcoef(scores_P[:, i], scores_X[:, j])[0, 1]

        fig_cross = go.Figure(data=go.Heatmap(
            z=cross_corr,
            x=[f'Thermo PC{i+1}' for i in range(n_pcs)],
            y=[f'Prob PC{i+1}' for i in range(n_pcs)],
            colorscale=GREEN,
            zmid=0,
            colorbar=dict(title='Correlation')
        ))
        fig_cross.update_layout(
            title="Cross-correlation between PCA Spaces",
            showlegend=False,
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            margin=dict(l=40, r=40, t=40, b=40),
            height=400
        )
        st.plotly_chart(fig_cross, use_container_width=True)

    # Footer with information
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
