import numpy as np
from scipy.stats import gmean
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Constants
kB = 8.617E-5       # Boltzmann constant in eV/K

def clr_transform(X):
    """Centered log-ratio transformation for compositional data"""
    X_safe = X + np.finfo(float).eps   # protect against log(0)
    gm = gmean(X_safe, axis=1)
    return np.log(X_safe / gm.reshape(-1, 1))

def calculate_boltzmann_statistics(E, T, thresh_p=1e-6):
    """Calculate Boltzmann statistics for given energy matrix and 
    temperature (in Kelvins)"""
    beta = 1 / (kB * T)

    # Partition function (sum over energy levels for each species)
    Z = np.sum(np.exp(-beta * E), axis=1)

    # Normalized Boltzmann probabilities 
    P = np.exp(-beta * E) / Z.reshape(-1, 1)

    # Average energy for each species 
    E_avg = np.sum(P * E, axis=1)

    # Entropy (eV/K)
    P_safe = np.maximum(P, np.finfo(float).eps)     # protect against log(0)
    S = -kB * sum(P * np.log(P_safe), axis=1)

    # Helmholtz free energy (eV)
    F = -kB * T * np.log(Z)

    # Percentage of inaccessible states
    pct_inaccessible = 100 * np.sum(P < thresh_p, axis=1) / E.shape[1]

    return P, Z, E_avg, S, F, pct_inaccessible

def perform_pca_analysis(data, feature_names):
    """Perform PCA analysis and return results"""
    # Standardize the data
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    # Perform PCA
    pca = PCA()
    scores = pca.fit_transform(data_scaled)

    # Get loadings
    loadings = pca.components_.T
    
    # Explained variance
    explained_var = pca.explained_variance_ratio_*100

    return scores, loadings, explained_var, pca

def create_3d_scatter(scores, color_values, color_name, title):
    """Create 3D scatter plot with PCA scores"""

    scatter_colorscale = [
        [0.0, "#00E8FF"],
        [0.2, "#2FBDF7"],
        [0.4, "#5E91EE"],
        [0.6, "#8D66E6"],
        [0.8, "#BC3ADD"],
        [1.0, "#EB0FD5"]
    ]

    fig = go.Figure(data=go.Scatter3d(
        x=scores[:, 0],
        y=scores[:, 1],
        z=scores[:, 2],
        mode='markers',
        marker=dict(
            size=8,
            color=color_values,
            colorscale='scatter_colorscale',
            showscale=dict(title=color_name)
        ),
        text=[f'Species {i+1}' for i in range(len(scores))],
        hovertemplate='<b>%{text}</b><br>' + 
                      'PC1: %{x:.3f}<br>' +
                      'PC2: %{y:.3f}<br>' +
                      'PC3: %{z:.3f}<br>' +
                      f'{color_name}: %{{marker.color:.3f}}<extra></extra>'
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='PC1',
            yaxis_title='PC2',
            zaxis_title='PC3',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
        ),
        height=600
    )

    return fig

def create_correlation_heatmap(scores, original_data, var_names, title):
    """Create correlation heatmap between PCA scores and orignal variables"""
    n_pcs = min(5, scores.shape[1])     # show up to 5 PCs
    n_vars = original_data.shape[1]

    correlation_colorscale = [
        [0.0, "#ccff33"],
        [1/7, "#9ef01a"],
        [2/7, "#70e000"],
        [3/7, "#38b000"],
        [4/7, "#008000"],
        [5/7, "#007200"],
        [6/7, "#006400"],
        [1.0, "#004b23"]
    ]

    corr_matrix = np.zeros((n_vars, n_pcs))
    for i in range(n_vars):
        for j in range(n_pcs):
            corr_matrix[i, j] = np.corrcoef(scores[:, j], original_data[:, i])[0, 1]

    fig = go.Figure(data=go.Heatmap(
        z=corr_matrix,
        x=[f'PC{i+1}' for i in range(n_pcs)],
        y=var_names,
        colorscale=correlation_colorscale,
        zmid=0,
        colorbar=dict(title='Correlation')
    ))

    fig.update_layout(
        title=title,
        xaxis_title='Principal Components',
        yaxis_title='Variables',
        height=400
    )

    return fig

def create_variance_plot(explained_var, title):
    """Create explained variance plot"""
    fig = go.Figure(data=go.Bar(
        x=[f'PC{i+1}' for i in range(min(10, len(exlained_var)))],
        y=explained_var[:10],
        marker_color='#ffea00'
    ))

    fig.update_layout(
            title=title,
            xaxis_title='Principal Componenet',
            yaxis_title='% Variance Explained',
            height=400
        )
    
    return fig