import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import gmean
from matplotlib.colors import LinearSegmentedColormap
import os

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
    