import streamlit as st
import base64
import os
from pathlib import Path

env = os.getenv('APP_ENV', 'remote')     # default to local
if env == 'remote':
    BASE_DIR = Path('.')
elif env == 'local':
    BASE_DIR = Path('/home/lshriver/projects/myPortfolio/pca_projects/thermo_pca/')
else:
    raise ValueError(f"Unknown enviornment: {env}")

STREAMLIT_CSS = BASE_DIR / "code/static/style.css"
BG_IMAGE = BASE_DIR / "code/static/images/wisp.jpg"

def load_custom_css():
    if STREAMLIT_CSS.exists():
        css = STREAMLIT_CSS.read_text()
        st.markdown(f"<style>{css}</style>", unsafe_allow_html=True)

def get_img_as_base64(fp: str) -> str:
    with open(fp, "rb") as f:
        return base64.b64encode(f.read()).decode()
    
def apply_background(image_path: str):
    if os.path.exists(image_path):
        img_b64 = get_img_as_base64(image_path)
        st.markdown(f"""
            <style>
            .stApp {{
                background-image: url("data:image/jpeg;base64,{img_b64}");
                background-size: cover;
                background-repeat: no-repeat;
                background-attachment: fixed;
                background-position: center;
            }}
            </style>
        """, unsafe_allow_html=True)