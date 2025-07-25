import streamlit as st
import base64
import os

st.write("Current working directory:", os.getcwd())
st.write("File exists:", os.path.exists("code/static/images/wisp.jpg"))

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

apply_background("code/static/images/wisp.jpg")
st.write("Test background")