#!/bin/bash
source venv/bin/activate
export APP_ENV=remote
pip install -r requirements.txt
python3 code/app.py