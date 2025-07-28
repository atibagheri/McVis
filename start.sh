#!/bin/bash

# Start Plumber R API in background
Rscript -e "pr <- plumber::plumb('backend/degviz_api/plumber.R'); pr$run(host='0.0.0.0', port=8000)" &

# Start Flask app
cd backend
python3 app.py
