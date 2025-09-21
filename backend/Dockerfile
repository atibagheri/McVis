# backend/Dockerfile
FROM python:3.12-slim

ENV VIRTUAL_ENV=/opt/venv \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1
RUN python -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# --- NEW: build tools & dev headers for wordcloud/lxml/pillow ---
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gcc g++ python3-dev \
      libxml2-dev libxslt1-dev \
      zlib1g-dev libjpeg-dev libpng-dev libfreetype6-dev \
      curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first (better layer cache)
COPY requirements.txt /app/
RUN pip install -r requirements.txt

# Copy the rest
COPY . /app

ENV FLASK_ENV=production \
    FLASK_PORT=5050 \
    R_SERVICE_URL=http://rapi:8000 \
    MOUSE_TPM_PATH=/app/degviz_api/data/MedianNorm_TPMs_order_Mouse.txt

EXPOSE 5050
CMD ["python", "app.py", "--host=0.0.0.0", "--port=5050"]
