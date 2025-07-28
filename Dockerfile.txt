FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    r-base \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    git \
    curl \
    npm

# Install Plumber
RUN R -e "install.packages('plumber', repos='https://cloud.r-project.org/')"

# Create working directory
WORKDIR /app

# Copy entire repo
COPY . .

# Install Flask dependencies
RUN pip3 install -r backend/requirements.txt

# Build React app
WORKDIR /app/degviz
RUN npm install && npm run build

# Expose ports
EXPOSE 5050 8000

# Startup script
WORKDIR /app
RUN chmod +x start.sh
CMD ["./start.sh"]
