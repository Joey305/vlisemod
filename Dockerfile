# Use Miniconda3 base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Copy environment file
COPY environment.yml .

# Create environment
RUN conda env create -f environment.yml

# Activate environment and make it default
RUN echo "conda activate base" >> ~/.bashrc

# Copy your project files
COPY . .

# Default command
# CMD ["python", "COPYapp.py"]
