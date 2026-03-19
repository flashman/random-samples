# Use the official Python image as a base
FROM python:3.11-slim

# Set environment variables
ENV POETRY_VERSION=1.8.3
ENV POETRY_VIRTUALENVS_CREATE=false

# Install Poetry
RUN pip install "poetry==$POETRY_VERSION"

# Set the working directory
WORKDIR /app

# Copy pyproject.toml and poetry.lock files
COPY pyproject.toml poetry.lock* /app/

# Install the dependencies
RUN poetry install --no-root --only main

# Copy the rest of the application code
COPY . /app

# Set the entrypoint to the CLI command
ENTRYPOINT ["python", "workflow_runner/main.py"]
