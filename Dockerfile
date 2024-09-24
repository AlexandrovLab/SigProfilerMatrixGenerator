# Use an official Python runtime as a parent image
FROM python:3.8-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install SigProfilerMatrixGenerator from PyPI
RUN pip install SigProfilerMatrixGenerator==1.2.29

# Create a non-root user named 'spm_user'
RUN useradd -m -s /bin/bash spm_user

# Change the ownership of the /usr/src/app directory and its contents to the new non-root user
RUN chown -R spm_user:spm_user /usr/src/app

# Switch to the non-root user for subsequent commands and when running the container
USER spm_user
