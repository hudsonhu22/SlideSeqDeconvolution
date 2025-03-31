# Start from the Bioconductor base image
FROM bioconductor/bioconductor_docker:RELEASE_3_16

# Install STdeconvolve
RUN R -e 'BiocManager::install("STdeconvolve", update=FALSE, ask=FALSE)'

# Verify installation
RUN R -e 'library(STdeconvolve); cat("STdeconvolve successfully installed!\n")'

# Create a demo script
RUN echo 'library(STdeconvolve); data(mOB); cat("mOB dataset loaded with dimensions:", dim(mOB$counts)[1], "spots and", dim(mOB$counts)[2], "genes\n")' > /usr/local/bin/test_stdeconvolve.R

# Set the working directory
WORKDIR /data

# Command to run when the container starts
CMD ["R"]