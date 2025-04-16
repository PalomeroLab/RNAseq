FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Create local bin directory
RUN mkdir -p /home/rstudio/.local/bin

# Install Micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
        tar -xvj -C /home/rstudio/.local/bin --strip-components=1 bin/micromamba

# Install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
        unzip awscliv2.zip && \
        ./aws/install -i /home/rstudio/.local/aws-cli -b /home/rstudio/.local/bin && \
        rm -rf aws awscliv2.zip

# Add local bin to PATH
ENV PATH="/home/rstudio/.local/bin:${PATH}"
