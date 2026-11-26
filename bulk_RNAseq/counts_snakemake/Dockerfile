FROM snakemake/snakemake:v8.14.0

RUN apt-get update
RUN dpkg --configure -a
RUN apt-get install -y unzip vim screen

ARG TARGETPLATFORM
ARG BUILDPLATFORM
ARG TARGETARCH

ARG R_VERSION=4.4.0

RUN echo "Running: $BUILDPLATFORM, building: $TARGETPLATFORM, arch: $TARGETARCH" > /docker_build.log

RUN if ["$TARGETPLATFORM" = "linux/arm64"]; then \
  curl https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip -o awscliv2.zip; \
  curl -O https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb; \
    else \
  curl https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip -o awscliv2.zip; \
  curl -O https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb; \
  fi

RUN unzip awscliv2.zip \
  && ./aws/install \
  && rm -rf aws awscliv2.zip

RUN apt-get -y install gdebi-core
RUN gdebi -n r-${R_VERSION}_1_amd64.deb \
  && ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R \
  && ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

RUN apt-get install -y libssl-dev libpng-dev
COPY scripts/ /scripts/
RUN Rscript /scripts/install.R 

COPY envs /envs
RUN mamba install -y bioconda::snakemake-storage-plugin-s3
RUN mamba env create --name salmon --file /envs/salmon.yaml
RUN mamba env create --name htseq --file /envs/htseq.yaml 

WORKDIR /workflow
COPY workflow/ /workflow/
COPY config/ /config/

COPY script.sh /script.sh
RUN chmod +x /script.sh

ENTRYPOINT ["/script.sh"]

    





