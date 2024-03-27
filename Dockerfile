# How to build and deploy the Docker image:
# docker build --build-arg VER=1.2.3 --no-cache -t ecogenomic/gtdbtk:latest -t ecogenomic/gtdbtk:1.2.3 .
# docker run -v /host/gtdbtk_io:/data -v /host/release_data:/refdata ecogenomic/gtdbtk classify_wf --genome_dir /data/genomes --out_dir /data/output
# docker push ecogenomic/gtdbtk:latest && sudo docker push ecogenomic/gtdbtk:1.2.3

FROM python:3.8-slim-bullseye

ARG VER

# ---------------------------------------------------------------------------- #
# --------------------- INSTALL HMMER, PYTHON3, FASTTREE---------------------- #
# ---------------------------------------------------------------------------- #
RUN apt-get update -y -m && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --no-install-suggests -y \
        wget \
        libgomp1 \
        libgsl25 \
        libgslcblas0 \
        build-essential \
        curl \
        hmmer=3.* \
        mash=2.2.* \
        prodigal=1:2.6.* \
        fasttree=2.1.* \
        unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/fasttreeMP /usr/bin/FastTreeMP && \
    curl https://sh.rustup.rs -sSf | sh -s -- -y

ENV PATH="/root/.cargo/bin:${PATH}"

# ---------------------------------------------------------------------------- #
# ----------------------------- INSTALL PPLACER ------------------------------ #
# ---------------------------------------------------------------------------- #
RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip -q && \
    unzip pplacer-linux-v1.1.alpha19.zip && \
    mv pplacer-Linux-v1.1.alpha19/* /usr/bin && \
    rm pplacer-linux-v1.1.alpha19.zip && \
    rm -rf pplacer-Linux-v1.1.alpha19

# ---------------------------------------------------------------------------- #
# ------------------------------ INSTALL SKANI ------------------------------- #
# ---------------------------------------------------------------------------- #

RUN wget https://github.com/bluenote-1577/skani/archive/refs/tags/v0.2.1.tar.gz
RUN tar -xvf v0.2.1.tar.gz
RUN cd skani-0.2.1 && cargo install --path . --root /usr
RUN chmod +x /usr/bin/skani
RUN cd ../
RUN rm -rf v0.2.1.tar.gz skani-0.2.1

# ---------------------------------------------------------------------------- #
# --------------------- SET GTDB-TK MOUNTED DIRECTORIES ---------------------- #
# ---------------------------------------------------------------------------- #
RUN mkdir /refdata && \
    mkdir /data
ENV GTDBTK_DATA_PATH="/refdata/"

# ---------------------------------------------------------------------------- #
# --------------------------- INSTALL PIP PACKAGES --------------------------- #
# ---------------------------------------------------------------------------- #
RUN python -m pip install --upgrade pip && \
    python -m pip install gtdbtk==${VER} \

# ---------------------------------------------------------------------------- #
# ---------------------------- SET THE ENTRYPOINT ---------------------------- #
# ---------------------------------------------------------------------------- #
ENTRYPOINT ["gtdbtk"]
CMD ["--help"]
