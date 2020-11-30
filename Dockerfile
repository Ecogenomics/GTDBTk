# How to build and deploy the Docker image:
# sudo docker build --no-cache -t ecogenomic/gtdbtk:latest -t ecogenomic/gtdbtk:1.3.0 .
# sudo docker push ecogenomic/gtdbtk:latest && sudo docker push ecogenomic/gtdbtk:1.3.0

FROM ubuntu:18.04

# ---------------------------------------------------------------------------- #
# --------------------- INSTALL HMMER, PYTHON3, FASTTREE---------------------- #
# ---------------------------------------------------------------------------- #
RUN apt-get update -y -m && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
                                               wget \
                                               libgomp1 \
                                               hmmer=3.1b2+dfsg-5ubuntu1 \
                                               mash=2.0-2 \
                                               prodigal=1:2.6.3-1 \
                                               fasttree=2.1.10-1 \
                                               unzip \
                                               python3.8 \
                                               python3-pip \
                                               python3-setuptools \
                                               python3-wheel && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------- #
# ------------------------------ ALIAS PYTHON3  ------------------------------ #
# ---------------------------------------------------------------------------- #
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.8 1 && \
    update-alternatives --set python /usr/bin/python3.8

# ---------------------------------------------------------------------------- #
# ----------------------------- INSTALL PPLACER ------------------------------ #
# ---------------------------------------------------------------------------- #
RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip -q && \
    unzip pplacer-linux-v1.1.alpha19.zip && \
    mv pplacer-Linux-v1.1.alpha19/* /usr/bin && \
    rm pplacer-linux-v1.1.alpha19.zip && \
    rm -rf pplacer-Linux-v1.1.alpha19

# ---------------------------------------------------------------------------- #
# ----------------------------- INSTALL FASTANI ------------------------------ #
# ---------------------------------------------------------------------------- #
RUN wget https://github.com/ParBLiSS/FastANI/releases/download/v1.32/fastANI-Linux64-v1.32.zip -q && \
    unzip fastANI-Linux64-v1.32.zip -d /usr/bin && \
    rm fastANI-Linux64-v1.32.zip

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
    python -m pip install dendropy>=4.1.0 && \
    python -m pip install numpy>=1.9.0 && \
    python -m pip install tqdm>=4.31.0 && \
    python -m pip install gtdbtk==1.4.0

# ---------------------------------------------------------------------------- #
# ---------------------------- SET THE ENTRYPOINT ---------------------------- #
# ---------------------------------------------------------------------------- #
ENTRYPOINT ["gtdbtk"]
CMD ["--help"]
