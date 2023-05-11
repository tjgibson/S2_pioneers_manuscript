FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="6b22ac7a4fa77efecf74f70ac59d3d572cc347e25d31b86fa3e0009aa8180908"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/ccff920ef2dca61d7d93dbd598aef221
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bowtie2 ==2.4.4  # Keep consistent with version specified in bowtie2/build
#     - samtools ==1.10
RUN mkdir -p /conda-envs/ccff920ef2dca61d7d93dbd598aef221
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/bowtie2/align/environment.yaml /conda-envs/ccff920ef2dca61d7d93dbd598aef221/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/16114a9972039cb62b9656ba636b27b0
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bowtie2 ==2.4.4  # Keep consistent with version specified in bowtie2/align
#     - samtools ==1.10
RUN mkdir -p /conda-envs/16114a9972039cb62b9656ba636b27b0
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/bowtie2/build/environment.yaml /conda-envs/16114a9972039cb62b9656ba636b27b0/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/macs2/callpeak/environment.yaml
#   prefix: /conda-envs/c2a8201b8d46efaba69ff3379bedbb7b
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - macs2>=2.2
RUN mkdir -p /conda-envs/c2a8201b8d46efaba69ff3379bedbb7b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/macs2/callpeak/environment.yaml /conda-envs/c2a8201b8d46efaba69ff3379bedbb7b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/sambamba/view/environment.yaml
#   prefix: /conda-envs/9a5b04e8a8a154f32771f0a56dc701bf
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - sambamba ==0.8.0
RUN mkdir -p /conda-envs/9a5b04e8a8a154f32771f0a56dc701bf
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/sambamba/view/environment.yaml /conda-envs/9a5b04e8a8a154f32771f0a56dc701bf/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/samtools/idxstats/environment.yaml
#   prefix: /conda-envs/9608721699f97513ba7f47bd4e3db24b
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - samtools ==1.10
RUN mkdir -p /conda-envs/9608721699f97513ba7f47bd4e3db24b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.1.0/bio/samtools/idxstats/environment.yaml /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.3.2/bio/hisat2/align/environment.yaml
#   prefix: /conda-envs/b33f755bffdf8829e59d1023684a2b44
#   channels:
#     - bioconda
#     - defaults
#   dependencies:
#     - hisat2 ==2.1.0
#     - samtools ==1.9
RUN mkdir -p /conda-envs/b33f755bffdf8829e59d1023684a2b44
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.3.2/bio/hisat2/align/environment.yaml /conda-envs/b33f755bffdf8829e59d1023684a2b44/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.3.2/bio/subread/featurecounts/environment.yaml
#   prefix: /conda-envs/d4a8e62782e6d0fd8b849811a698b9d8
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - subread =2.0
RUN mkdir -p /conda-envs/d4a8e62782e6d0fd8b849811a698b9d8
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.3.2/bio/subread/featurecounts/environment.yaml /conda-envs/d4a8e62782e6d0fd8b849811a698b9d8/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/DEseq2.yaml
#   prefix: /conda-envs/a683b7068d9e2ce45ba985ef66cc1c89
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - r-tidyverse = 1.3.1
#     - bioconductor-deseq2
#     - bioconductor-rtracklayer = 1.50.0
RUN mkdir -p /conda-envs/a683b7068d9e2ce45ba985ef66cc1c89
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/DEseq2.yaml /conda-envs/a683b7068d9e2ce45ba985ef66cc1c89/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/NGmerge.yaml
#   prefix: /conda-envs/6f7502c45f8f32d3c16f398e1e3cfbff
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - ngmerge = 0.3
RUN mkdir -p /conda-envs/6f7502c45f8f32d3c16f398e1e3cfbff
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/NGmerge.yaml /conda-envs/6f7502c45f8f32d3c16f398e1e3cfbff/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/curl.yaml
#   prefix: /conda-envs/165a6884838c8278606bcee5ee86ab20
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - curl=7.65.3
#     - pigz=2.6
RUN mkdir -p /conda-envs/165a6884838c8278606bcee5ee86ab20
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/curl.yaml /conda-envs/165a6884838c8278606bcee5ee86ab20/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/deeptools.yaml
#   prefix: /conda-envs/71fb0136399bd6f669cff3c445df6818
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - deeptools = 3.5.1
RUN mkdir -p /conda-envs/71fb0136399bd6f669cff3c445df6818
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/deeptools.yaml /conda-envs/71fb0136399bd6f669cff3c445df6818/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/samtools.yaml
#   prefix: /conda-envs/da0a0a3f54e418d0bf464e31b190514c
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - samtools=1.14
RUN mkdir -p /conda-envs/da0a0a3f54e418d0bf464e31b190514c
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/samtools.yaml /conda-envs/da0a0a3f54e418d0bf464e31b190514c/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/seqkit.yaml
#   prefix: /conda-envs/b9f098cecf23a67fb2bdcd5108c3a96a
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - seqkit=0.16.1
RUN mkdir -p /conda-envs/b9f098cecf23a67fb2bdcd5108c3a96a
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/seqkit.yaml /conda-envs/b9f098cecf23a67fb2bdcd5108c3a96a/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/sratools.yaml
#   prefix: /conda-envs/99ca6a3dbe53ce45720f1e4182c33767
#   name: sra-tools
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - sra-tools = 2.9.1
#     - pigz >= 2.6
#     - pbzip2 >= 1.1
#     - snakemake-wrapper-utils=0.3
RUN mkdir -p /conda-envs/99ca6a3dbe53ce45720f1e4182c33767
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/sratools.yaml /conda-envs/99ca6a3dbe53ce45720f1e4182c33767/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/zscore_normalize_bw.yaml
#   prefix: /conda-envs/3d4579fb17a22f4e47c0d2a5ea7fae3c
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - r-tidyverse = 1.3.1
#     - bioconductor-rtracklayer = 1.50.0
#     - bioconductor-genomicranges = 1.42.0
RUN mkdir -p /conda-envs/3d4579fb17a22f4e47c0d2a5ea7fae3c
ADD https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/envs/zscore_normalize_bw.yaml /conda-envs/3d4579fb17a22f4e47c0d2a5ea7fae3c/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/fqgrep.yaml
#   prefix: /conda-envs/eac6a95815f7ebe29516a1e1ef31fe7a
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fqgrep=1.0.1
RUN mkdir -p /conda-envs/eac6a95815f7ebe29516a1e1ef31fe7a
ADD https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/fqgrep.yaml /conda-envs/eac6a95815f7ebe29516a1e1ef31fe7a/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/seqkit.yaml
#   prefix: /conda-envs/03e0a7036026645725391bede1883a49
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - seqkit=2.3.1
RUN mkdir -p /conda-envs/03e0a7036026645725391bede1883a49
ADD https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/seqkit.yaml /conda-envs/03e0a7036026645725391bede1883a49/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/zscore_normalize_bw.yaml
#   prefix: /conda-envs/cb5325ac8805707ab7f5fcb32a17793f
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - r-base = 4.2.2
#     - r-tidyverse = 1.3.2
#     - bioconductor-rtracklayer = 1.58.0
#     - bioconductor-genomicranges = 1.50.0
RUN mkdir -p /conda-envs/cb5325ac8805707ab7f5fcb32a17793f
ADD https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/envs/zscore_normalize_bw.yaml /conda-envs/cb5325ac8805707ab7f5fcb32a17793f/environment.yaml

# Conda environment:
#   source: https://github.com/tjgibson/NGS-workflow-RNAseq/raw/main/workflow/envs/DEseq2.yaml
#   prefix: /conda-envs/6a44d6a5ba14fe77f1b0729bbe706dc6
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - r-base >= 4.2
#     - r-tidyverse >= 1.3.1
#     - bioconductor-deseq2
#     - bioconductor-rtracklayer >= 1.50.0
RUN mkdir -p /conda-envs/6a44d6a5ba14fe77f1b0729bbe706dc6
ADD https://github.com/tjgibson/NGS-workflow-RNAseq/raw/main/workflow/envs/DEseq2.yaml /conda-envs/6a44d6a5ba14fe77f1b0729bbe706dc6/environment.yaml

# Conda environment:
#   source: workflow/envs/meme.yaml
#   prefix: /conda-envs/b00e8b3d503bd0d2e4a89b646eb00d01
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - meme = 5.4.1
RUN mkdir -p /conda-envs/b00e8b3d503bd0d2e4a89b646eb00d01
COPY workflow/envs/meme.yaml /conda-envs/b00e8b3d503bd0d2e4a89b646eb00d01/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/ccff920ef2dca61d7d93dbd598aef221 --file /conda-envs/ccff920ef2dca61d7d93dbd598aef221/environment.yaml && \
    mamba env create --prefix /conda-envs/16114a9972039cb62b9656ba636b27b0 --file /conda-envs/16114a9972039cb62b9656ba636b27b0/environment.yaml && \
    mamba env create --prefix /conda-envs/c2a8201b8d46efaba69ff3379bedbb7b --file /conda-envs/c2a8201b8d46efaba69ff3379bedbb7b/environment.yaml && \
    mamba env create --prefix /conda-envs/9a5b04e8a8a154f32771f0a56dc701bf --file /conda-envs/9a5b04e8a8a154f32771f0a56dc701bf/environment.yaml && \
    mamba env create --prefix /conda-envs/9608721699f97513ba7f47bd4e3db24b --file /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml && \
    mamba env create --prefix /conda-envs/b33f755bffdf8829e59d1023684a2b44 --file /conda-envs/b33f755bffdf8829e59d1023684a2b44/environment.yaml && \
    mamba env create --prefix /conda-envs/d4a8e62782e6d0fd8b849811a698b9d8 --file /conda-envs/d4a8e62782e6d0fd8b849811a698b9d8/environment.yaml && \
    mamba env create --prefix /conda-envs/a683b7068d9e2ce45ba985ef66cc1c89 --file /conda-envs/a683b7068d9e2ce45ba985ef66cc1c89/environment.yaml && \
    mamba env create --prefix /conda-envs/6f7502c45f8f32d3c16f398e1e3cfbff --file /conda-envs/6f7502c45f8f32d3c16f398e1e3cfbff/environment.yaml && \
    mamba env create --prefix /conda-envs/165a6884838c8278606bcee5ee86ab20 --file /conda-envs/165a6884838c8278606bcee5ee86ab20/environment.yaml && \
    mamba env create --prefix /conda-envs/71fb0136399bd6f669cff3c445df6818 --file /conda-envs/71fb0136399bd6f669cff3c445df6818/environment.yaml && \
    mamba env create --prefix /conda-envs/da0a0a3f54e418d0bf464e31b190514c --file /conda-envs/da0a0a3f54e418d0bf464e31b190514c/environment.yaml && \
    mamba env create --prefix /conda-envs/b9f098cecf23a67fb2bdcd5108c3a96a --file /conda-envs/b9f098cecf23a67fb2bdcd5108c3a96a/environment.yaml && \
    mamba env create --prefix /conda-envs/99ca6a3dbe53ce45720f1e4182c33767 --file /conda-envs/99ca6a3dbe53ce45720f1e4182c33767/environment.yaml && \
    mamba env create --prefix /conda-envs/3d4579fb17a22f4e47c0d2a5ea7fae3c --file /conda-envs/3d4579fb17a22f4e47c0d2a5ea7fae3c/environment.yaml && \
    mamba env create --prefix /conda-envs/eac6a95815f7ebe29516a1e1ef31fe7a --file /conda-envs/eac6a95815f7ebe29516a1e1ef31fe7a/environment.yaml && \
    mamba env create --prefix /conda-envs/03e0a7036026645725391bede1883a49 --file /conda-envs/03e0a7036026645725391bede1883a49/environment.yaml && \
    mamba env create --prefix /conda-envs/cb5325ac8805707ab7f5fcb32a17793f --file /conda-envs/cb5325ac8805707ab7f5fcb32a17793f/environment.yaml && \
    mamba env create --prefix /conda-envs/6a44d6a5ba14fe77f1b0729bbe706dc6 --file /conda-envs/6a44d6a5ba14fe77f1b0729bbe706dc6/environment.yaml && \
    mamba env create --prefix /conda-envs/b00e8b3d503bd0d2e4a89b646eb00d01 --file /conda-envs/b00e8b3d503bd0d2e4a89b646eb00d01/environment.yaml && \
    mamba clean --all -y

# Step 3: install base R
RUN mamba install -c conda-forge r-base=4.2.1 && \
    mamba clean --all -y

# Step 4: install renv
ENV RENV_VERSION 0.17.3
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Step 5: install necessary R packages
WORKDIR /project
COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library
RUN R -e "renv::restore()"