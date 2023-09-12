mkdir ./bin
cd ./bin
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip gatk-4.4.0.0.zip
rm -rf gatk-4.4.0.0.zip
cd
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

git clone -b raredisease --single-branch https://github.com/nf-core/test-datasets.git
# mamba create -c conda-forge -c bioconda -n snakemake snakemake
# mamba env create -f first_step_mamba.yml
