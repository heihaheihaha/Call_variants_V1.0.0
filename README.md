# Call_variants_V1.0.0
Preliminary test passed on Codespace.

Please modify the ___config.yaml___ befor using.

## Prepare the environment
For the limitation of file size upload to Github repo, you need to download the test files and GATK
```bash
mkdir ./bin
cd ./bin
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip gatk-4.4.0.0.zip
rm gatk-4.4.0.0.zip
```

## Useful commands 
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Silent installation may cause some path issues, so manual installation is recommended.

## Test data
You can use the data of Raredisease pipeline of nf-core, which contain the testdata and reference gonome.
```bash
git clone -b raredisease --single-branch https://github.com/nf-core/test-datasets.git
```
## Run test
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
snakemake --cores 2 --use-conda --conda-create-envs-only
snakemake --cores 2 --use-conda 
```
