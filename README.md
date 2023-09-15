# Call_variants_V1.0.0
Preliminary test passed on Codespace.

Please modify the ___config.yaml___ befor using.

## Bugs to fix:
1. mamba environment need to be bulid previously, use the first_step_mamba.yml to bulid the env and then add the path to the env to the config.ymal
2. To run the whole workflow, a sentry file has been added to the process. Run `snakemake --cores 2 ../results/variants/.sentinel`
3. How to deal with many files at a time

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
## Run test on codespace
```bash 
bash BULID.sh
cd
##
bash Miniconda3-latest-Linux-x86_64.sh
conda init
```
Start a new terminal
```bash
cd
bash Mambaforge-Linux-x86_64.sh
conda init
```
Start a new terminal
```bash
mamba init
```
```
mamba env create -f first_step_mamba.yml
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
snakemake --cores 2 --use-conda
```
