# Call_variants_V1.0.0
Preliminary test passed on Codespace.

Please modify the ___config.yaml___ befor using.

## Bugs to fix:
1. mamba environment need to be bulid previously, use the first_step_mamba.yml to bulid the env and then add the path to the env to the config.ymal
2. To run the whole workfolw, a sentry file has been added to the process. Run `snakemake --cores 2 ../results/variants/.sentinel`

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

## Test data
You can use the data of Raredisease pipeline of nf-core.
```bash
git clone -b raredisease --single-branch https://github.com/nf-core/test-datasets.git
```
