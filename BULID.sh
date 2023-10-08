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
perl /data/wangty/WES_pipeline_V2.2/lib/get_result.pl \
                /data/wangty/seq_data/Xiongyin/result/variants/Xiongyin_L02.strict.xls \
                /data/wangty/seq_data/Xiongyin/result/variants/Xiongyin_L02.clinvar_HGMD.xls \
                /data/wangty/seq_data/Xiongyin/result/variants/Xiongyin_L02.loose.xls \
                /data/wangty/seq_data/Xiongyin/result/variants/Xiongyin_L02.final.xls \
                /data/wangty/seq_data/Xiongyin/result/variants/Xiongyin_L02.stat.xlsx