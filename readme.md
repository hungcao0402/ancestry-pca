# Project Name
> Ancestry PCA

## Description
Ancestry analysis and evaluate the similarity or difference between populations or superpopulations based on 1000 genomes data by using principal component analysis

method - PCA in plink tool to reduce dimension (10 components) and select 2 principal components that maximizes the variance explained to visualize results.

- Data: Chromosome 20 genotypes (vcf file) and mapping population (csv file). Default directory: database/

- Information: The 1kGP data provide 2504 samples which were collected with specific populations as follows:

    | Superpopulation | Population | Details |
    | --- | --- | --- |
    | EAS | CHB | Han Chinese |
    |     | JPT | Japanese |
    |     | CHS | Southern Han Chinese |
    |     | CDX | Dai Chinese |
    |     | KHV | King Vietnamese |
    |     | CHD | Denver Chinese |
    | --- | --- | --- |
    | EUR | CEU | CEPH |
    |     | TSI | Tuscan |
    |     | GBR | British |
    |     | FIN | Finnish |
    |     | IBS | Spanish |
    | --- | --- | --- |
    | AFR | YRI | Yoruba |
    |     | LWK | Luhya |
    |     | GWD | Gambian |
    |     | MSL | Mende |
    |     | ESN | Esan |
    | --- | --- | --- |
    | AMR | ASW | African-American SW |
    |     | ACB | African-Caribbean |
    |     | MXL | Mexican-American |
    |     | PUR | Puerto Rican |
    |     | CLM | Colombian |
    |     | PEL | Peruvian |
    | --- | --- | --- |
    | SAS | GIH | Gujarati |
    |     | PJL | Punjabi |
    |     | BEB | Bengali |
    |     | STU | Sri Lankan |
    |     | ITU | Indian |


- Output: 2D plot with first and second variance explained components. Default directory: output_figure/

## Usage:
1. Create conda environment with required packages:
```bash
conda env create -f ancestry.yml
conda activate ancestry
```
2. Run ancestry_pca.py:
```bash
# examples:

# with -s (superpopulation)
python ancestry_pca.py -s "EUR EAS AFR SAS AMR"
# with -p (population)
python ancestry_pca.py -p "ASW ACB KHV JPT STU CEU ESN"
# check percentile ancestry of samples, example: HG01595 HG01607 HG03702 HG03078
python ancestry_pca.py -s "EUR EAS AFR SAS AMR" -u "HG01595 HG01607 HG03702 HG03078"

# Output .png file in default directory 'output_figure'
```

- Flags:
```bash
    # Note: choose only population (-p) or superpopulation (-s)
    
    # choose population to plot
    -p --population
    # choose superpopulation to plot
    -s --superpopulation
    # check percentile ancestry of samples
    -u --user
    
    # view all flags
    python ancestry_pca.py --help
```