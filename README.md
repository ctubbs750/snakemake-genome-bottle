# snakemake-genome-bottle

A Snakemake module for installing a reference genome along with its supporting files.

# Notes

This module is designed to be incorporated into an existing Snakemake pipeline for integrative analysis. To use it, you can point towards its Snakefile raw url during module import:

`https://raw.githubusercontent.com/ctubbs750/genome/main/workflow/Snakefile`

 By default, it will install the hg38 human reference genome into the existing project directory under:

 `/resources/data/genome/{build}`

 For this module to function, it will look for parameters under the `GENOME` tag in config file.  Copy and paste the config parameters found in `/config/config.yaml` into the desired project's config file and edit as needed.
