# *Plasmodium vivax* AmpliSeq Analysis Pipeline - BETA

## Background

### Purpose
This git repo hosts a **BETA** version of code for analysis of *Plasmodium vivax* data generated using an AmpliSeq design that amplifies 496 targets covering 415 unique targets across all 14 chromosomes of the *P. vivax* genome ([BED File](./AmpliSeq_Design.bed)). Please read through this README for instructions on running this workflow. The analysis goal is the following:
    
1) Identify 'same-strain' infections of *P. vivax* (in other words, closely related clusters),

2) Prediction of *P. vivax* samples to a geographic region,

3) Identify SNPs in proposed drug-resistance markers in *P. vivax*.

Please cite the following paper when referencing this workflow:

`Barratt et al 2024. Genetic characterization of Plasmodium vivax linked to autochthonous malaria transmission in the US (2023) using Illumina AmpliSeq technology: a genetic epidemiology study by Dr Joel Barratt. Submitted to The Lancet Microbe`

### Contact
Bioinformatic support: David Jacobson (quh7@cdc.gov)

*P. vivax* genotyping inquiries: malarialab@cdc.gov

General questions CDC's Domestic Parasite Surveillance: Joel Barratt (nsk9@cdc.gov)


## Requirements for Running the Pipeline

### Software Dependencies
**Nextflow**

Our analysis workflow uses nextflow to execute commands, see [here](https://www.nextflow.io/docs/latest/index.html) for more information on installing and running nextflow. The pipeline is currently running on nextflow version `23.10.0` and uses [nextflow DSL 2](https://www.nextflow.io/blog/2020/dsl2-is-here.html), which was introduced in version `20.07.1`.

The nextflow config file is in the `pvivax_ampliseq_analysis/nextflow_configFiles` folder. It contains the default parameters (e.g., RAM, threads, haplotype calling depth) and paths to a variety of reference files and python/R scripts. The reference parameters should not need to be changed if this repo is cloned from GitHub.

One line of the nextflow.config file may need to be changed to fit running in your working environment: `singularity.runOptions = "-B /scicomp:/scicomp"`. You will likely need to either comment this line out with `//` or change this line to match your HPC. The code should run fine if the line is commented out but contact your system administrator if you have questions/issues.


**Singularity**

All other software dependencies are in a series of singularity containers, see [here](https://docs.sylabs.io/guides/latest/user-guide/). Please note that the containers were built with singularity version `3.8` and have not been tested on more recent singularity versions (e.g. versions `4.0` and above).

The [nextflow config](nextflow_configFiles/nextflow.config) file uses the `withLabel` syntax to specify the container used for different processes. The default option is to use downloaded versions of the container; however, nextflow should be able to download the singularity container hosted on sylabs when the nextflow script is first executed. See example below, where the `balkClassifer` label is used for a specific SIF container and the option for downloading the container from a hosted site is commented out with `//`

    withLabel: balkClassifier {
        container = "/scicomp/groups-pure/Projects/CycloDeploy/Plasmodium_Ampliseq_Nextflow/singularityContainers/balkClassifier.sif"
        //container = "library://davejacobson12/pvivax_containers/balkclassifier:latest"
      }

If you want to download the containers and store local versions, please download the containers from sylabs: https://cloud.sylabs.io/library/davejacobson12. There are options for pulling from the command line or downloading from the browser.

**If you download the containers you must update the nextflow.config file to point to the containers you have downloaded. The full path to the local containers must be specified.**


**Additional Scripts**

Python and R scripts used in the pipeline are stored in the `pvivax_ampliseq_analysis/python_R_scripts` folder. Users should not need to access/edits these scripts.

**Data Format**

FASTQ files must be paired, gzip compressed, and end in R1_001.fastq.gz / R2_001.fastq.gz for forward and reverse reads. The default location to place FASTQ files for entry into the workflow is the `TEST_DATA`` folder.

### Additional Installation Steps

**Train Geographic Prediction Model**

*TBD*

**Reference Database Creation**

A handful of reference files need to be downloaded/unzipped/made for running the pipeline. Follow the directions below.

Navigate to the `pvivax_ampliseq_analysis/REFERENCES/pv_geo_refs/` subfolder and run the following code. This will create a bowtie2 database for each target amplicon.

    ls *fasta | awk -F "." '{print$1}' | while read line; do bowtie2-build $line.fasta $line\_bt2 ; done

Navigate to the `pvivax_ampliseq_analysis/REFERENCES` folder and download and then unzip a human reference genome (instructions borrowed from [metagenomics.wiki](https://www.metagenomics.wiki/tools/short-read/remove-host-sequences)).

    wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip

    unzip GRCh38_noalt_as.zip

Navigate to the `pvivax_ampliseq_analysis/HAPLOTYPE_CALLER_PVIVAX/` folder and unzip the REF_SEQS folder

    unzip REF_SEQS.zip

Navigate to the `pvivax_ampliseq_analysis/Eukaryotpying-Python-main/` folder and unzip the haplotype reference database.

    unzip PVIVAX_REF_GENOTYPES.zip

## Running the Workflow
### Options

Within the `main_pvivx_full.nf` there are different workflows that combine various elements of the pipeline.

**Option 1**.  Run the complete workflow: read in FASTQs, drug resistance surveillance (MaRS), geographic prediction, and haplotype calling + clustering all at once.

    nextflow run main_pvivax_full.nf  -entry everything -c nextflow_configFiles/nextflow.config -profile singularity -with-report nextflow_logFiles/report.html -with-trace nextflow_logFiles/trace.txt


**Option 2**. Mars and GeoPrediction only - Call SNPs at MaRS markers for *P. vivax* (results in `MaRS_output`). Predict geographic origin for P. vivax. Region/country prediction will be merged with metadata (results in `GeoPrediction_output/predictedOut`)

    nextflow run main_pvivax_full.nf  -entry P_vivax_MaRS_GeoPrediction -c nextflow_configFiles/nextflow.config -profile singularity -with-report nextflow_logFiles/report.html -with-trace nextflow_logFiles/trace.txt

**Option 3A**. Haplotype calling without clustering. Can be used if you want to batch process the haploytpeing steps and then run distance matrix + clustering at the very end.

    nextflow run main_pvivax_full.nf -entry Workflow1_v2_vivax -c nextflow_configFiles/nextflow.config -profile singularity -with-report nextflow_logFiles/report.html -with-trace nextflow_logFiles/trace.txt

**Option 3B**. Generate haplotype sheeting without hap calling. If you need to generate a haplotype sheet from a specific set of samples and/or references. Use the `--hapSheet_specifcGenotypes` flag to give path to a folder with specific samples to incluede in the haplotype sheet. Use the `--hapSheet_specifcReferences` flag to give path to a folder with specific references to incluede in the haplotype sheet.

    nextflow run main_pvivax_full.nf  -entry hapSheetOnly -c nextflow_configFiles/nextflow.config  -profile singularity -with-report nextflow_logFiles/report.html -with-trace nextflow_logFiles/trace.txt

**Option 3C**. Distance matrix + clustering without haplotype calling. Use if you do not need to perform haplotype calling on any additional specimens and only need to perform clustering. By default, it will use the most recently generated haplotype sheet in the `haplotype_sheets`. Use the `--wf2HapSheet` flag to analyze a specific haplotype sheet.

    nextflow run main_pvivax_full.nf  -entry Workflow2and3_Pvivax -c nextflow_configFiles/nextflow.config -profile singularity -with-report nextflow_logFiles/report.html -with-trace nextflow_logFiles/trace.txt

### Test the Workflow Installation

Contact David Jacobson (quh7@cdc.gov) or the CDC malaria genotyping inbox (malarialab@cdc.gov) to request data for testing your workflow installation.

### Workflow Outputs

The primary workflow outputs are below

**Haplotype Calling and Clustering**

1) A genotype for each specimen will be in the `HAPLOTYPE_CALLER_PVIVAX/SPECIMEN_GENOTYPES`
   
2) A haplotype sheet with genotypes for all specimens in the `SPECIMEN GENOTYPES` and `PVIVAX_REF_GENOTYPES` folder
   
3) A pairwise distance matrix will be in the `ensemble_matrices` folder
   
4) Cluster memberships for each specimen will be in the `clusters_detected` folder

**Geographic Prediction**

1)  A predicted geographic result will be in the `GeoPrediction_output/predictedOut/` folder. We strongly recommend using the predictions in the file with *gatkRegion* in the title.
   
2)  VCF for each sample at all SNPs used in geographic prediction will be in the `GeoPrediction_output/variantCalls/` folder

**Drug Resistance Surveillance**

1) A variety of processing files will be output as part of the MaRS drug resistance surveillance pipeline - all of which will be in the `MaRS_output/` folder
   
2) The drug resistance primary outputs of interest will be in `MaRS_output/Summary` folder 

**Reports and Trees**

The workflow is not currently set up to generate reports/trees using this GitHub Repo. However, there are some scripts for report/tree generation in the `python_R_scripts/reporting_scripts/` folder. These scripts are not maintained.


**Check the readme.txt files**

The following folders are pre-populated with txt files describing the contents of that folder. Read those files to ascertain the contents of each folder when the workflow runs.


    clusters_detected
    ensemble_matrices
    ensemble_matrices_perChromosome
    GeoPrediction_Out (all subdirs)
    MaRS_out no read me
    SPECIMEN_GENOTYPES
    cleanReads
    haplotype_sheets
    nextflow_logFiles
    reportFolder
    treeFolder

**Template for clearance: This project serves as a template to aid projects in starting up and moving through clearance procedures. To start, create a new repository and implement the required [open practices](open_practices.md), train on and agree to adhere to the organization's [rules of behavior](rules_of_behavior.md), and [send a request through the create repo form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) using language from this template as a Guide.**

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Access Request, Repo Creation Request

* [CDC GitHub Open Project Request Form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) _[Requires a CDC Office365 login, if you do not have a CDC Office365 please ask a friend who does to submit the request on your behalf. If you're looking for access to the CDCEnt private organization, please use the [GitHub Enterprise Cloud Access Request form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUQjVJVDlKS1c0SlhQSUxLNVBaOEZCNUczVS4u).]_




I need to update the report generation aspects of the code in a way that make sense for newest versions of code

pvivax metadata xlsx sheet, maybe have an empty format?
    pvivax_metadata.xlsx is the default in the config. Uploaded example under pvivax_metadata_github.xlsx name

Explain that region prediction results are best for pvivax geo prediction

Explain where results will go

## Code examples

Copy options from qsub script?



## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

## Overview

Describe the purpose of your project. Add additional sections as necessary to help collaborators and potential collaborators understand and use your project.
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
