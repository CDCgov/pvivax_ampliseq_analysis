# *Plasmodium vivax* AmpliSeq Analysis Pipeline - BETA

## Background

This git repo hosts a BETA version of code for analysis of *Plasmodium vivax* data generated using an AmpliSeq design that amplifies 495 targets covering 415 unique targets across all 14 chromosomes of the *P. vivax* genome ([BED File](./AmpliSeq_Design.bed)). The analysis goal is the following:
    
1) Identify 'same-strain' infections of *P. vivax* (in other words, closely related clusters),

2) Prediction of *P. vivax* samples to a geographic region,

3) Identify SNPs in proposed drug-resistance markers in *P. vivax*.

Please cite the following paper when referencing this workflow:

`Genetic characterization of Plasmodium vivax linked to autochthonous malaria transmission in the US (2023) using Illumina AmpliSeq technology: a genetic epidemiology study by Dr Joel Barratt. Barratt et al 2024. Submitted to The Lancet Microbe`



## Requirements for Running Pipeline

### Software Dependencies
nextflow
    
    update config

singularity

    Access to containers

### Additional Installation Steps

**Train Geographic Prediction Model**

Add in steps for this

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

**Check the readme.txt files**

The following folders are pre-populated with txt files describing the contents of that folder. Read those files to ascertain the contents of each folder when the workflow runs.



**Template for clearance: This project serves as a template to aid projects in starting up and moving through clearance procedures. To start, create a new repository and implement the required [open practices](open_practices.md), train on and agree to adhere to the organization's [rules of behavior](rules_of_behavior.md), and [send a request through the create repo form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) using language from this template as a Guide.**

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Access Request, Repo Creation Request

* [CDC GitHub Open Project Request Form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) _[Requires a CDC Office365 login, if you do not have a CDC Office365 please ask a friend who does to submit the request on your behalf. If you're looking for access to the CDCEnt private organization, please use the [GitHub Enterprise Cloud Access Request form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUQjVJVDlKS1c0SlhQSUxLNVBaOEZCNUczVS4u).]_




readme files in folders



run test data?

probably need to retrain a classifier but not sure? Would need a beta tester. Can provide instructions for creating model

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
