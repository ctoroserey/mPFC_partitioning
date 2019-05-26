## General information
This repo contains all the scripts (and some data) for the following paper: 

**Individual heterogeneity in the functional network topography of medial prefrontal cortex (2019)**

The paper can be reproduced in pdf format by running `Manuscript.Rmd` in R. The markdown file is organized so that code for the analyses precedes each written section. The original analyses were run in R 3.4, but also work in R 3.5.

## Reproduce the preprint
If you would like to actually run the file, you can download a partial version of this repo that contains the data and essential files from: [https://osf.io/4x6z8/](https://osf.io/4x6z8/). Once the directory ("mPFC_partitioning") is downloaded, you can simply run `Manuscript.Rmd`. **Note that the data are somewhat heavy, and the first run through will be slow. A cache will be produced the first time the paper file is run (this will add multiple GBs to the directory), so subsequent runs will be much faster.** As of now you would have to install the required packages, but a Docker image is in the works so a container with the proper package versions can be built.

## Description and order of scripts
The downloaded data are summaries of the already-run individualized partitionings, *not preprocessed HCP time series*. The repo contains the scripts we used to analyze most of the data (preprocessing not included). The order the scripts were run was:

### Meta-analysis

- `metacoordsVol2Surf` will project study-specific peak coordinates onto the Glasser et al., 2016 surface parcellation. To run it, use `sh metacoordsVol2Surf <coordDir>`, where "coordDir" is either metacoordsDN or metacoordsSV (for each literature separately). Note that this requires FSL, Freesurfer, and AFNI, and will take hours to run.

### Individualized community detection

- `graphAnalysis.R` computes all the reported analyses for an individual, and is estimated to take at least a couple of hours to run. The outputs are the summary files found in the osf repository. This version assumes that:

    - The individual underwent comprehensive preprocessing
    - The time series have been transformed to csv
    - Data for each session are stored in independent csv files
    - A non-demeaned version of each run exists

    Naming convensions can be found within the "Load HCP Data" section of the script.

The remainder of the analyses, including the meta-analytic permutation and inter-individual analyses, are run within `Manuscript.Rmd`.

## Try it yourself!

We have put together a short script that will run both Spectral Partitioning and Modularity for a single individual's data (according to the specifications mentioned in the paper): `communityDetection.R`. While more straightforward, this script still assumes that the fMRI time series are available in a csv file, and comply with HCP surface conventions (at a density of 32K, specifically).

If your data are in a 32k CIFTI format, you can use the HCP's `wb_command` utilities to create a csv file with just cortical time series.

```
wb_command -cifti-convert -to-text <subj.dtseries.nii> <temp_subj_timeSeries.csv>

awk 'NR <= 59412 { print }' temp_subj_timeSeries.csv >> subj_timeSeries.csv
```
Then, to run `communityDetection.R`, just use `Rscript communityDetection.R <subj>` (the file and script must be in the same directory). The script will output a simple summary, as well as SP and modularity text files ready to be converted to CIFTI. You can use `vector2Cifti.sh` to transform the resulting text files to CIFTI .dscalar.nii files that can be viewed on the HCP's `wb_view`.
