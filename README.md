## General information
This repo contains all the scripts (and some data) for the following paper: 

**Spectral partitioning identifies individual heterogeneity in the functional network topography of ventral and anterior medial prefrontal cortex (2020). ([NeuroImage](https://authors.elsevier.com/a/1ZzuD3lc~r74YX))**

The paper can be reproduced in pdf format by running `Manuscript.Rmd` in R. The markdown file is organized so that code for the analyses precedes each written section. The original analyses were run in R 3.4, but also work in R 3.5.

## Reproduce the preprint
If you would like to actually run the file, you can download a partial version of this repo that contains the data and essential files from: [https://osf.io/4x6z8/](https://osf.io/4x6z8/). Once the directory ("mPFC_partitioning") is downloaded, you can simply run `Manuscript.Rmd` (provided all the required libraries are installed). **Note that the data are somewhat heavy, and the first run through will be slow. A cache will be produced the first time the paper file is run (this will add multiple GBs to the directory), so subsequent runs will be much faster.** 

### Docker environment

Rather than manually installing the R packages I used, you can load a Docker image with the necessary dependencies to run `Manuscript.Rmd`. This has the advantage of ensuring that package versions are compatible with the code we wrote. If you don't know what [Docker](https://www.docker.com/) is, and how you can share computing environments for reproducibility, go check it out! And take a look at [this tutorial](https://ropenscilabs.github.io/r-docker-tutorial/) for an easy way to set up your own images. 

If you have Docker installed, here are instructions on how to set the environment up:

- Store the directory from OSF wherever you want (I tend to leave it on the desktop to make it easier)

- Run the following command, where `<yourpath>` is the path to the saved folder from OSF (this will download a ~3GB image to your machine):

```
docker run --rm -p 8787:8787 -e PASSWORD=mPFC -v </yourpath>/mPFC_partitioning:/home/rstudio/ ctoroserey/mpfc_partitioning:preprintenv
```

- This will allow you to load the data from that directory, and everything produced within the Docker container will be stored locally within the mPFC_partitioning dir. Now, to load Rstudio go to your browser and type this on your URL bar: `<yourIPaddress>:8787`. Rstudio will now load (user will be "rstudio" and password "mPFC"), and you will see the local files loaded along with it.

- Alternatively, you can download the Docker image in a tar file (mpfc_partitioning.tar from OSF), and use `docker load --input mpfc_partitioning.tar` to access the environment.

You will now be able to play with the data you downloaded from OSF, and knit the preprint pdf. Just beware that running the script can be slow, as it's tied to your local hardware's performance. Note that, for the sake of storage, this environment does not include the software to run the preprocessing or volumetric projection from the meta-analysis. It is mainly provided to allow anyone to knit the preprint and reproduce the group-level statistics reported in the paper. If you would like to test spectral partitioning, see further  below.

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

We have put together a short script that will run both Spectral Partitioning and Modularity for a single individual's data (according to the specifications mentioned in the paper): `communityDetection.R`. You'll just have to install `igraph` in R.

While more straightforward, this script still assumes that the fMRI time series are available in a csv file, and comply with HCP surface conventions (at a density of 32K, specifically).

If your data are in a 32k CIFTI format, you can use the HCP's `wb_command` utilities to create a csv file with just cortical time series.

```
wb_command -cifti-convert -to-text <subj.dtseries.nii> <temp_subj_timeSeries.csv>

awk 'NR <= 59412 { print }' temp_subj_timeSeries.csv >> subj_timeSeries.csv
```
Then, to run `communityDetection.R`, just use `Rscript communityDetection.R <subj>` (the file and script must be in the same directory). The script will output a simple summary, as well as SP and modularity text files ready to be converted to CIFTI. You can use `vector2Cifti.sh` to transform the resulting text files to CIFTI .dscalar.nii files that can be viewed on the HCP's `wb_view`.
