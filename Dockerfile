# file created based on https://colinfay.me/docker-r-reproducibility/

# rocker houses callable versions of base R
FROM rocker/r-ver:3.5.2

# use this if you'd like to enter a desired data every time the image is called
# ARG WHEN

# otherwise set the specific date
WHEN="2019-05-20/"

# dir to house container
RUN mkdir /home/paper

# instructions to load CRAN packages from that $WHEN moment in time
RUN R -e "options(repos = \
  list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/${WHEN}')); \
  install.packages(c('tidyverse', 'mcclust', 'knitr'))"

# take the manuscript file to the container dir
COPY Manuscript.Rmd /home/paper/Manuscript.Rmd

# and source it
CMD R -e "source('/home/paper/Manuscript.Rmd')"

## how to use
# to build (name: paper)
#docker build --build-arg WHEN=2019-01-06 -t paper . # WHEN only if undefined previously

# and launch
#docker run paper 