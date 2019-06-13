FROM rocker/verse:3.6.0

RUN install2.r --error \
  mcclust \
  lme4 \
  corrplot \
  DescTools \
  igraph \
  gridExtra
