FROM rocker/verse:3.6.0

RUN install2.r --error \
  --deps TRUE \
  mcclust \
  lme4 \
  parallel \
  corrplot \
  DescTools \
  knitr