# Applying statistical thinking to palaeo data through GAMs

This repo contains the slide deck and R code for my session keynote at INQUA 2023 in Roma.

[![DOI](https://zenodo.org/badge/665432680.svg)](https://zenodo.org/badge/latestdoi/665432680)

Please cite as: Simpson, Gavin L. (2023, July 13-20) _Applying statistical thinking to palaeo data through GAMs_ [Session keynote]. XXI Congress of the International Union for Quaternary Research (INQUA), Rome, Italy. [10.5281/zenodo.8186141](https://doi.org/10.5281/zenodo.8186141)

### Gavin Simpson

### Aarhus University, Denmark

Palaeoecological and palaeoenvironmental stratigraphic data exhibit a number of
features that make them complex to analyse; for example they

1. are autocorrelated in time (and, potentially, also in space),
2. are generally irregularly spaced in time, and
3. have complex variance properties due to time averaging and sampling
processes.

This complexity has some unfortunate consequences. All too often, data are not
subject to any form of statistical analysis, and, when they are, often that
analysis uses inappropriate statistical methods. Furthermore, analysts will
often invent their own *ad hoc* approaches to get answers to their specific
questions. The statistical properties of such *ad hoc* approaches are rarely
known, however.

A number statistical methods have been developed over recent years and decades
that can handle the complex needs of palaeo data. Among these methods is the
generalized additive model (GAM) and related extensions. The use of GAMs to
model temporal trends in palaeo data has increased steadily in recent years,
thanks, in part, to the availability of high quality, open source software like
the R package mgcv. GAMs are relatively simple models that will be familiar to
anyone who has taken a graduate course in ecological statistics, and yet GAMs
are incredibly versatile models that can be adapted to suit many pressing
questions that palaeo scientists wish to address.

In this talk I will briefly introduce GAMs and show how penalised splines, which underpin the GAM, work, before moving on to cover some recent developments from
my own work that extend GAMs to new areas and questions. In particular, I will
discuss;

i. models that go beyond the mean to investigate trends in other moments, such
as variances of palaeo time series,
ii. using derivatives of splines fitted in a GAM as estimates of rates of
change in palaeo time series, and
iii. copula GAMs that can estimate how the correlation between two variables
itself changes over time.

In each case I present a relevant motivating example and supplement the talk
with fully worked examples in R that are available on GitHub.
