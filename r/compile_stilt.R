# STILT makefile for first-time installation of STILT-R version 2, DW, 06/28/2020

# install dependencies for STILT
if (!require('devtools')) install.packages('devtools')
devtools::install_github('benfasoli/uataq')
uataq::stilt_init('stilt_hysplit', branch = 'hysplit-merge')
