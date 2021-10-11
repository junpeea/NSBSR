########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_ImportLibrary

# Date : 201811xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

### Install Fundamental Libraries
if(!require(spectralGP))   {install.packages("spectralGP");require(spectralGP)}
if(!require(RandomFields)) {install.packages("RandomFields");require(RandomFields)}
if(!require(geoR))         {install.packages("geoR");require(geoR)}
if(!require(coda))         {install.packages("coda");require(coda)}
if(!require(fftwtools))    {install.packages("ffwtools");require(fftwtools)}
if(!require(spam))         {install.packages("spam");require(spam)}
if(!require(matrixcalc))   {install.packages("matrixcalc");require(matrixcalc)}
if(!require(MASS))         {install.packages("MASS");require(MASS)}
if(!require(mvtnorm))      {install.packages("mvtnorm");require(mvtnorm)}
if(!require(fields))       {install.packages("fields");require(fields)}
if(!require(dplyr))        {install.packages("dplyr");require(dplyr)}



