Fitting of diffractive data

# Plan
* Crosscheck the 2006 DPDF (i.e. in NLO FFNS)
* Repeat the fit to inclusive data with VFNS (NLO)
* Repeat the fit to inclusive data with VFNS (NNLO)
* Include jet data to the fit
* Do combined fit at NLO vs at NNLO (chi2 difference?)


# Crosscheck of the 2006 DPDF dit

## The check of the chi2 with published F2, FL
Agrees within 1 unit of chi2  
chi2 /ndf = 157.043 / 184 = 0.853497 (fit A) : org val 158.0154  
chi2 /ndf = 163.366 / 184 = 0.887861 (fit B) : org val 164.5480   
The small difference can be caused by rounding errors of the provided published values from the text file.
Or some uncertainty is slightly different.

## Check of the evolution
Evolution checked against QCDNUM.   
Notice the poor z-grid interpolation at the high-z in the original H1 2006 fits.
The input parametrisation must include steps (i.e. not the analytical formula) to get the evolution consistent.


## Check of the convolution
Convolution checked using QCDNUM ZMSTF and HQSTF packages.  
The difference observed in F2c and F2l at the higher scales, why?
