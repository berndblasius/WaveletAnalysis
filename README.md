# WaveletAnalysis
Functions for time continuous wavelet transform and 
bivariate wavelet analysis, implemented in Julia

This code was used in a forthcoming paper: Blasius et al. (2019)


## Time continuous wavelet transorm

This is based on a straightforward Julia translation of the
1d Wavelet transform from Torrence and Compo

There is no intention to have a fancy, performant, or elegant "translation"
it just should be able to "do stuff"
The copyright remains by Torrence and Compo (see below)
This code is provided without any warranties
by Bernd Blasius (bernd.blasius@gmail.com)

 HERE GOES THE ORIGINAL COPYRIGHT
---------------------------------------------------------------------------
  Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo

  This software may be used, copied, or redistributed as long as it is not
  sold and this copyright notice is reproduced on each copy made. This
  routine is provided as is without any express or implied warranties
  whatsoever.

Notice: Please acknowledge the use of the above software in any publications:
   ``Wavelet software was provided by C. Torrence and G. Compo,
     and is available at URL: http://paos.colorado.edu/research/wavelets/''.

Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
           Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.

Please send a copy of such publications to either C. Torrence or G. Compo:
 Dr. Christopher Torrence               Dr. Gilbert P. Compo
 Research Systems, Inc.                 Climate Diagnostics Center
 4990 Pearl East Circle                 325 Broadway R/CDC1
 Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
 E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
---------------------------------------------------------------------------


## Bivariate wavelet analysis

This is code to perform phase analysis of two (or more) time series

