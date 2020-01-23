# anelastic_boussinesq_2D
A simple 2D spectral Anelastic/Boussineq model developed for the Horinouchi et al. (2020) paper on convective bursts and gravity waves

## Author

Takeshi Horinouchi
horinout __at__ 
ees.hokudai.ac.jp

## Description

This is the 2D spectral Anelastic/Boussinesq model developed for the paper
Horinouchi et al. (2020) (see Reference below). It is a simple FFT-based
spectral model on the 2D horizontal-vertical domain. The model is based
on the anelastic approximation but it can be a Boussinesq model.
In this model, the flow is driven by diabatic heating, which is
either prescribed or parameterized (switched by option). See the
reference for what one can do with this model.

## Reference

If you publicate your results by using this model, please refer
the following original paper:

Horinouchi, T., U. Shimada, A. Wada, 2020,
Convective Bursts With Gravity Waves in Tropical Cyclones: Case Study
With the Himawari-8 Satellite and Idealized Numerical Study.
Geophysical Research Letters, art num. GRL60158, doi:10.1029/2019GL086295

## Programing language

Ruby

## Installation

To use the model, you need the class library GPhys 
(http://ruby.gfd-dennou.org/products/gphys/doc/),
which you can install by

```gem install gphys
```

## Usage

```$ ruby boussinesq2D_CBmodel.rb [options]
```

See the source code for available options.

## Source code

* boussinesq2D.rb : Boussinesq version of the model. It defines the class
  NumRu::Boussinesq2D
* anelastic2D.rb : anelastic version of the model It defines the class
  NumRu::Anelastic2D, which is defined by inheriting NumRu::Boussinesq2D
* boussinesq2D_CBmodel.rb : the driver
* fftw_ext.rb : extension of NumRu::FFTW3

## Remark

Initially, the comments in the source code is in Japanese,
but it will be revised to English. Sorry about the inconvenience.
