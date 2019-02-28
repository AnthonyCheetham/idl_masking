# idl_masking
IDL pipeline for processing and analysing Sparse Aperture Masking data.

This pipeline was developed by a lot of people over many years, and this represents one of many forks made over the years. Credit goes to the major developers, including Peter Tuthill, Mike Ireland and John Monnier.

## How do I use it?

Example scripts to process data are included in the main directory. The basic reduction is split into 4 steps:

1. Data cleaning. Performs sky-subtraction, centring, flat fielding, bad pixel / cosmic ray correction, and gets relevant information from the headers. Relevant routines are fizeau/sphere/qbe_sphere.pro , fizeau/conica/qbe_conica.pro, fizeau/gpi/qbe_gpi.pro etc.
2. Measurement of interferometric quantities (closure phase, square-visibility, with limited support for bispectral amplitudes and complex visibilities). This is done through direct measurement from the Fourier transform of each image with a matched-filter approach. Relevant routines are calc_bispect.pro and calc_bispect_gpi.pro (for 4D data).
3. Calibration of interferometric quantities. For traditional SAM data, calibrate_v2_cp.pro and calibrate_v2_cp_gpi.pro (for 4D data) are used. These take weighted averages of calibrator quantities and divide (for amplitudes) or subtract (for phases) from those measured on the target.


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.