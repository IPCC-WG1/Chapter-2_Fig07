##########################################################################
# ---------------------------------------------------------------------------------------------------------------------
# This is [insert software name (e.g., Python)]* code to produce IPCC AR6 WGI Figure/s [Figure 2.7 ] *
# Creator: [Mark Weber@uni-bremen.de; Vitali.fioletov@canada.ca; melissa.gomis@universite-paris-saclay.fr ] 
# Contact: [Mark Weber@uni-bremen.de]*
# Last updated on: 01.11.2020
# --------------------------------------------------------------------------------------------------------------------
#
# - Code functionality:  “bams_oz_zonalmean_v2_1965.pro, plots zonal mean total ozone columns in 6 regions. The produced plots have subsequently been modified for graphical presentation in Illustrator, without modification of data.”] *
# - Input data: [gb_1964-2020_za_final.txt (WOUDC, 1964-2020, Vitali Fioletov <vitali.fioletov@ec.gc.ca>), sbuv.v87.mod_v12.70-20.za.txt (NASA, 1970-2020, Stacey Frith <Stacey.M.Frith@nasa.gov>), total_o3_7820.dat (NOAA, 1978-2020, Jeannette Wild <jeannette.wild@noaa.gov>), GSG_merged_zonalmean.dat (GSG, 1995-2020, Mark Weber <weber@uni-bremen.de>), GTO-ECV_CCI_OMI_ext.zonalmean.1995-2020.dat (GTO, 1995-2020, Melanie Coldewey-Egbers <melanie.coldewey-egbers@dlr.de>), toc_zonal_monthly_mean_1970_2020_msr2.txt (MSR2, 1970-2020, Ronald van der A <avander@knmi.nl>)*
# - Output variables:  The code plots the figure as in the report.
# - spatial scale: global
# - source data type: observations, monthly mean zonal mean total ozone column (5 degs step)
# - plot data type: Near global 60S-60N;Northern Hemisphere (35N-60N);Tropics (20N-20S), Southern Hemisphere (35 S-60 S), Polar Northern Hemisphere (60-90 N), Polar Southern Hemisphere (60S-90S)
# ----------------------------------------------------------------------------------------------------
# Information on  the software used
# - Software Version: [IDL8.51 and libraries]
# - Landing page to access the software:  https://www.l3harrisgeospatial.com/Software-Technology/IDL
# - Operating System: [OS independent]*
# - Environment required to compile and run: [environement independent]*
#  ----------------------------------------------------------------------------------------------------
#
#  License: [ Default license: Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)]*
#
# ----------------------------------------------------------------------------------------------------
# How to cite: https://doi.org/10.5281/zenodo.6353844
# When citing this code, please include both the code citation and the following citation for the related report component:
Weber M., Fioletov, V., van der A, R., Coldewey-Egbers, M., Frith, S. M., Wild, J. (2021); bams_oz_zonalmean_v2_1964.pro, IDL program to produce graph in BAMS State of the Climate 2020 report: https://doi.org/10.1175/2020BAMSStateoftheClimate_Chapter8.1 
########################################################################
