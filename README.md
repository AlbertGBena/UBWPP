# Alternative raw processing for Wind Profiler, PCL1300 - UBWPP.py

UBWPP is a novel Wind Profiler processing methodology developed for Wind Profiler model PCL1300, manufactured by Degreane. UBWPP works with the compresseed spectral power from the wind profiler manufacturer files (.dat files). 
UBWPP produces a number of output fields which include equivalent reflectivity (Ze), Vertical Doppler speed and derived parameters such as spectral width, skewness, kurtosis, horizontal wind components (if the device has more than one beam), and a simplified precipitation type classification (rain, mixed, snow and unknown) with 2 approaches (Atlas et al. 1973) and (Ralph et al. 1996). <br/>All output files stored in a netcdf file.<br/><br/>

**Note1:  More information about the processing of UBWPP can be found on the article:  Garcia-Benadi A, Bech J, Udina M, Campistron B, Paci A. Multiple Characteristics of Precipitation Inferred from Wind Profiler Radar Doppler Spectra. Remote Sens. 2022, 14, 5023. https://doi.org/10.3390/rs14195023**<br/><br/>

## Versions and dependences

The main script is called UBWPP.py and it is available in Python 3.8. or later. The following libraries are necessary:

	numpy, version 1.19.2 or later
	scipy, version 1.5.2. or later 
	astropy, version 4.3.1. or later 	
	netCDF4, version 1.5.4 or later 

	Other necessary libraries are cftime, calendar, datetime, math, os, glob, struct and sys (make sure they are installed or install them otherwise).
	
The libraries can be installed with pip, using these sentences:

	pip install numpy
	pip install scipy
	pip install netCDF4
	pip install astropy
	pip install cftime
	

**As mentioned above the script works with the windprofiler raw data files (.dat).**

## How to cite

If you use this script for your publication, please cite as:<br/>
Garcia-Benadi A, Bech J, Udina M, Campistron B, Paci A. Multiple Characteristics of Precipitation Inferred from Wind Profiler Radar Doppler Spectra. Remote Sens. 2022, 14, 5023. https://doi.org/10.3390/rs14195023


## Outputs
The script produces the following outputs from UBWPP netcdf data in High and Low mode (in the next variable names the term 'mode' will be 'high' or 'low' depending on the mode of the raw data processed):<br />

**W_mode:** Vertical speed (m s<sup>-1</sup>)<br />
**Spectral width_mode:** spectral width of the vertical velocity distribution (m s<sup>-1</sup>)<br />
**Skewness_mode:** skewness of the vertical velocity distribution<br />
**Kurtosis_mode:** kurtosis of the vertical velocity distribution<br />
**Type_mode:** Hydrometeor type (unknown[20], rain [10], mixed [0], snow [-10]) following the Atlas approach.[doi:10.1029/RG011i001p00001]<br />
**Type_2_mode:** Hydrometeor type (unknown[20], rain [10], mixed [0], snow [-10]) following the Ralph approach.[doi:10.1175/1520-0477(1995)076<1717:USMDFN>2.0.CO;2] <br />
**Liquid Water Content_mode:** Liquid water content (g m<sup>-3</sup>)<br />
**RR_mode:** Rain rate (mm h<sup>-1</sup>)<br />
**Z_mode:** Equivalent reflectivity considering only liquid drops (dBZ)<br />
**U_mode:** Horizontal wind component, zonal velocity  (m s<sup>-1</sup>). Available if the device has more than one beam.<br />
**V_mode:** Horizontal wind component, meridional velocity  (m s<sup>-1</sup>). Available if the devices ha more than one beam.<br />
**shape parameter_mode:** shape parameter from the Drop Size Distribution <br />
**slope parameter_mode:** slope parameter from the Drop Size Distribution (m<sup>-1</sup>)<br />
**log10(N0)_mode:** Drop Size Distribution (log10(m<sup>-4-mu</sup>)) supposing that all hydrometeors are in liquid phase<br />
**SNR_mode:** Signal noise relation from signal (dB)<br />
**VKEF_mode:** Vertical kinetic energy (g s<sup>-3</sup>)<br />
**HKEF_mode:** Horizontal kinetic energy (g s<sup>-3</sup>)<br />
**C2n_mode:** refractive index (log(s<sup>(-2/3)</sup>)<br />

**time_utc:** Time vector in UTC format

If the user applies -Pot parameters and additional menu appears:
Signal without processing in three dimensions (Time, Height and bins):
**P1_mode:** Power returned in whole spectra in beam 1 (W)<br />
**P2_mode:** Power returned in whole spectra in beam 2 (W), if the beam is avalaible<br />
**P3_mode:** Power returned in whole spectra in beam 3 (W), if the beam is avalaible<br />
**P4_mode:** Power returned in whole spectra in beam 4 (W), if the beam is avalaible<br />
**P5_mode:** Power returned in whole spectra in beam 5 (W), if the beam is avalaible<br />


<br />


## How to execute the script
The script can be executed from a command line at the system prompt (see MS-Windows example):<br />
<br />
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)
<br />
at the directory where UBWPP.py has been copied:
```
python UBWPP.py

```
The script has some additional command line execution options. Please note that their use implies a substantial increase of the netcdf output file (see below). <br />The possible command line arguments available are (more than one is possible, in any order): <i>-Pot</i>  <i>-xxx</i>  <i>-hxxxx</i>.<br /> 
<i>-Pot</i>: the script saves the values of the Power for each radar beam, in high and low mode. The netcdf size increases about 4 times.<br />
<i>-NameFile</i>:The script asks the user the name of the netcdf file output, if this option isn't chosen the file name will be the year, month and day in format yyyymmdd.nc<br />
<i>-hxxxx</i>: forces the antenna height is at xxx meters above sea level.  (for example -100.3 would mean the antenna is at 100.3 m above sea level). This parameter is important to determine the hydrometeor terminal speed in function of height.<br />
<i>-cHxxxx</i>: Gives the calibration constant in high mode, if it isn't included the value used will be 116.5 dB.<br />
<i>-cLxxxx</i>: Gives the calibration constant in low mode, if it isn't included the value used will be 111.3 dB.<br />

The syntax of these options are:

```
python UBWPP.py -Pot

```
```
python UBWPP.py -NameFile

```
```
python UBWPP.py -Pot -NameFile

```
```
python UBWPP.py -Pot -NameFile -h100.3

```
```
python UBWPP.py -Pot -NameFile -h100.3 -cH120.5 -cL105.8

```

The script asks the directory where the dat files to be processed are located (it will process all the RWP dat files of the selected directory), for example:
```
C:\RWP\data\
```
**NOTE 1: The path must end with \\ in Windows or a / in Linux**<br />
**NOTE 2:  Please avoid blank spaces and special characters in your file path**<br />
**NOTE 3: In macOS, depending on your path environment configuration, it may not be necessary to indicate the complete path so "./" may be enough**<br />

The script indicates the number of dat files in directory, the UHF location (longitude and latitude), the number of modes, the radar frequency and its wavelenght, and the number of beams. The next line indicates the processing time estimation, and the next line shows the % of file processed.

The result is stored in a netcdf file.


## Contact
If you have any question, please contact with Albert at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
