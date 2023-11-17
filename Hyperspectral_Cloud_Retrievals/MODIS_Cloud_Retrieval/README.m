----------------------------------------------
------------ MODIS-Cloud-Retrieval------------
----------------------------------------------

The set of functions in this folder will preform MODIS-like operations on Level 0/1 data products in order to calculate cloud properties such as cloud optical thickness, effective cloud droplet radius and thermodynamic phase.


INTRO
=====
The MODIS algorithm assumes a plane-parallel, liquid water cloud in all of its retrievals.  The DISORT RTE solver automatically uses a plane-parallel atmosphere. To define water clouds in uv_spec, you must include the data file 'wc_file' (for ice clouds, 'ic_file'). 

MODIS uses a scan mirror to take images in the across track direction (perpendicular to the direction of motion). 


DESCRIPTION OF MODIS DATA AND VARIABLES
=======================================

There are several data sets of interest to this code. An incomplete list is described below
	- Level 0: Raw instrument packets
	- Level 1A: scans of raw radiances in counts
	- Level 1B: Calibrated radiances
	- Geolocation: Set of lat/longs for each pixel on the ground and the Sun's position


Each on of these data have a different HDF file and variable suite. Below is an incomplete list of the relevant variables for each Level.

	
	- Level 1A HDF output:
		- SDS (Science Data Sets): multidimensional arrays used to store science data
			- EV_250m: Earth-viewing data at 250 meter resolution. Data stored as (track, band, frame), where track is the number number of scans times the number of detectors, band gives the data at some particular MODIS band, and frame is the number of frames times the number of samples. This SDS covers spectral bands 1 and 2
			- EV_500m: Earth-viewing data at 500 meter resolution. Data stored as (track, band, frame). This SDS covers spectral bands 3,4,5,6,7
			- EV_1km_day: Earth-viewing data at 1 kilometer resolution. Data stored as (track, band, frame). This SDS covers spectral bands 8-19, and 26
			- EV_1km_night: Earth-viewing data at 1 kilometer resolution. Data stored as (track, band, frame). This SDS covers spectral bands 20-25, and 27-36
		- vdata (fixed length data)
		- file attributes
	- Level 1B HDF output:

		- VGROUP: 
			

	- GeoLocation Data: Contains the values for sensor azimuth and zenith, solar azimuth and zenith, sensor height and range with respect to each pixel on the ground

		- Solar Azimuth: 
		- Height: Height of the sensor above the Earth ellipsoid
		- Range: range to the satellite from... 


