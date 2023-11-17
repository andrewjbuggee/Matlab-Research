-----------------------------------------------------------------------------
------------------- How to Read AVIRIS Data files ---------------------------
-----------------------------------------------------------------------------
-----------------------------------------------------------------------------
------------------------- By Andrew J. Buggee -------------------------------
-----------------------------------------------------------------------------


INTRO:
------
	AVIRIS classic stores data files in the ENVI binary image format. Each image is a binary image cube with a detached ASCII header. There are 3 levels of data. Level 0: raw data, Level 1: calibrated radiances, Level 2: atmospherically corrected surface radiances. The level 1 calibrated data us usually orthorectified, meaning the data is mapped to an (x,y) coordinate system. This mapping is stored in the geometry file .IGM


AVIRIS Data Products:
---------------------
	Level 0 (raw data):
	-------------------


	Level 1 (calibrated radiance): This is orthorectified, calibrated, at-sensor radiance data.
	------------------------------
	It is usually band interleave by pixel. This is almost always the data product we want to
	work with. 

	Level 2 (reflectance):
	----------------------



HOW TO READ DATA USING MATLAB: 
------------------------------

	Method 1 - multibandread():  This built-in function is designed to read hyper spectral data
	---------------------------
	cubes of ENVI format. This is also known as band-interleaved data stored as a binary file.   
	Band here means spectral bands. The inputs are as follows:
	[dataCube] = multibandread(filename,dataDim,dataType,headerOffset,interleave,byteOrder)
	
		filename: a bit self explanatory. you need the full folder and file path.
		---------

		dataDim: data dimensions are found in the ort_img.hdr (header) file. The inputs are
		---------
		[lines,samples,bands], which equates to [rows,columns,wavelenghts]

		dataType: The number precision used. The most common type is 'int16' or 16 bit 
		---------
		signed integer, for the AVIRIS data. The other options can be seen in the 
		read_AVIRIS_classic.m file

		headerOffset: The number of bytes that make up the header of a given file. Often the 	
		-------------
		header offset is 0, and the header information is simply stored in a separate file.
		This information is stored in the ort_img.hdr file

		interleave: This tells matlab how the band data is distributed throughout the data	
		-----------
		cube. Band data can be band-sequential ('bbq') where each line is followed by a line
		of the next band. There is also band-interleaved-by-pixel, ('bip') where the first 
		pixel at every band is stored, followed by the second pixel at every band. Finally
		there is band-interleaved-by-line ('bil') which stores the first band in the first
		line, and the second band in the second line, etc.


		byteOrder: The order of bytes. This isn't a helpful description, but it will have
		----------
		to do. This is some sort of numeric formatting thing. There can be 'native' or 
		'host' formatting, or, more commonly, there can be 'ieee-be' or 'ieee-le', which
		is known as a network format. 

