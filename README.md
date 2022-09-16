# TDOA tracking package

A Matlab package for multi-target tracking of Time-Differences-Of-Arrivals (TDOAs) of signals between two sensors from towed hydrophone array recordings.

For details see Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.

Copyright (c) 2021, Pina Gruden


## How to use

Before running the package specify the array, paths and parameters for your application by modifying the following scripts: 
- Array_info.csv - If needed add the array information in a new row, specifying sensor separations for your array. Note- the name you give your array should then match the one you specify in the Specify_Parameters4Xcorr.m.
- Specify_Parameters4Xcorr.m - This is where you specify your array, encounter information and settings for the cross-correlogram. Change any parameters in the sections labeled “ CHANGABLE:” as needed. 
- Specify_Parameters4Tracking.m - This is where you specify parameters for measurement extraction and tracking. Change any parameters in the sections labeled “ CHANGABLE:” as needed. 
- Specify_Paths.m - Specify folders where data is located and where results should be saved to. The expected data format are .wav files, and you can process either the full encounter or individual files. It is expected that you will be processing one encounter at the time (or if individual files are processed, that these files are from the same encounter). The expected name for .wav files is 'xxx_yyyyMMdd_HHmmss_SSS.wav', where 'xxx_' can be any string (or none), and 'yyyyMMdd_HHmmss_SSS' part of the name specifies the date (year, month, day) and time (hours, minutes, seconds, and milliseconds (SSS)).


Then run the package by running:
1) A1_Compute_CrossCorrelograms.m - This computes and saves the cross-correlogram based on your audio data. IMPORTANT: specify what signal type you want to be processing for (line 16)- choose either "clicks", "whistles", or "both", depending on your application.
2) A2_ExtractMeasurements_and_TrackTDOAs.m - This extracts measurements and tracks TDOAs based on that and returns the extracted TDOA tracks. Results are also displayed as plots.


## Output

The package outputs cross-correlograms, extracted TDOA tracks, and plots the results. 

Specifically, the outputs that are saved in a .mat file are:
1) From A1_Compute_CrossCorrelograms.m part of the processing:
- Cross-correlograms ('Rxy_envelope_ALL') (based on envelopes of the generalized cross-correlation function)- matrix with a dimension of NxM, where N is the number of time steps and M is the number of TDOAs. 
- Range of possible TDOAs ('lags') - a vector of possible TDOAs for a given sensor separation.
- Time vector ('t_serialdate') - a time vector in datetime format (for more info type help datetime in Matlab's command prompt)
- Time vector ('t') - a time vector in seconds from the beginning of the file/encounter.
- Parameters used in the processing ('parameters')- these are parameters specifying array & encounter information, parameters used in cross-correlogram computation, measurement extraction.
- Parameters used specifically in processing of clicks and/or whistles ('param_signal')

2) From A2_ExtractMeasurements_and_TrackTDOAs.m part of the processing:
- Extracted TDOA tracks ('Tracks') - 1xM structure, where M is number of targets. For each target there are three fields: 'time', 'time_local','tdoa', where 'time' refers to time is seconds from the start of the file/encounter, 'time_local' refers to time in timedate format.
- Models used in GM-PHD-SA filter ('model')- a structure array containing models and parameters for the filter.
- Parameters used for normalization, measurement extraction and track extraction ('parameters_measure_tracking')
- Parameters used in the processing ('parameters')- these are parameters specifying array & encounter information, parameters used in cross-correlogram computation, measurement extraction. (Comes from A1 part of the processing)
- Parameters used specifically in processing of clicks and/or whistles ('parameters_clicks' and 'parameters_whistles' or 'param_signal') (Comes from A1 part of the processing)
- Scalars used for normalization of the cross-correlograms ('scalar_clicks'  and 'scalar_whistles' or 'scalar').






 
