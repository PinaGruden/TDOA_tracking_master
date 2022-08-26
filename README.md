# TDOA tracking package

A Matlab package for multi-target tracking of Time-Differences-Of-Arrivals (TDOAs) of signals between two sensors from towed hydrophone array recordings.

For details see Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.

Copyright (c) 2021, Pina Gruden


## How to use

Before running the package specify the array, paths and parameters for your application by modifying the: 
- Array_info.csv - if needed add the array information in a new row, specifying sensor separations for your array. Note- the name you give your array should then match the one you specify in the Specify_Parameters.m.
- Specify_Parameters.m - if needed change any parameters in the sections labeled “ CHANGABLE:” 
- Specify_Paths.m - specify folders where data is located and where results should be saved to. The expected data format are .wav files, and you can process either the full encounter or individual files. The expected name for .wav files is 'xxx_yyyyMMdd_HHmmss_SSS.wav', where 'xxx_' can be any string (or none), and 'yyyyMMdd_HHmmss_SSS' part of the name specifies the date (year, month, day) and time (hours, minutes, seconds, and milliseconds (SSS)).


Then run the package by running either:
- RUN_TDOA_TRACKING_mixedsignaltypes.m - this is the case where the species produces clicks and whistles and you want to track based on both of these signals.
- RUN_TDOA_TRACKING_onesignaltype.m - this is the case where the species produces a single signal type.  


## Output

The package outputs extracted TDOA tracks, as well as measurements and cross-correlograms, and plots the results.

Specifically, the outputs that are saved in a .mat file are:
- Normalized cross-correlograms ('Rxy_envelope_scaled_clicks' and/or 'Rxy_envelope_scaled_whistles') - matrix with a dimension of NxM, where N is the number of time steps and M is the number of TDOAs. Additional outputs are also the scalars used for normalization of the cross-correlograms ('scalar_clicks' and/or 'scalar_whistles').
- Range of possible TDOAs ('lags') - a vector of possible TDOAs for a given sensor separation.
- Time vector ('t_serialdate') - a time vector in datetime format (for more info type help datetime in Matlab's command prompt)
- Time vector ('t') - a time vector in seconds from the beginning of the file/encounter.
- Parameters used in the processing ('parameters')- these are parameters specifying array & encounter information, parameters used in cross-correlogram computation, measurement extraction.
- Parameters used specifically in processing of clicks and/or whistles ('parameters_clicks','parameters_whistles')
- Measurements on which the tracking was performed ('measure')- a structure containing measurements ('measure.Z') and number of time steps ('measure.T'). 'measure.Z' is a 1xN cell array, where N is number of time steps, and each cell contains measurements consisting of TDOAs and amplitudes of the cross-correlation information for that time step.
- Models used in GM-PHD-SA filter ('model')- a structure array containing models and parameters for the filter.
- Extracted TDOA tracks ('Tracks') - 1xM structure, where M is number of targets. For each target there are three fields: 'time', 'time_local','tdoa', where 'time' refers to time is seconds from the start of the file/encounter, 'time_local' refers to time in timedate format.



  