# TDOA tracking package

A Matlab package for multi-target tracking of Time-Differences-Of-Arrivals (TDOAs) of signals between two sensors from towed hydrophone array recordings.

For details see Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.

Copyright (c) 2021, Pina Gruden


## How to use

Before running the package specify the paths and parameters for your application by modifying the: 
- Specify_Parameters.m - if needed change any parameters in the sections labeled “ CHANGABLE:” 
- Specify_Paths.m - specify folders where data is located and where results should be saved to. The expected data format are .wav files, and you can process either the full encounter or individual files. 


Then run the package by running either:
- RUN_TDOA_TRACKING_mixedsignaltypes.m - this is the case where the species produces clicks and whistles and you want to track based on both of these signals.
- RUN_TDOA_TRACKING_onesignaltype.m - this is the case where the species produces a single signal type.  


## Output

The package outputs extracted TDOA tracks, as well as measurements and cross-correlograms, and plots the results. Specifically, the outputs are:
- Normalized cross-correlograms ('Rxy_envelope_scaled_clicks' and/or 'Rxy_envelope_scaled_whistles') - matrix with a dimension of NxM, where N is the number of time steps and M is the number of TDOAs
  