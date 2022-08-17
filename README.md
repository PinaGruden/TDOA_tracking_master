# TDOA_tracking_master
Package for tracking Time-Difference-Of-Arrivals from two sensors

TDOA_tracking_master is a Matlab package for tracking Time-Differences-Of-Arrivals (TDOAs) based on two sensors from towed hydrophone array.

For details see Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.


Before running the package specify the paths and parameters for your application by modifying the: 
- Specify_Parameters.m - if needed change any parameters in the sections labeled “ CHANGABLE:” 
- Specify_Paths.m - specify folders where data is located and where results should be saved to.


Then run the package by running either:
- RUN_TDOA_TRACKING_mixedsignaltypes.m - this is the case where the species produces clicks and whistles and you want to track based on both of these signals.
- RUN_TDOA_TRACKING_onesignaltype.m - this is the case where the species produces a single signal type.  
