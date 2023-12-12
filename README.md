# TDOA tracking package

This is a Matlab based package that performs tracking of Time-Difference-Of-Arrivals (TDOAs) from multiple sources based on information from two moving sensors in a linear array. It was written and tested for towed arrays, but could also be used with stationary arrays.

The TDOAs are obtained from cross-correlograms, which are computed based on either clicks, whistles, or both. In addition to TDOAs, amplitudes are also obtained from the cross-correlograms. Both TDOAs and amplitudes represent measurements based on which the tracking is carried out. The tracking uses GM-PHD-SA filter (Gruden et al, 2021) to track (i.e., connect all TDOAs from one source or closely spaced group of sources into a continuous trajectory over time). See Gruden et al (2021) for more details on the cross-correlogram computation, extraction of measurements, and tracking.

_References_:

Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.

Copyright (c) 2021, Pina Gruden


## 1.   Required Matlab version and toolboxes

This package was developed and tested with **Matlab version 2022a (9.12)**. It uses the following Matlab toolboxes:
- *Signal Processing Toolbox*
- *Statistics and Machine Learning Toolbox*
- *Curve Fitting Toolbox*

## 2.  Package contents
Package contains folders and functions that execute pre-processing, measurement extraction, and TDOA tracking.

The package contains the following functions, scripts, and files in the main folder:

- `A1_Compute_CrossCorrelograms.m` - main script that computes cross-correlograms.
- `A2_ExtractMeasurements_and_TrackTDOAs.m` - main script that extracts measurements from cross-correlograms and extracts TDOA tracks from them.
- `Array_Info.csv` - a table specifying array spacing.
- `plot_results.m` - plots TDOA tracking package results.
- `README.md` - readme document specifying package usage.
- `Specify_Parameters4Tracking.m` - specifies parameters for measurement extraction and tracking.
- `Specify_Parameters4Xcorr.m`  - specifies parameters for cross-correlogram computation.
- `Specify_Paths.m`  - specifies paths to folders where data is located and results saved to.

The package contains the following folders in the main folder:

1) **./Preprocess_Extract_Measurements/** - contains code to pre-process audio files, to construct cross-correlograms, and to extract measurements for tracking. Functions included are:

- `compute_crosscorrelogram.m` - computes a cross-correlogram between two channels in audio data.
- `extract_measure_crosscorrelogram.m` -  extracts TDOAs and amplitudes of the cross-correlation information from a cross-correlogram and saves it as a random finite set (RFS).
- `extract_peaks_crosscorr.m`- extracts peaks from the cross-correlogram.
- `gcc.m` - computes the generalized cross-correlation (GCC) of two signals. 
- `getStrDateTime.m` - gets date and time from a file header.
- `norm_background_crosscorr.m` - normalizes background noise to have std=1 under the Rayleigh assumption.
- `preprocess.m` - applies de-clicking of the signal. 

2) **./GMPHD_SA/** - contains code to extract TDOA tracks from cross-correlograms using GM-PHD-SA filter. Functions and files included are:

- `birth_process_amplitude.m` - generates newborn targets based on the measurements.
- `gen_tdoa_tracking_models.m` - generates models required 
           for TDOA tracking with GMPHD-SA filter.
- `gmphd_adaptive_amplitude.m` - tracks multiple targets
         using the Gaussian Mixture Probability Hypothesis Density- SA filter.
- `gmphd_limit.m` - eliminates Gaussian components based on their
        weights.
- `gmphd_merge.m` - merges Gaussian components that are close together.
- `gmphd_prune.m` - performs pruning of Gaussian components.
- `kalman_predict_multiple.m` - uses Kalman filter to predict target states according to system model.
- `kalman_update_multiple.m` - uses Kalman filter to update target states according to measurement model.
- `kalman_update_single.m` - uses Kalman filter to update a single target state according to measurement model.
- `postprocess.m` - removes short tracks and smooths the tracks.
- `tracktarget_tdoa_labels.m` - finds targets with the same label and collects them into tracks.
- `BayesOptimization_GMPHDParams_GTchunked_lambda4_5_DiffBirth.mat` - optimized parameters.
- `BirthVelocityPrior_AllTrainData.mat` - birth velocity prior from training data.

3) **./Test_example/** - contains data to use for testing the package works correctly.



## 3.  How to use

In order to use this package you must first specify the array used to record the measurements and parameters that you want to use for tracking - see Section  3.1. After this, you can run the measurements extraction and tracking - see Section  3.2. 

###  3.1. Modify

First, specify paths, the array, and parameters for your application by modifying the following scripts: 
- `Specify_Paths.m` - Specify folders where data is located and where results should be saved to. The expected data format are .wav files, and you can process either the full encounter or individual files. It is expected that you will be processing one encounter at the time (or if individual files are processed, that these files are from the same encounter). The expected name for .wav files is 'xxx_yyyyMMdd_HHmmss_SSS.wav', where 'xxx_' can be any string (or none), and 'yyyyMMdd_HHmmss_SSS' part of the name specifies the date (year, month, day) and time (hours, minutes, seconds, and milliseconds (SSS)).
- `Array_info.csv` This is where your hydrophone array sensor separation is specified. If adding a new array, enter the array information and sensor separation in a new row. Note, the name you give your array should then match the one you specify in the `Specify_Parameters4Xcorr.m`.

   **Important:** The hydrophone sensor spacing should be specified in a continuous manner in meters from the first to the last hydrophone. For example: Assume your array consists of two sub-arrays that are separated by a 20 m separator, and each sub-array has three sensors that are separated by 1 m. Further assume that your first sensor in each sub-array is placed right at the very beginning of a given sub-array, and assume the last sensor in each sub-array is placed right at the very end of a given sub-array. Then you enter the following values for hydrophones (each in a corresponding column): 0, 1, 2, 22, 23, 24. In practice, there is typically some separation between the beginning of the sub-array and the first sensor, and between the last sensor and the end of the sub-array (so measure and note sensor positions in a continuous manner). You will also have to enter the distance at which the array is towed behind the boat- this is a distance of your tow cable- from the back of the boat to the beginning of the first sub-array (or array if you only have one section with sensors). Carefully measure and note all sensor positions in a continuous manner. Insert -999 for fields where spacing is unknown (but the channel exists in the recordings), or if your number of sensors is less than 6. 
- `Specify_Parameters4Xcorr.m` - This is where you specify your array, encounter information, and settings for the cross-correlogram. Scroll down to change any parameters in the sections labeled "CHANGABLE:" as needed.  
- `Specify_Parameters4Tracking.m` - This is where you specify parameters for measurement extraction and tracking. Scroll down to change any parameters in the sections labeled "CHANGABLE:" as needed. 

###  3.2. Run

Then run the package by running:
1) `A1_Compute_CrossCorrelograms.m` - This computes and saves the cross-correlogram based on your audio data. IMPORTANT: specify what signal type you want to be processing for (line 16)- choose either "clicks", "whistles", or "both", depending on your application.
2) `A2_ExtractMeasurements_and_TrackTDOAs.m` - This extracts measurements and tracks TDOAs based on that and returns the extracted TDOA tracks. Results are also displayed as plots.


## 4. Output

The package outputs cross-correlograms, extracted TDOA tracks, plots, and saves the results.

Specifically, the outputs that are saved in a .mat file are:
1) From `A1_Compute_CrossCorrelograms.m`:
- Cross-correlograms ('Rxy_envelope_ALL') (based on envelopes of the generalized cross-correlation function)- matrix with a dimension of NxM, where N is the number of time steps and M is the number of TDOAs. 
- Range of possible TDOAs ('lags') - a vector of possible TDOAs for a given sensor separation.
- Time vector ('t_serialdate') - a time vector in datetime format (for more info type help datetime in Matlab's command prompt)
- Time vector ('t') - a time vector in seconds from the beginning of the file/encounter.
- Parameters used in the processing ('parameters')- these are parameters specifying array & encounter information, and parameters used in cross-correlogram computation.
- Parameters used specifically in processing of clicks and/or whistles ('param_signal')

2) From `A2_ExtractMeasurements_and_TrackTDOAs.m`:
- Extracted TDOA tracks ('Tracks') - 1xM structure, where M is number of targets. For each target there are three fields: 'time', 'time_local','tdoa', where 'time' refers to time is seconds from the start of the file/encounter, 'time_local' refers to time in timedate format.
- Models used in GM-PHD-SA filter ('model')- a structure array containing models and parameters for the filter.
- Parameters used for normalization, measurement extraction, and track extraction ('parameters_measure_tracking').
- Parameters used in the cross-correlogram computation ('parameters'), (Comes from processing A1).
- Parameters used specifically in processing of clicks and/or whistles ('parameters_clicks' and 'parameters_whistles' or 'param_signal'), (Comes from processing A1).
- Scalars used for normalization of the cross-correlograms ('scalar_clicks' and 'scalar_whistles' or 'scalar').


