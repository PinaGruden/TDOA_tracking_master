{\rtf1\ansi\ansicpg1252\cocoartf2580
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 TDOA_tracking_master is a Matlab package for tracking Time-Differences-Of-Arrivals (TDOAs) based on two sensors from towed hydrophone array.\
\
For details see Gruden, P.,  Nosal, E.-M. and Oleson, E. (2021). Tracking time differences of arrivals of multiple sound sources in the presence of clutter and missed detections. The Journal of the Acoustical Society of America  150(5): 3399--3416.\
\
\
Before running the package specify the paths and parameters for your application by modifying the: \
- Specify_Parameters.m - if needed change any parameters in the sections labeled \'93 CHANGABLE:\'94 \
- Specify_Paths.m - specify folders where data is located and where results should be saved to.\
\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf0 Then run the package by running either:\
- RUN_TDOA_TRACKING_mixedsignaltypes.m - this is the case where the species produces clicks and whistles and you want to track based on both of these signals.\
- RUN_TDOA_TRACKING_onesignaltype.m - this is the case where the species produces a single signal type.  }