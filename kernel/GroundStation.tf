KPL/FK
 
   FILE: GroundStation.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.2.0 --- September 6, 2016
   PINPOINT RUN DATE/TIME:    2019-04-23T21:35:36
   PINPOINT DEFINITIONS FILE: def.txt
   PINPOINT PCK FILE:         pck00010.tpc
   PINPOINT SPK FILE:         GroundStation.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'GROUNDSTATION'
   NAIF_BODY_CODE                      += 399123
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame GROUNDSTATION_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame GROUNDSTATION_TOPO is centered at the
      site GROUNDSTATION, which has Cartesian coordinates
 
         X (km):                 -0.5018181404560E+04
         Y (km):                  0.3532655423430E+04
         Z (km):                 -0.1731653379930E+04
 
      and planetodetic coordinates
 
         Longitude (deg):       144.8555559999836
         Latitude  (deg):       -15.8583330674403
         Altitude   (km):         0.4010634573776E-03
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  6.3781366000000E+03
         Polar radius      (km):  6.3567519000000E+03
 
      All of the above coordinates are relative to the frame ITRF93.
 
 
\begindata
 
   FRAME_GROUNDSTATION_TOPO            =  1399123
   FRAME_1399123_NAME                  =  'GROUNDSTATION_TOPO'
   FRAME_1399123_CLASS                 =  4
   FRAME_1399123_CLASS_ID              =  1399123
   FRAME_1399123_CENTER                =  399123
 
   OBJECT_399123_FRAME                 =  'GROUNDSTATION_TOPO'
 
   TKFRAME_1399123_RELATIVE            =  'ITRF93'
   TKFRAME_1399123_SPEC                =  'ANGLES'
   TKFRAME_1399123_UNITS               =  'DEGREES'
   TKFRAME_1399123_AXES                =  ( 3, 2, 3 )
   TKFRAME_1399123_ANGLES              =  ( -144.8555559999836,
                                            -105.8583330674403,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file def.txt
--------------------------------------------------------------------------------
 
begindata
 
         SITES         = ( 'GROUNDSTATION')
 
         GROUNDSTATION_CENTER = 399
         GROUNDSTATION_FRAME  = 'ITRF93'
         GROUNDSTATION_IDCODE = 399123
         GROUNDSTATION_XYZ    = ( -5018.18140456, 3532.65542343, -1731.65337993 )
         GROUNDSTATION_UP     = 'Z'
         GROUNDSTATION_NORTH  = 'X'
 
 
begintext
 
[End of definitions file]
 
