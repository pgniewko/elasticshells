>Notice: This is research code that will not necessarily be maintained in the future.
>The code is under development so make sure you are using the most recent version.
>We welcome bug reports and PRs but make no guarantees about fixes or responses.

USAGE
=====
* ```observers.config``` - configures observers that are gauged during the simulation or that will be executed by the ```analyze.sh``` script.
Each observer is configured in a separate line. For details on configuring the observers, please refer to the original source-code or contact Pawel Gniewek at gniewko.pablo@gmail.com 

* ```schedule.config``` - configures the simulation protocol. Each line in this file creates a separate schedule. 
The columns in a schedule line refer to:
    1. Number of steps to be executed with a given schedule
    2. Time interval at which the schedule is executed
    3. Average change ```dx``` of the X-dimension of the box
    4. Average change ```dy``` of the Y-dimension of the box
    5. Average change ```dz``` of the Z-dimension of the box
    6. Range of the ```dx``` size change randomness. The value should be gt 0. 0 if ```dx``` is meant to always be the same magnitude
    7. Range of the ```dy``` size change randomness. The value should be gt 0. 0 if ```dy``` is meant to always be the same magnitude
    8. Range of the ```dz``` size change randomness. The value should be gt 0. 0 if ```dz``` is meant to always be the same magnitude
    9. The upper limit of the volume fraction - above which the schedule is no longer executed.

COPYRIGHT NOTICE
================
Copyright (C) 2014-2018, Pawel Gniewek  
Email : gniewko.pablo@gmail.com  
All rights reserved.  
License: BSD 3  
