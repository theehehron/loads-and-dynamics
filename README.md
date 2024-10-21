# Loads And Dynamics
## Repository of MATLAB and Python Code Related to Loads and Dynamics

### Purpose
This code was written as I worked through Tom Irvine's tutorials from [vibrationdata.com](http://www.vibrationdata.com). I started this repository with MATLAB, but I will primarily be adding to the python directory from now on (to sharpen my python skills).

### SDOF reponse to base input
The code can be used to calculate/plot the response of an SDOF system to a base input. An example is shown below (this plot was created using the MATLAB code):

![SDOF_response MATLAB](sdof_responses.png)

You can also calculate/plot the vibration response spectrum for a base vibration input as shown below (plot created in python):

![vrs python](python/vrs_python.png)

The "equivalent" shock response spectrum to a random vibration test can be calulated given the random vibration PSD and a test duration. The "equivalent" SRS is assumed to be equal to the n$\sigma$ vibration response spectrum, where n $= \sqrt2\ln{(f_nT)}$ where $f_n$ is the natural frequency and $T$ is the test duration. The plot below shows the equivalent SRS for MIL-STD-1540C ATP with test duration of 180 seconds.

![vrs srs equivalence](python/vrs_nsigma_python.png)