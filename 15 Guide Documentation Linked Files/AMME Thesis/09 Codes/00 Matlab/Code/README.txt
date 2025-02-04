README

XPlane - MATLAB interactions
Jeremy Cox
Brendan Waters

Updated 29/04/2023
-'udp' deprecated, superceded by 'udpport'
- Introduced function to calculate and store geodetic position from state-space

To Do:
- Override displyed control inputs
- Override cockpit indicators

There are three main modes of operation that have been assessed or developed:
-Driving XPlane as a visual display
    -From data within a MATLAB script
    -From Simulink Desktop Real-Time
-Getting flight data recording out of XPlane
	- NOT RECOMMENDED receiving flight data in a MATLAB script in real time

REQUIREMENTS:
-MATLAB R2021a or later recommended (for ned2lla)
-Simulink Desktop Real-Time and dependencies
-Instrument Control Toolbox and dependencies
-Potentially other toolboxes

===========================================================================
Setup notes
MATLAB
1) 	Create a 3 x nSteps matrix to store [lat (rad) ; lon (rad) ; alt (m)] each timestep
2) 	Create a 3 x nSteps matrix to store [phi (rad) ; theta (rad) ; heading (rad)] each timestep.  Heading should be true heading but psi works
3) 	In the timeloop call the function computeGeodeticPosition to get integrated Lat, Lon, Alt, i.e:
	
		lla(:,i) = computeGeodeticPosition(u,v,w,phi,theta,psi,dt,lla_o);

		Where:
		      - u, x-axis velocity [m/s] (body axis)
		      - v, y-axis velocity [m/s] (body axis)
      	  	      - w, z-axis velocuty [m/s] (body axis)

		      - phi, roll angle [rad]
  		      - theta, pitch angle [rad]
		      - psi, yaw angle [rad]

 	 	     - dt, temporal step [s]

 	   	     - lla_o, 3x1 vector containing [lat; lon; alt] from previous timestep i.e. lla(:,i-1)

4) 	Store pth(:,i) = X(7:9);  (should be true heading not psi, psi used here as example)
5) 	Outside the timeloop, convert lla and ptp matricies from radians to degrees
6) 	Send UDP packets to xplane using xplane_Visualise(t,lla,ptp);


X-Plane 11
1)	Start Up X-Plane
2)	Load Flight:
	(a) Choose a start location near the commanded [lat; lon; alt] (from matlab).
		i.e. if matlab is starting from  33°52′11.44″ South, 151°12′29.83 East (Sydney)
		     then select the Sydney Airport as starting location
	(b) Set time of day/weather to appropriate selections (dark + gloomy is hard to see)
	(c) Set view point 
		- SHIFT + 8 is chase view, can use ',' key '.' key and arrows to adjust position
		- SHIFT + 7 is the pilot POV
3) 	Start Flight 
4)	Run MATLAB

===========================================================================
Operation notes

Driving XPlane is done using a range of UDP packets. There are a multitude of different packets that can be sent, as in the documentation for XPlane 
(file TXT.rtf in the folder ~\X-Plane 11\Instructions\X-Plane SPECS from Austin\Exchanging Data with X-Plane.rtfd). For our purposes, VEHX, DREF and 
DATA have been implemented. DREF should in theory allow nearly any interaction with the aircraft, and setting nearly any of the aircraft states 
(including the state of cockpit gauges, switches and controls).

VEHX packets will stop XPlane's built-in flight model from running. The simplest way to restart this is to restart XPlane.

When sending VEHX packets, XPlane can sometimes glitch out and fail to respond immediately, or try to draw the aircraft in some kind of limbo if you suddenly cease sending VEHX packets.

MATLAB's 'udp' function has been used, and is set to be deprecated (I started scripting in a version of MATLAB where it was not yet flagged as to be deprecated). 
It should not be difficult to migrate to the newer udp function 'udpport', as the functionality is essentially unchanged.

===========================================================================
Driving XPlane as a visual display from within a MATLAB script

xplane_Visualise_test.m shows a minimum working example of this.

A recorded history of aircraft states can be displayed in XPlane. To do this, the aircraft states must be formatted as times, lat/lon/altitude and phi/theta/psi in 3xn matrices. 

The aircraft states must be ordered chronologically with increasing time. Very large time histories may reduce performance (this could be overcome but is not an issue for practical 
size flight state histories). Where aircraft states are provided at very high rates, not all states may be displayed due to graphics framerate, and MATLAB will skip sending aircraft 
states if it cannot send them quickly enough (MATLAB will *probably* be able to send the states faster than your graphics framerate is capable of rendering them).

xplane_UDP_Demo.m shows how to send VEHX, DREF and DATA packets manually.

===========================================================================
Driving XPlane as a visual display from Simulink Desktop Real-Time

Display_Demo.slx shows a minimum working example of this. You will need to run DREF_init.m first which sets up some large arrays for this to run properly.

The VEHX packet writes the position and attitude of the aircraft to XPlane.

The DREF packet attempts to write to the airspeed indicator in the cockpit. Depending on the way a given aircraft has been setup, writing to the gauges in the virtual cockpit may work 
differently, or be overriden by behaviours related to the implementation of that cockpit.

In the example of the Cessna Skyhawk, a lua script is used by this aircraft to update the indicated airspeed gauge. This script needs to be hamstrung by commenting out the parts where 
it writes these values, otherwise the lua script overwrites any IAS that is sent by simulink real time. It's likely that connecting other gauges so they work properly involves similar
work to figure out which DREF drives each gauge, and figuring out if there is any aircraft-specific way that the gauge is implemented that prevents the gauge from being written.

===========================================================================
Getting flight data recording out of XPlane

To do this, run XPlane and go to data output settings, general data output. Select the checkboxes in the Disk (data.txt File) column that you require. It's handy to show them in the 
cockpit so you can see what you are getting. Adjust the desired data rate with the slider on the right. Adjusting to 'max' 99.9 samples/second seems to instead write once per graphics 
frame, even if graphics frames are rendered faster than 99.9/sec (physics sub-steps per graphics frame can be configured, and is by default set to 4).

data.txt is in the XPlane directory. data.txt is overwritten when XPlane is started, so previous sessions will be deleted unless you save a copy. data.txt seems to only get populated
 when XPlane is closed.

MATLAB's import interface readily interprets these files.

===========================================================================
NOT RECOMMENDED Getting flight data into MATLAB in real time

Data_Logger.m will receive data from XPlane in real time. MATLAB requests a stream of data. RPOS data gives aircraft flight state (spatial position and velocity, angular position and velocity). 

Other information (like control surface states) can be sent by checking 'Network via UDP' in the same interface as logging to data.txt. A suitable network configuration needs to be setup for 
this in the same page of settings. It's possible (but untested) to configure which DATA packets are sent by sending instructions over UDP from MATLAB (essentially allowing programmatic setup 
of which data is sent and how fast).

This is NOT RECOMMENDED because matlab was not able to receive these packets at the full rate they were sent, and the packet loss rate varied from none to sporadicly dropping many packets in a
row. Because of this limitation, the script is incomplete, packets are not parsed (only check for length to find RPOS packets, DATA packets would need a switch statement to parse them correctly 
depending on what data has been sent). It's possible the newer MATLAB function udpport will alleviate the dropped packet issue so that packets can be received at full rate.