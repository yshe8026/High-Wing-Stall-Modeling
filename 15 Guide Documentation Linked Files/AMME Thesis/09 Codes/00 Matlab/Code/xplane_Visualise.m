function xplane_Visualise(t,lla,ptp)
%Visualise a lat/lon/alt and phi/theta/psi flight history in XPlane.
%t - 1D vector of length n containing the timestamps for the lla and ptp
%data.t must be monotonically increasing, but does not need to be uniformly
%spaced. If the data sample rate is higher than MATLAB can send to XPlane, 
%or higher than x plane's framerate then some data points may be skipped.
%lla - 3xn matrix containing lat/lon/alt. altitude is positive up, lat/lon
%are in degrees. Can be found from NED Cartesian coordinates using ned2lla
%(requires 2021a or later).
%ptp - 3xn matrix containing phi/theta/psi in degrees
%Note: if this function is interrupted, then there is a udp port which does
%not get closed correctly. These will accumulate over time, consuming small
%amounts of memory and CPU. Restarting MATLAB will resolve any accumulated
%open ports. The performance impact should be minimal.
uByte = udpport("byte");

%VEHX preamble
preamble = ['VEHX', 0x00];

%Make the time vector start at 0
t = t-t(1);
t_max = t(end);

datetime_start = datetime();
n = -1;

while (datetime()-datetime_start)*86400 < t_max
    %Find which data entry we currently need to send
    n_prvs = n;
    n = find(t-(datetime()-datetime_start)*86400>0,1)-1;
    if n == n_prvs
        %If we've handled this data already, don't clog the udp ports
        continue;
    end
    %assemble and send the data
    data = [typecast(int32(0),'uint8'), ...
        typecast(double(lla(1,n)),'uint8') typecast(double(lla(2,n)),'uint8') typecast(double(lla(3,n)),'uint8') ...
        typecast( single(ptp(3,n)),'uint8') typecast( single(ptp(2,n)),'uint8') typecast( single(ptp(1,n)),'uint8')];
    write(uByte,[preamble data],'127.0.0.1',49000);
end

%At the end, return the aircraft to the start position after a second
pause(1);
n = 1;
data = [typecast(int32(0),'uint8'), ...
    typecast(double(lla(1,n)),'uint8') typecast(double(lla(2,n)),'uint8') typecast(double(lla(3,n)),'uint8') ...
    typecast( single(ptp(3,n)),'uint8') typecast( single(ptp(2,n)),'uint8') typecast( single(ptp(1,n)),'uint8')];
 write(uByte,[preamble data],'127.0.0.1',49000);

clear uByte

%VEHX format:
%VEHX 0x00 lat lon ele psi the phi