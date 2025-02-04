clear all; close all;

%Set up some variables for writing datarefs

VEHX_ID_str = [uint8('VEHX') 0x00];
DREF_ID_str = [uint8('DREF') 0x00];
DATA_ID_str = [uint8('DATA') 0x00];



uByte = udpport("byte");

%% write VEHX
lla = [-30, 151, 1000]; %lat long altitude, deg deg metres respectively. lat/lon are apparently in WGS84
ptp = [-10, 10, 45]; %phi theta psi, deg deg deg respectively

preamble = ['VEHX', 0x00];
data = [typecast(int32(0),'uint8'), ...
        typecast(double(lla(1)),'uint8') typecast(double(lla(2)),'uint8') typecast(double(lla(3)),'uint8') ...
        typecast( single(ptp(3)),'uint8') typecast( single(ptp(2)),'uint8') typecast( single(ptp(1)),'uint8')]; %note the order here is psi, theta, phi but in AERO3560 we normally order these as phi, theta, psi
write(uByte,[preamble data],'127.0.0.1',49000);

%% write DREF

DREF_override_ias = uint8('sim/operation/override/override_ias');
DREF_override_ias(end+2:500) = 0x20; %packet must be padded with ascii space character (not null)
write(uByte,[DREF_ID_str typecast(single(1),'uint8') DREF_override_ias],'127.0.0.1',49000);

DREF_calibrated_airspeed_kts_pilot = uint8('sim/cockpit2/gauges/indicators/calibrated_airspeed_kts_pilot');
DREF_calibrated_airspeed_kts_pilot(end+2:500) = 0x20;
write(uByte,[DREF_ID_str typecast(single(100),'uint8') DREF_calibrated_airspeed_kts_pilot],'127.0.0.1',49000);

DREF_airspeed_kts_pilot = uint8('sim/cockpit2/gauges/indicators/airspeed_kts_pilot');
DREF_airspeed_kts_pilot(end+2:500) = 0x20;
write(uByte,[DREF_ID_str typecast(single(100),'uint8') DREF_airspeed_kts_pilot],'127.0.0.1',49000);

%% write DATA
%Write to the airspeeds in the general data output packet #3
packet_ID = 3;
IAS = 100;
write(uByte,[DATA_ID_str typecast(int32(packet_ID),'uint8') ...
    typecast(single(IAS),'uint8') typecast(single(IAS),'uint8') typecast(single(IAS),'uint8') typecast(single(IAS),'uint8') ...
    typecast(single(IAS),'uint8') typecast(single(IAS),'uint8') typecast(single(IAS),'uint8') typecast(single(IAS),'uint8')],'127.0.0.1',49000);

clear uByte