clear all; close all;

%Set up some variables for packet ID

VEHX_ID_str = [uint8('VEHX') 0x00];
DREF_ID_str = [uint8('DREF') 0x00];
DATA_ID_str = [uint8('DATA') 0x00];

%% DREF strings
DREF_override_ias = uint8('sim/operation/override/override_ias');
DREF_override_ias(end+2:500) = 0x20; %packet must be padded with ascii space character (not null)

DREF_calibrated_airspeed_kts_pilot = uint8('sim/cockpit2/gauges/indicators/calibrated_airspeed_kts_pilot');
DREF_calibrated_airspeed_kts_pilot(end+2:500) = 0x20;

DREF_airspeed_kts_pilot = uint8('sim/cockpit2/gauges/indicators/airspeed_kts_pilot');
DREF_airspeed_kts_pilot(end+2:500) = 0x20;