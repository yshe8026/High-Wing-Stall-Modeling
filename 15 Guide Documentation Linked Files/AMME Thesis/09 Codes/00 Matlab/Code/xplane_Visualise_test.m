%% Brendon Example
% clear all; close all;
% %test the x plane visualise function
% n = 1000;
% t = linspace(0,10,n);
% 
% lat = linspace(-33,-34,n);
% lon = linspace(151,152,n);
% alt = linspace(1000,20000,n);
% 
% phi = linspace(-180,180,n);
% theta = linspace(0,1,n);
% psi = linspace(-0,1,n);
% 
% lla = [lat;lon;alt];
% ptp = [phi;theta;psi];
% 
% xplane_Visualise(t,lla,ptp);


%% Richard To Outter Space
% clear all; close all;
% %test the x plane visualise function
% n = 10000;
% t = linspace(0,100,n);
% 
% lat = linspace(-33.856159,-33.856159,n);
% lon = linspace(151.215256,151.215256,n);
% alt = linspace(1000,10000000,n);
% alt = (alt./100).^2; 
% alt = alt + 1000;
% 
% phi = linspace(0,0,n);
% theta = linspace(-90,-90,n);
% psi = linspace(0,0,n);
% 
% lla = [lat;lon;alt];
% ptp = [phi;theta;psi];
% 
% xplane_Visualise(t,lla,ptp);


%% Richard To Outter Space
clear all; close all;
% Load only specific variables from the .mat file
load('air_race_t.mat', 't');
load('air_race_X_mat.mat', 'X_mat');
% Create Euler angle matrix
quat = X_mat(7:10,:);
euler =  q2e(quat);
euler = (180/pi).*euler;
% Create a 3 x nSteps matrix to store [lat (rad) ; lon (rad) ; alt (m)] each timestep


% Accelerate time
% t = t./1;
%test the x plane visualise function
n = length(t);
% t = linspace(0,5,n);

lat0 = -33.856159;
lon0 = 151.215256;
local_cosine = cos(lat0*(pi/180));
earth_circumference = 12756000;
local_circumference = earth_circumference * local_cosine;
lat_m2deg = 360/earth_circumference;
lon_m2deg = 360/local_circumference;

lat = X_mat(11,:) * lat_m2deg;
lat = lat + lat0;
lon = X_mat(12,:) * lon_m2deg;
lon = lon + lon0;
alt = -X_mat(13,:);
alt = alt + 200;

phi = euler(1,:);
theta = euler(2,:);
psi = euler(3,:);

lla = [lat;lon;alt];
ptp = [phi;theta;psi];

xplane_Visualise(t,lla,ptp);