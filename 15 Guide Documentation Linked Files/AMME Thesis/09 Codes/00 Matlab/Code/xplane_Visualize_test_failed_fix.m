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
quat_mat = X_mat(7:10,:);
euler_mat =  q2e(quat_mat);
euler_mat = (180/pi).*euler_mat;

% Create a 3 x nSteps matrix to store [lat (rad) ; lon (rad) ; alt (m)] each timestep

% Accelerate time
% t = t./1;
%test the x plane visualise function
n = length(t);
% t = linspace(0,5,n);

% % Geodesic coordinate of sydney harbour bridge
% lat0 = -33.856159;
% lon0 = 151.215256;
% alt0 = -X_mat(13,1) + 400;

% Geodesic coordinate of Burj Khalifacc
lat0 = deg2rad(25.1972);
lon0 = deg2rad(55.2744);
alt0 = -X_mat(13,1) + 400;

% initial 3x1 vector containing [lat0; lon0; alt0]
lla = [lat0; lon0; alt0];
% time step size
dt = 0.01;

for i = 1:n

lla_o = lla;
u = X_mat(1,i);
v = X_mat(2,i);
w = X_mat(3,i);
phi = euler_mat(1,i);
theta = euler_mat(2,i);
psi = euler_mat(3,i);
[lla] = computeGeodeticPosition(u,v,w,phi,theta,psi,dt,lla_o);
lla_mat(:,i) = lla;

end

% lat = X_mat(11,:) * lat_m2deg;
% lat = lat + lat0;
% lon = X_mat(12,:) * lon_m2deg;
% lon = lon + lon0;
% alt = -X_mat(13,:);
% alt = alt + 200;

phi_series = euler_mat(1,:);
theta_series = euler_mat(2,:);
psi_series = euler_mat(3,:);

lla_mat(1,:) = rad2deg(lla_mat(1,:));
lla_mat(2,:) = rad2deg(lla_mat(2,:));
ptp_mat = [phi_series;theta_series;psi_series];

xplane_Visualise(t,lla_mat,ptp_mat);