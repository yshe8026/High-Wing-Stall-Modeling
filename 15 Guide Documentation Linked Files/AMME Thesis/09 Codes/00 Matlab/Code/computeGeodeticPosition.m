% Function to compute geodetic position from aircrafts state

% Reference: 
% E. Johnson, F. Lewis, and B. Stevens: "Aircraft control and simulation : 
%      dynamics, controls design, and autonomous systems", 3rd Ed, Wiley-Blackwell, 2015 
% Navigation Equations: pp. 111
% Geodetic Equations: pp. 23-29

% Inputs
%   - u, x-axis velocity [m/s] (body axis)
%   - v, y-axis velocity [m/s] (body axis)
%   - w, z-axis velocuty [m/s] (body axis)

%   - phi, roll angle [rad]
%   - theta, pitch angle [rad]
%   - psi, yaw angle [rad]

%   - dt, temporal step [s]

%   - lla_o, 3x1 vector containing [lat; lon; alt] from previous timestep

function [lla] = computeGeodeticPosition(u,v,w,phi,theta,psi,dt,lla_o)

    %% Calculate Euler angle and position derivatives for state vector  
    % i.e. using flat-earth, body axis, 6-DOF navigation equations:
    %   - dot(d_n) = v_n    = Uc𝜃c𝜓 + V (−c𝜙s𝜓 + s𝜙s𝜃c𝜓) + W(s𝜙s𝜓 + c𝜙s𝜃c𝜓)
    V_N   = u*cos(theta)*cos(psi) + ...
            v*(-cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi) ) + ...
            w*( sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi) );
    
    %   - dot(d_e) = v_e    = Uc𝜃s𝜓 + V(c𝜙c𝜓 + s𝜙s𝜃s𝜓) + W(−s𝜙c𝜓 + c𝜙s𝜃s𝜓)
    V_E   = u*cos(theta)*sin(psi) + ...
            v*( cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi) ) + ...
            w*(-sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi) );

    %   - dot(d_v) = h_dot  = Us𝜃 − Vs𝜙c𝜃 − Wc𝜙c𝜃
    h_dot = u*sin(theta) - v*sin(phi)*cos(theta) - w*cos(phi)*cos(theta);       %!!!! DOUBLE CHECK CORRECT DIRECTION !!!!!!
    

    %% Calculate Geodetic Lat, Lon, Alt
    % Constants
    a = 6378137.00;             % Semi-major Axis (m)
    b = 6356752.00;             % Semi-minor Axis (m)
    e = sqrt(a^2 - b^2)/a;      % First Eccentricity
    
    % Calculated Variables
    M = (a*(1-e^2))/((1-e^2*(sin(lla_o(1))^2))^(3/2)); % Meridian Radius of Curvature (m)
    N = a/sqrt(1 - e^2 * (sin(lla_o(1))^2));           % Prime Vertical Radius of Curvature (m)

    % Derivative of latitude -> dot(phi) = V_N/(M + h)
    lat_dot = V_N/(M+lla_o(3));  

    % Derivative of longitude -> dot(lambda) = V_E/(N + h)cos(phi)
    lon_dot = V_E/((N + lla_o(3))*cos(lla_o(1)));

    %% Euler Integration
    lat = lla_o(1) + lat_dot*dt;
    lon = lla_o(2) + lon_dot*dt;
    alt = lla_o(3) + h_dot*dt;
    
    % bind longitude @ north/south pole
    lon = mod(lon+pi, 2*pi) - pi;

    
    %% Output
    lla = [lat; lon; alt];

end