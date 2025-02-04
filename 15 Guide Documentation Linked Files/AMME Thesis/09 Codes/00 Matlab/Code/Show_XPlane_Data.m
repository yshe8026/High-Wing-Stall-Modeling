clear all; close all;

%load data.txt using the matlab import tool, default settings are fine.
%Rename the resulting table to Data, so that the references to this
%variable below work.
load data.mat

figure()
tiledlayout(4,2)
ax = [];

%% airspeed
ax = [ax nexttile()];
plot(Data.real_time, Data.Vtrue_ktas);
ylabel('V_{true} (kts)');
xlabel('Time (s)');

%% aerodynamic angles
ax = [ax nexttile()];
plot(Data.real_time, Data.alpha__deg); hold on;
plot(Data.real_time, Data.beta__deg); hold on;
legend('\alpha','\beta');
ylabel('Aerodynamic Angles (deg)');
xlabel('Time (s)');

%% angular rates
ax = [ax nexttile()];
plot(Data.real_time, Data.Prads*180/pi); hold on;
plot(Data.real_time, Data.Qrads*180/pi); hold on;
plot(Data.real_time, Data.Rrads*180/pi); hold on;
legend('p','q','r');
ylabel('Angular Rates (deg/s)');
xlabel('Time (s)');

%% Attitude
ax = [ax nexttile()];
plot(Data.real_time, Data.roll__deg); hold on;
plot(Data.real_time, Data.pitch__deg); hold on;
plot(Data.real_time, Data.hding_true); hold on;
legend('\phi','\theta','\psi');
ylabel('Euler Angles (deg)');
xlabel('Time (s)');

%% altitude
ax = [ax nexttile()];
plot(Data.real_time, Data.altftmsl);
ylabel('Altitude (ft)');
xlabel('Time (s)');

%% controls
ax = [ax nexttile()];
plot(Data.real_time, Data.elev_surf); hold on;
plot(Data.real_time, Data.ailrn_surf); hold on;
plot(Data.real_time, Data.ruddr_surf); hold on;
legend('elevator','aileron','rudder');
ylabel('Control Deflections (normalised)');
xlabel('Time (s)');

%% moments
ax = [ax nexttile()];
plot(Data.real_time, Data.L_ftlb); hold on;
plot(Data.real_time, Data.M_ftlb); hold on;
plot(Data.real_time, Data.N_ftlb); hold on;
%The following are the *aerodynamic* moments (above is the *total* moments,
%which may include moments from the propulsion system and possibly
%gyroscopic moments created by moving the angular momentum vector of the
%rotating parts of the propulsion system)
% plot(Data.real_time, Data.Llbft); hold on;
% plot(Data.real_time, Data.Mlbft); hold on;
% plot(Data.real_time, Data.Nlbft); hold on;
legend('L','M','N');
ylabel('Moment (ftlbs)');
xlabel('Time (s)');

%% aerodynamic forces
ax = [ax nexttile()];
plot(Data.real_time, Data.drag___lb); hold on;
plot(Data.real_time, Data.side___lb); hold on;
plot(Data.real_time, Data.lift___lb); hold on;
legend('D','Y','L');
ylabel('Aerodynamic Forces (lb)');
xlabel('Time (s)');



linkaxes(ax,'x');