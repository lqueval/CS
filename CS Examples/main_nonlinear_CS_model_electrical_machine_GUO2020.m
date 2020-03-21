%---------------------------------------------------
%  NAME:      main_nonlinear_CS_model_electrical_machine_GUO2020.m
%  WHAT:      Nonlinear current sheet of an electrical machine (https://doi.org/10.1109/TMAG.2019.2950614)
%  REQUIRED:  CSmodel toolbox 20200321
%  AUTHOR:    20200321, J. Guo, L. Quéval (loic.queval@gmail.com), F. Trillaud, B. Roucaries
%  COPYRIGHT: 2020, Loïc Quéval, BSD License (http://opensource.org/licenses/BSD-3-Clause)
%----------------------------------------------------

clear all, close all, clc

addpath('..\CS Core') %add CSmodel toolbox functions to the path

%% Constants

mu0 = 4*pi*1e-7;

%% Variables

isnonlin = 1; %0 for linear case, 1 for nonlinear case
h_min = 1; %min harmonic number to consider (must be odd)
h_max = 13; %max harmonic number to consider

%% Define machine

% Machine geometry
P = 6;             %number pole pair []
r1 = 1.320;        %rotor inner radius[m]
r2 = 1.470;        %rotor outer radius [m]
r3 = 1.546;        %radius of field coil[m]
rTe = 1.619;       %radius for estimation of Te [m]
r4 = 1.683;        %radius of amature coil [m]
r5 = 1.750;        %stator inner radius[m]
r6 = 2.000;        %Stator outer radius[m]

% Rotor winding
Nf = 100;           %field coil turns []
hf = 0.057;         %field coil height [m]
theta_1f = 0.163;   %field coil width angle [elec. rad]
theta_2f = 2.703;   %field coil aperture angle [elec. rad]
wf = 0.042;         %field coil width [m]
alfa_m_deg_ini = -15;	%initial rotor angle [mec. deg]
alfa_m_rad_ini = alfa_m_deg_ini*pi/180;   %rotor angle [mec. rad]
i_f = 5.03e3;      %field winding current [A]

% Stator winding
Na= 120;            %amature coil turns []
ha = 0.041;         %armature coil height[m]
theta_1a = 0.692;   %armature coil width angle [elec. rad]
theta_2a = 0.664;   %armature coil aperture angle [elec. rad]
wa = 0.194;         %armature coil width [m]
i_a = -1.53e3;     %phase a current [A]
i_b = 2.465e3;     %phase b current [A]
i_c = -0.935e3;    %phase c current [A]

% iron
fct_mur_B([]); %Define iron (and plot)

% build machine geom
r_l = [r1,r2,r3,r4,r5,r6]; %vector of current sheet radius
is_currentsheet = [0,0,1,1,0,0]; %0 for nothing, 1 for current sheet
mur_l = [1,1200,1,1,1,1200,1]; %vector of domain (initial) relative permeability (air, rotor iron, airgap, airgap, airgap, stator iron, air)
is_nonlinear = [0,1,0,0,0,1,0]*isnonlin; %0 for linear, 1 for nonlinear domain

% plot machine geometry
fig1 = figure();
box on, grid on, hold on
fct_plot_geom_generic(r_l,is_currentsheet)
xlabel('x [m]'), ylabel('y [m]')

%% Solution

tic

% Initialize current sheets
K_lh_s = zeros(length(r_l),h_max);
K_lh_c = zeros(length(r_l),h_max);

%get current sheet rotor winding
[K_lh_c(3,:),K_lh_s(3,:)] = fct_get_K3(h_min,h_max,Nf,i_f,theta_1f,theta_2f,alfa_m_rad_ini,P,wf);

%get current sheet stator winding
[K_lh_c(4,:),K_lh_s(4,:)] = fct_get_K4(Na,i_a,i_b,i_c,h_min,h_max,theta_1a,theta_2a,wa);

%get (a,b,c,d) in the nonlinear case
[a_lh,b_lh,c_lh,d_lh,mur_l] = fct_get_XX_nonlin(h_min,h_max,P,r_l,mur_l,mu0,K_lh_s,K_lh_c,is_nonlinear);

toc

%% Postprocessing and export

% Plot 2D map of NORMB
[THETA,RHO] = ndgrid((0:0.5:60)*pi/180,0:r6/100:2*r6); %define field points in sector
[RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA); % Get (BR,BTH) at field points (RHO,THETA)

figure(); box on, grid on, hold on
colormap jet;
[X,Y] = pol2cart(THETA,RHO);
nb_contour = linspace(0,5,100);  %100 contours between 0 and 5T
    [~,h] = contourf(X,Y,NORMB(:,:,1),nb_contour); %plot NORMB
set(h,'edgecolor','none');
    fct_plot_geom_generic(r_l,is_currentsheet) %plot geometry
axis tight; xlim([0,1.2*r_l(end)]); ylim([0,1.2*r_l(end)]);
xlabel('x [m]'), ylabel('y [m]')
cbh = colorbar; set(get(cbh,'Title'),'String','|B| [T]')
caxis([0,5])

% Plot 1D air-gap magnetic field BR and BTH
[THETA,RHO] = ndgrid((-180:0.5:180)*pi/180,rTe); %define field points at rRTe
[RHO,THETA,BR,BTH,NORMB] = fct_get_B(a_lh,b_lh,c_lh,d_lh,h_min,h_max,P,r_l,RHO,THETA); % Get (BR,BTH) at field points (RHO,THETA)

figure();
subplot(2,1,1), box on, grid on, hold on % Plot Br
    plot(THETA*180/pi,BR,'-k')
ylabel('B_r [T]')
xlim([0,60]), ylim([-3,2.5])
subplot(2,1,2), box on, grid on, hold on % Plot Bth
    plot(THETA*180/pi,BTH,'-k');
ylabel('B_{\theta} [T]')
xlim([0,60]), ylim([-3,2.5])
xlabel('\theta [mec. deg]')