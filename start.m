%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 0;           % draw figures of kinematic analysis if 1
fig_kin_check = 1;        
%fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r11 = 22;
r12 = sqrt(185);
r13 = sqrt(185);
r2 = 4;
r3 = 5;
r4 = 13;
r5 = 10;
r6 = 10;
r7 = 13;
r8 = 5;
r9a = 14;
r9b = 6;
r10a = 14;
r10b = 6;


phi12 = convert_radial(323.9726266);
phi13 = convert_radial(216.0273734);
phi11 = 0;
% % dynamic parameters, defined in a local frame on each of the bars.
% X2 = r2/2;               % X coordinates of cog (centre of gravity)
% X3 = r3/2;
% X4 = r4/2;
% 
% Y2 = 0;                  % Y coordinates of cog
% Y3 = 0.0102362;
% Y4 = 0;
% 
% m2 = r2*1.76;
% m3 = r3*1.76;
% m4 = r4*0.54;
% 
% J2 = m2*r2^2/12;
% J3 = m3*r3^2/12;
% J4 = m4*r4^2/12;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% position analysis
phi3_init = convert_radial(180+45);    % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi4_init = convert_radial(180+75);  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi5_init = convert_radial(30);
phi6_init = convert_radial(160);
phi7_init = convert_radial(280);
phi8_init = convert_radial(180+75);
phi9_init = convert_radial(35);
phi10_init = convert_radial(150);

t_begin = 0;                   % start time of simulation
t_end = 2*pi;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector
%disp("t: " + t);
% initialization of driver
omega = 2;
A = 1;
phi2= omega*t;
dphi2=omega* ones(size(t));
ddphi2 = zeros(size(t));
% calculation of the kinematics (see kin_4bar.m)

[phi3,phi4,phi5, phi6, phi7, phi8,phi9,phi10,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,ddphi3,ddphi4, ddphi5, ddphi6, ddphi7, ddphi8, ddphi9, ddphi10] = kinematics_4bar(r11, r12, r13,r2,r3,r4,r5,r6,r7,r8,r9a, r9b,r10a, r10b, phi11, phi12, phi13,phi2,dphi2,ddphi2, phi3_init, phi4_init, phi5_init, phi6_init, phi7_init, phi8_init, phi9_init, phi10_init,t,fig_kin_4bar);

 disp("Done Kinematics");
 
 VarNames = {'phi2', 'phi3', 'phi4', 'phi5', 'phi6', 'phi7', 'phi8', 'phi9', 'phi10'};
T = table(convert_to_degree(phi2), phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10, 'VariableNames',VarNames);
 %disp(T);
 
 % verify kinematics (see kin_check.m)
[diffphi3,diffphi4,diffphi5,diffphi6,diffphi7,diffphi8,diffphi9,diffphi10,ddiffphi3,ddiffphi4,ddiffphi5,ddiffphi6,ddiffphi7,ddiffphi8,ddiffphi9,ddiffphi10 ] = ...
 check_kinematics(phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,Ts, fig_kin_check, t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
% [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
%   m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
%  load fourbar_movie Movie
%  movie(Movie)

function ang = convert_radial(angle) 
    ang = ((2*pi)/360) * angle;
end

function ang = convert_to_degree(angle) 
    ang = (360/(2*pi))*angle;
end

