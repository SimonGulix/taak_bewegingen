
clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 0;           % draw figures of kinematic analysis if 1
fig_kin_check = 0;        
fig_dyn_4bar = 0;        % draw figures of dynamic analysis if 1
fig_dyn_check = 0;

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
g = 9.81;


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
density=8; %kg/dm^3  && each beam Cross-section of 1dm x 1dm
m11 = r11*density;
m12 = r12*density;
m13 = r13*density;
m2 = r2*density;
m3 = r3*density;
m4 = r4*density;
m5 = r5*density;
m6 = r6*density;
m7 = r7*density;
m8 = r8*density;
m9 = (r9a+r9b)*density;
m10 = (r10a+r10b)*density;


J2 = m2*r2^2/12;
J3 = m3*r3^2/12;
J4 = m4*r4^2/12;
J5 = m5*r5^2/12;
J6 = m6*r6^2/12;
J7 = m7*r7^2/12;
J8 = m8*r8^2/12;
J9 = m9*(r9a+r9b)^2/12;
J10 = m10*(r10a+r10b)^2/12;


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
Ts = 0.1;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector
%disp("t: " + t);
% initialization of driver

thau = t/t_end;
omega = 1;
A = 1;
phi2= omega*t;
dphi2=omega* ones(size(t));
ddphi2 = zeros(size(t));


% phi2 = 10*(3*thau.^2 - 2*thau.^3);
% dphi2 = 10*(6*thau/t_end - 6*thau.^2 /t_end);
% ddphi2 = 10*(6/t_end^2 - 12*thau/t_end^2);

% phi2 = 10*(thau - 1/(2*pi) * sin(2*pi*thau));
% 
% dphi2 = 10*(1/t_end - (1/t_end * cos(2*pi*thau)));
% ddphi2 =10*(1/t_end^2 *2*pi* sin(2*pi*thau));
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
[F_A_x,F_A_y,F_I_x,F_I_y,F_G_x,F_G_y,F_C_x,F_C_y,F_B_x,F_B_y,F_J_x,F_J_y,F_H_x,F_H_y,F_D9_x,F_D9_y,F_D7_x,F_D7_y,F_D8_x,F_D8_y,F_F10_x,F_F10_y,F_F9_x,F_F9_y,F_F2_x,F_F2_y,F_E4_x,F_E4_y,F_E3_x,F_E3_y,F_E10_x,F_E10_y,M_C,...
    vel_2x, vel_2y, vel_3x, vel_3y, vel_4x, vel_4y, vel_5x, vel_5y, vel_6x, vel_6y, vel_7x, vel_7y, vel_8x, vel_8y, vel_9x, vel_9y, vel_10x, vel_10y,...
    acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,acc_9x,acc_9y,acc_10x,acc_10y,...
    omega2, omega3, omega4, omega5, omega6, omega7, omega8, omega9, omega10,...
    alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8,alpha9,alpha10] = ...
    dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,r2,r3,r4,r5,r6,r7,r8,r9a,r9b,r10a,r10b, ...
  m2,m3,m4,m5,m6,m7,m8,m9,m10,J2,J3,J4,J5,J6,J7,J8,J9,J10,t,fig_dyn_4bar);


dyn_check(vel_2x,vel_2y,vel_3x,vel_3y,vel_4x,vel_4y,vel_5x,vel_5y,vel_6x,vel_6y,vel_7x,vel_7y,vel_8x,vel_8y,vel_9x,vel_9y,vel_10x,vel_10y,...
                   acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,acc_9x,acc_9y,acc_10x,acc_10y,...
                   M_C,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,m2,m3,m4,m5,m6,m7,m8,m9,m10,...
                   J2,J3,J4,J5,J6,J7,J8,J9,J10,t,fig_dyn_check,g);
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

