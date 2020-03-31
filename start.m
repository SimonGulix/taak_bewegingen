
clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 1;           % draw figures of kinematic analysis if 1
fig_kin_check = 1; 
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1
fig_dyn_check = 1;
fig_dyn_check_shaking=1;

% kinematic parameters (link lengths)
r11 = 2*22*(10)^(-2);
r12 = 2*sqrt(185)*(10)^(-2);
r13 = 2*sqrt(185)*(10)^(-2);
r2 = 2*4*(10)^(-2);
r3 = 2*5*(10)^(-2);
r4 = 2*13*(10)^(-2);
r5 = 2*10*(10)^(-2);
r6 = 2*10*(10)^(-2);
r7 = 2*13*(10)^(-2);
r8 = 2*5*(10)^(-2);
r9a =2* 14*(10)^(-2);
r9b =2* 6*(10)^(-2);
r10a =2* 14*(10)^(-2);
r10b =2* 6*(10)^(-2);
CA_x=2*(-11)*(10)^(-2);%vector CA
CA_y=2*8*(10)^(-2);
CB_x=2*11*(10)^(-2);%vector CB
CB_y=2*8*(10)^(-2);

g = 9.81;

phi12 = convert_radial(323.9726266);%calculated via Pythagoras
phi13 = convert_radial(216.0273734);
phi11 = 0;

density=8000; %kg/m3
Asurf=(pi/4)*0.03^2; %3cm diameter
m11 = r11*density*Asurf; 
m12 = r12*density*Asurf;
m13 = r13*density*Asurf;
m2 = r2*density*Asurf;
m3 = r3*density*Asurf;
m4 = r4*density*Asurf;
m5 = r5*density*Asurf;
m6 = r6*density*Asurf;
m7 = r7*density*Asurf;
m8 = r8*density*Asurf;
m9 = (r9a+r9b)*density*Asurf;
m10 = (r10a+r10b)*density*Asurf;


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

% % % position analysis PHI2_INIT= 225
phi3_init = convert_radial(180+45);    % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi4_init = convert_radial(180+75);  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi5_init = convert_radial(30);
phi6_init = convert_radial(160);
phi7_init = convert_radial(280);
phi8_init = convert_radial(180+75);
phi9_init = convert_radial(35);
phi10_init = convert_radial(150);


% % % ALTERNATIVE START PHI2=0
% phi3_init = convert_radial(290);    
% phi4_init = convert_radial(280);  
% phi5_init = convert_radial(55);
% phi6_init = convert_radial(125);
% phi7_init = convert_radial(285);
% phi8_init = convert_radial(350);
% phi9_init = convert_radial(25);
% phi10_init = convert_radial(170);



t_begin = 0;                   % start time of simulation
t_end = 2*pi;                    % end time of simulation
Ts = 0.1;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector
omega = 1;
A = 1;
% phi2= omega*t;
% dphi2=omega* ones(size(t));
% ddphi2 = zeros(size(t));

phi2= pi*(1-cos((t/2))); %%uit theorie van nokken gehaald( Harmonic H5-lecture 7, slide 45)
dphi2=(pi/2)*sin(t/2);
ddphi2 =(pi/4)*cos(t/2);


% thau = t/t_end;
% 
% phi2 = 2*pi*(3*thau.^2 - 2*thau.^3);
% dphi2 = 2*pi*(6*thau/t_end - 6*thau.^2 /t_end);
% ddphi2 = 2*pi*(6/t_end^2 - 12*thau/t_end^2);

% phi2 = 2*pi*(thau - 1/(2*pi) * sin(2*pi*thau));
% dphi2 = 2*pi*(1/t_end - (1/t_end * cos(2*pi*thau)));
% ddphi2 =2*pi*(1/t_end^2 *2*pi* sin(2*pi*thau));

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
    alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, alpha8,alpha9,alpha10,...
    CF_vec,AE_vec,EI_vec,IG_vec,JH_vec,DJ_vec,BD_vec,GD_vec,HE_vec] = ...
    dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10,r2,r3,r4,r5,r6,r7,r8,r9a,r9b,r10a,r10b, ...
  m2,m3,m4,m5,m6,m7,m8,m9,m10,J2,J3,J4,J5,J6,J7,J8,J9,J10,g,t,fig_dyn_4bar);


dyn_check(vel_2x,vel_2y,vel_3x,vel_3y,vel_4x,vel_4y,vel_5x,vel_5y,vel_6x,vel_6y,vel_7x,vel_7y,vel_8x,vel_8y,vel_9x,vel_9y,vel_10x,vel_10y,...
                   acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,acc_9x,acc_9y,acc_10x,acc_10y,...
                   M_C,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,m2,m3,m4,m5,m6,m7,m8,m9,m10,...
                   J2,J3,J4,J5,J6,J7,J8,J9,J10,t,fig_dyn_check,g);
               
dyn_check_shaking(vel_2x,vel_2y,vel_3x,vel_3y,vel_4x,vel_4y,vel_5x,vel_5y,vel_6x,vel_6y,vel_7x,vel_7y,vel_8x,vel_8y,vel_9x,vel_9y,vel_10x,vel_10y,...
                   acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,acc_9x,acc_9y,acc_10x,acc_10y,...
                   M_C,F_A_x,F_A_y,F_C_x,F_C_y,F_B_x,F_B_y,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,m2,m3,m4,m5,m6,m7,m8,m9,m10,...
                   J2,J3,J4,J5,J6,J7,J8,J9,J10,CF_vec,AE_vec,EI_vec,IG_vec,JH_vec,DJ_vec,BD_vec,GD_vec,HE_vec,CA_x,CA_y,CB_x,CB_y,t,fig_dyn_check_shaking,g);
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
    ang = (360/(2*pi)) * angle;
end

