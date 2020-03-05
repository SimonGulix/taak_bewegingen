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


function [F_P_x,F_Q_x,F_R_x,F_S_x,F_P_y,F_Q_y,F_R_y,F_S_y,M_P] = ...
dynamics_4bar(phi2,phi3,phi4,dphi2,dphi3,dphi4,ddphi2,ddphi3,ddphi4,r2,r3,r4, ...
    m2,m3,m4,X2,X3,X4,Y2,Y3,Y4,J2,J3,J4,t,fig_dyn_4bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
cog2_C_x = -r2/2 * cos(phi2);
cog2_C_y = -r2/2* sin(phi2);
cog2_F_x = r2/2*cos(phi2);
cog2_F_y = r2/2*sin(phi2);

cog3_A_x = -r3/2*cos(phi3);
cog3_A_y = -r3/2*cos(phi3);
cog3_E_x = r3/2*cos(phi3);
cog3_E_y = r3/2*cos(phi3);

cog4_E_x = -r4/2*cos(phi4);
cog4_E_y = -r4/2*sin(phi4);
cog4_I_x = r4/2*cos(phi4);
cog4_I_y = r4/2*sin(phi4);

cog5_I_x = -r5/2*cos(phi5);
cog5_I_y = -r5/2*sin(phi5);
cog5_G_x = r5/2*cos(phi5);
cog5_G_y = r5/2*sin(phi5);

cog6_J_x = -r6/2*cos(phi6);
cog6_J_y = -r6/2*sin(phi6);
cog6_H_x = r6/2*cos(phi6);
cog6_H_y = r6/2*sin(phi6);

cog7_D_x = -r7/2*cos(phi7);
cog7_D_y = -r7/2*sin(phi7);
cog7_J_x  = r7/2*cos(phi7);
cog7_J_y = r7/2*sin(phi7);

cog8_B_x = -r8/2*cos(phi8);
cog8_B_y = -r8/2*sin(phi8);
cog8_D_x = r8/2*cos(phi8);
cog8_D_y = r8/2*sin(phi8);

cog9_G_x = -(r9a+r9b)/2*cos(phi9);
cog9_G_y = -(r9a+r9b)/2*sin(phi9);
cog9_D_x = (r9a+r9b)/2*cos(phi9);
cog9_D_y = (r9a+r9b)/2*sin(phi9);

cog10_H_x = -(r10a+r10b)/2*cos(phi10);
cog10_H_y = -(r10a+r10b)/2*sin(phi10);
cog10_E_x = (r10a+r10b)/2*cos(phi10);
cog10_E_y = (r10a+r10b)/2*sin(phi10);

% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
omega8 = [zeros(size(phi2)) zeros(size(phi2)) dphi8];
omega9 = [zeros(size(phi2)) zeros(size(phi2)) dphi9];
omega10 = [zeros(size(phi2)) zeros(size(phi2)) dphi10];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
alpha8 = [zeros(size(phi2)) zeros(size(phi2)) ddphi8];
alpha9 = [zeros(size(phi2)) zeros(size(phi2)) ddphi9];
alpha10 = [zeros(size(phi2)) zeros(size(phi2)) ddphi10];

% 3D model vectors
C_cog2_vec = [-cog2_C_x    -cog2_P_y    zeros(size(phi2))];
A_cog3_vec = [-cog3_A_x    -cog3_A_y    zeros(size(phi2))];
E_cog4_vec = [-cog4_E_x    -cog4_E_y    zeros(size(phi2))];
I_cog5_vec = [-cog5_I_x    -cog5_I_y    zeros(size(phi2))];
J_cog6_vec = [-cog6_J_x    -cog6_J_y    zeros(size(phi2))];
D_cog7_vec = [-cog7_D_x    -cog7_D_y    zeros(size(phi2))];
B_cog8_vec = [-cog8_B_x    -cog8_B_y    zeros(size(phi2))];
G_cog9_vec = [-cog9_G_x    -cog9_G_y    zeros(size(phi2))];
H_cog10_vec = [-cog10_H_x    -cog10_H_y    zeros(size(phi2))];

%PQ_vec = [r2*cos(phi2) r2*sin(phi2) zeros(size(phi2))];
CF_vec = [r2*cos(phi2) r2*sin(phi2) zeros(size(phi2))];
AE_vec = [r3*cos(phi3) r3*sin(phi3) zeros(size(phi2))];
EI_vec = [r4*cos(phi4) r4*sin(phi4) zeros(size(phi2))];
IG_vec = [r5*cos(phi5) r5*sin(phi5) zeros(size(phi2))];
JH_vec = [r6*cos(phi6) r6*sin(phi6) zeros(size(phi2))];
DJ_vec = [r7*cos(phi7) r7*sin(phi7) zeros(size(phi2))];
BD_vec = [r8*cos(phi8) r8*sin(phi8) zeros(size(phi2))];
GD_vec = [(r9a+r9b)*cos(phi9) (r9a+r9b)*sin(phi9) zeros(size(phi2))];
HE_vec = [(r10a+r10b)*cos(phi10) (r10a+r10b)*sin(phi10) zeros(size(phi2))];

% acceleration vectors
acc_2 =       cross(omega2,cross(omega2,C_cog2_vec))+cross(alpha2,C_cog2_vec);
acc_3 =       cross(omega3,cross(omega3,A_cog3_vec))+cross(alpha3,A_cog3_vec);
acc_E =       cross(omega3,cross(omega3,AE_vec))+cross(alpha3,AE_vec);
acc_4 = acc_E+cross(omega4,cross(omega4,E_cog4_vec))+cross(alpha4,E_cog4_vec);
acc_I = acc_E+cross(omega4,cross(omega4,EI_vec))+cross(alpha4,EI_vec    );
acc_5 = acc_I+cross(omega5,cross(omega5,I_cog5_vec))+cross(alpha5,I_cog5_vec);

acc_8 =      cross(omega8,cross(omega8,B_cog8_vec))+cross(alpha8,B_cog8_vec);
acc_D =      cross(omega8,cross(omega8,BD_vec))+cross(alpha8,BD_vec);
acc_7 =acc_D+cross(omega7,cross(omega7,D_cog7_vec))+cross(alpha7,D_cog7_vec);
acc_J =acc_D+cross(omega7,cross(omega7,DJ_vec))+cross(alpha7,DJ_vec);
acc_6 =acc_J+cross(omega6,cross(omega6,J_cog6_vec))+cross(alpha6,J_cog6_vec);

acc_G =acc_I+cross(omega5,cross(omega5,IG_vec))+cross(alpha5,IG_vec);
acc_9 =acc_G+cross(omega9,cross(omega9,G_cog9_vec))+cross(alpha9,G_cog9_vec);
acc_H =acc_J+cross(omega6,cross(omega6,JH_vec))+cross(alpha5,IG_vec);
acc_10 =acc_H+cross(omega10,cross(omega10,H_cog10_vec))+cross(alpha10,H_cog10_vec);


acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);
acc_8x = acc_8(:,1);
acc_8y = acc_8(:,2);
acc_9x = acc_9(:,1);
acc_9y = acc_9(:,2);
acc_10x = acc_10(:,1);
acc_10y = acc_10(:,2);


% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
F_A_x = zeros(size(phi2));
F_A_y = zeros(size(phi2));
F_I_x = zeros(size(phi2));
F_I_y = zeros(size(phi2));
F_G_x = zeros(size(phi2));
F_G_y = zeros(size(phi2));
F_C_x = zeros(size(phi2));
F_C_y = zeros(size(phi2));
F_B_x = zeros(size(phi2));
F_B_y = zeros(size(phi2));
F_J_x = zeros(size(phi2));
F_J_y = zeros(size(phi2));
F_H_x = zeros(size(phi2));
F_H_y = zeros(size(phi2));

F_D9_x = zeros(size(phi2));
F_D9_y = zeros(size(phi2));
F_D7_x = zeros(size(phi2));
F_D7_y = zeros(size(phi2));
F_D8_x = zeros(size(phi2));
F_D8_y = zeros(size(phi2));

F_F10_x = zeros(size(phi2));
F_F10_y = zeros(size(phi2));
F_F9_x = zeros(size(phi2));
F_F9_y = zeros(size(phi2));
F_F2_x = zeros(size(phi2));
F_F2_y = zeros(size(phi2));

F_E4_x = zeros(size(phi2));
F_E4_y = zeros(size(phi2));
F_E3_x = zeros(size(phi2));
F_E3_y = zeros(size(phi2));
F_E10_x = zeros(size(phi2));
F_E10_y = zeros(size(phi2));

M_C = zeros(size(phi2));

% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
     
  A = [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
      0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;
      0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0;
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0;
      ];
    
  B = [m2*acc_2x(k);
      m2*acc_2y;
      1;
      m3*acc_3x;
      m3*acc_3y;
      1;
      m8*acc_8x;
      m8*acc_8y;
      1;
      m7*acc_7x;
      m7*acc_7y;
      1;
      m4*acc_4x;
      m4*acc_4y;
      1;
      m5*acc_5x;
      m5*acc_5y;
      1;
      m6*acc_6x;
      m6*acc_6y;
      1;
      m9*acc_9x;
      m9*acc_9y;
      1;
      m10*acc_10x;
      m10*acc_10y;
      1;
      0;
      0;
      0;
      0;
      0;
      0;];
    
    x = A\B;
    
    % save results
    F_P_x(k) = x(1);
    F_P_y(k) = x(2);
    F_Q_x(k) = x(3);
    F_Q_y(k) = x(4);
    F_R_x(k) = x(5);
    F_R_y(k) = x(6);
    F_S_x(k) = x(7);
    F_S_y(k) = x(8);
    M_P(k)   = x(9);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_P_x,F_P_y),grid
    xlabel('F_P_x [N]')
    ylabel('F_P_y [N]')
    axis tight
    subplot(222)
    plot(F_Q_x,F_Q_y),grid
    xlabel('F_Q_x [N]')
    ylabel('F_Q_y [N]')
    axis tight
    subplot(223)
    plot(F_R_x,F_R_y),grid
    xlabel('F_R_x [N]')
    ylabel('F_R_y [N]')
    axis tight
    subplot(224)
    plot(F_S_x,F_S_y),grid
    xlabel('F_S_x [N]')
    ylabel('F_S_y [N]')
    axis tight
    
    figure
    plot(t,M_P)
    ylabel('M_P [N-m]')
    xlabel('t [s]')
    
end


