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


function F=loop_closure_eqs(r11, r12, r13,r2,r3,r4,r5,r6,r7,r8,r9a, r9b,r10a, r10b, phi11, phi12, phi13,phi2, phi3_init, phi4_init, phi5_init, phi6_init, phi7_init, phi8_init, phi9_init, phi10_init) 

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi3 = phi3_init;
phi4 = phi4_init;
phi5 = phi5_init;
phi6 = phi6_init;
phi7 = phi7_init;
phi8 = phi8_init;
phi9 = phi9_init;
phi10 = phi10_init;

% loop closure equations:
F(1) = r3*cos(phi3) + r4*cos(phi4) + r5*cos(phi5) + r9b*cos(phi9) + r2*cos(phi2) - r12*cos(phi12);
F(2) = r3*sin(phi3) + r4*sin(phi4) + r5*sin(phi5) + r9b*sin(phi9) + r2*sin(phi2) - r12*sin(phi12);
F(3) = r8*cos(phi8) + r7*cos(phi7) + r6*cos(phi6) + r10b*cos(phi10) + r2*cos(phi) - r13*cos(phi13);
F(4) = r8*sin(phi8) + r7*sin(phi8) + r6*sin(phi6) + r10b*sin(phi10) + r2*sin(phi) - r13*sin(phi13);
F(5) = r3*cos(phi3) - r10a*cos(phi10) + r2*cos(phi2) -r12*cos(phi12);
F(6) = r3*sin(phi3) - r10a*sin(phi10) + r2*sin(phi2) - r12*sin(phi12);
F(7) = r8*cos(phi8) + r9a*cos(phi9) + r2*cos(phi2) - r13*cos(phi13);
F(8) = r8*sin(phi8) + r9a*sin(phi9) + r2*sin(phi2) - r13*sin(phi13);



