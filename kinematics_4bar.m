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

function [phi3,phi4,phi5, phi6, phi7, phi8,phi9,phi10,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10] = kinematics_4bar(r11, r12, r13,r2,r3,r4,r5,r6,r7,r8,r9a, r9b,r10a, r10b, phi11, phi12, phi13,phi2,dphi2, ddphi2, phi3_init, phi4_init, phi5_init, phi6_init, phi7_init, phi8_init, phi9_init, phi10_init,t,fig_kin_4bar)

% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
phi9 = zeros(size(t));
phi10 = zeros(size(t));

dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
dphi9 = zeros(size(t));
dphi10 = zeros(size(t));

ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));
ddphi9 = zeros(size(t));
ddphi10 = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off');

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi3_init phi4_init phi5_init phi6_init phi7_init phi8_init phi9_init phi10_init]',optim_options,r11, r12, r13,r2,r3,r4,r5,r6,r7,r8,r9a, r9b,r10a, r10b, phi11, phi12, phi13, phi2(k) );
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
    % save results of fsolve
    phi3(k)=x(1);
    phi4(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi7(k)=x(5);
    phi8(k)=x(6);
    phi9(k)=x(7);
    phi10(k)=x(8);
    
    
    % *** velocity analysis ***
%     
%     A = [-r3*sin(phi3(k)),  r4*sin(phi4(k));
%          r3*cos(phi3(k)), -r4*cos(phi4(k))];
%     B = [ r2*sin(phi2(k))*dphi2(k);
%          -r2*cos(phi2(k))*dphi2(k)];
       A = [-r3*sin(phi3(k)),-r4*sin(phi4(k)),-r5*sin(phi5(k)),0,0,r8*sin(phi8(k)),-(r9a+r9b)*sin(phi9(k)),0;
           r3*cos(phi3(k)),r4*cos(phi4(k)),r5*cos(phi5(k)),0,0,-r8*cos(phi8(k)),(r9a+r9b)*cos(phi9(k)),0;
           r3*sin(phi3(k)),0,0,-r6*sin(phi6(k)),-r7*sin(phi7(k)),-r8*sin(phi8(k)),0,-(r10a+r10b)*sin(phi10(k));
           -r3*cos(phi3(k)),0,0,r6*cos(phi6(k)),r7*cos(phi7(k)),r8*cos(phi8(k)),0,(r10a+r10b)*cos(phi10(k));
           -r3*sin(phi3(k)),0,0,0,0,0,0,r10a*sin(phi10(k));
           r3*cos(phi3(k)),0,0,0,0,0,0,-r10a*cos(phi10(k));
           0,0,0,0,0,-r8*sin(phi8(k)),r9a*sin(phi9(k)),0;
           0,0,0,0,0,r8*cos(phi8(k)),-r9a*cos(phi9(k)),0;];
       B= [0;
           0;
           0;
           0;
           -r2*sin(phi2(k))*dphi2(k);
           r2*cos(phi2(k))*dphi2(k);
           -r2*sin(phi2(k))*dphi2(k);
           r2*cos(phi2(k))*dphi2(k);];
       
     x = A\B;
     
     % save results
     dphi3(k) = x(1);
     dphi4(k) = x(2);
     dphi5(k) = x(3);
     dphi6(k) = x(4);
     dphi7(k) = x(5);
     dphi8(k) = x(6);
     dphi9(k) = x(7);
     dphi10(k) = x(8);
    
    
    % *** acceleration analysis ***
    
%     A = [-r3*sin(phi3(k)),  r4*sin(phi4(k));
%          r3*cos(phi3(k)), -r4*cos(phi4(k))];
%     B = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4*cos(phi4(k))*dphi4(k)^2;
%          r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4*sin(phi4(k))*dphi4(k)^2];
A=[-r3*sin(phi3(k)),-r4*sin(phi4(k)),-r5*cos(phi5(k)),0,0,r8*sin(phi8(k)),-(r9a+r9b)*sin(phi9(k)),0;
r3*cos(phi3(k)),r4*cos(phi4(k)),r5*cos(phi5(k)),0,0,-r8*cos(phi8(k)),(r9a+r9b)*cos(phi9(k)),0;
r3*sin(phi3(k)),0,0,-r6*sin(phi6(k)),-r7*sin(phi7(k)),-r8*sin(phi8(k)),0,-(r10a+r10b)*sin(phi10(k));
-r3*cos(phi3(k)),0,0,r6*cos(phi6(k)),r7*cos(phi7(k)),r8*cos(phi8(k)),0,(r10a+r10b)*cos(phi10(k));
-r3*sin(phi3(k)),0,0,0,0,0,0,r10a*sin(phi10(k));
r3*cos(phi3(k)),0,0,0,0,0,0,-r10a*cos(phi10(k));
0,0,0,0,0,-r8*sin(phi8(k)),r9a*sin(phi9(k)),0;
0,0,0,0,0,r8*cos(phi8(k)),-r9a*sin(phi9(k)),0;];

B=[r3*cos(phi3(k))*dphi3(k)^2+r4*cos(phi4(k))*dphi4(k)^2+r5*cos(phi5(k))*dphi5(k)^2-r8*cos(phi8(k))*dphi8(k)^2+(r9a+r9b)*cos(phi9(k))*dphi9(k)^2;
r3*sin(phi3(k))*dphi3(k)^2+r4*sin(phi4(k))*dphi4(k)^2+r5*sin(phi5(k))*dphi5(k)^2-r8*sin(phi8(k))*dphi8(k)^2+(r9a+r9b)*sin(phi9(k))*dphi9(k)^2;
-r3*cos(phi3(k))*dphi3(k)^2+r6*cos(phi6(k))*dphi6(k)^2+r7*cos(phi7(k))*dphi7(k)^2+r8*cos(phi8(k))*dphi8(k)^2+(r10a+r10b)*cos(phi10(k))*dphi10(k)^2;
-r3*sin(phi3(k))*dphi3(k)^2+r6*sin(phi6(k))*dphi6(k)^2+r7*sin(phi7(k))*dphi7(k)^2+r8*sin(phi8(k))*dphi8(k)^2+(r10a+r10b)*sin(phi10(k))*dphi10(k)^2;
r3*cos(phi3(k))*dphi3(k)^2-r10a*cos(phi10(k))*dphi10(k)^2-r2*sin(phi2(k))*ddphi2(k)-r2*cos(phi2(k))*dphi2(k)^2;
r3*sin(phi3(k))*dphi3(k)^2-r10a*sin(phi10(k))*dphi10(k)^2+r2*cos(phi2(k))*ddphi2(k)-r2*sin(phi2(k))*dphi2(k)^2;
r8*cos(phi8(k))*dphi8(k)^2-r9a*cos(phi9(k))*dphi9(k)^2-r2*sin(phi2(k))*ddphi2(k)-r2*cos(phi2(k))*dphi2(k)^2;
r8*sin(phi8(k))*dphi8(k)^2-r9a*sin(phi9(k))*dphi9(k)^2+r2*cos(phi2(k))*ddphi2(k)-r2*sin(phi2(k))*dphi2(k)^2;];
%      A = [-r3*sin(phi3(k))  , -r4*sin(phi4(k))  , -r5*cos(phi5(k))  , 0                 , 0                 , r8*sin(phi8(k))   , -(r9a+r9b)*sin(phi9(k))   , 0;
%           r3*cos(phi3(k))   , r4*cos(phi4(k))   , r5*cos(phi5(k))   , 0                 , 0                 , -r8*cos(phi8(k))  , (r9a + r9b)*cos(phi9(k))  , 0;
%           r3*sin(phi3(k))   , 0                 , 0                 , -r6*sin(phi6(k))  , -r7*sin(phi7(k))  , -r8*sin(phi8(k))  , 0                         , -(r10a+r10b)*sin(phi10(k));
%           -r3*cos(phi3(k))  , 0                 , 0                 , r6*cos(phi6(k))   , r7*cos(phi7(k))   , r8*cos(phi8(k))   , 0                         , (r10a+r10b)*cos(phi10(k));
%           -r3*sin(phi3(k))  , 0                 , 0                 , 0                 , 0                 , 0                 , 0                         , r10a*sin(phi10(k));
%           r3*cos(phi3(k))   , 0                 , 0                 , 0                 , 0                 , 0                 , 0                         , -r10a*cos(phi10(k));
%           0                 , 0                 , 0                 , 0                 , 0                 , -r8*sin(phi8(k))  , r9a*sin(phi9(k))          ,0;
%           0                 , 0                 , 0                 , 0                 , 0                 , r8*cos(phi8(k))   , -r9a*sin(phi9(k))         , 0;];
     
%      B = [r3*cos(phi3(k))*dphi3(k)^2 + r4*cos(phi4(k))*dphi4(k)^2 + r5*cos(phi5(k))*dphi5(k)^2 -r8*cos(phi8(k))*dphi8(k)^2 + (r9a+r9b)*cos(phi9(k))*dphi9(k)^2;
%          r3*sin(phi3(k))*dphi3(k)^2 +r4*sin(phi4(k))*dphi4(k)^2 + r5*sin(phi5(k))*dphi5(k)^2 - r8*sin(phi8(k))*dphi8(k)^2 + (r9a+r9b)*sin(phi9(k))*dphi9(k)^2;
%          -r3*cos(phi3(k))*dphi3(k)^2 + r6*cos(phi6(k))*dphi6(k)^2 +r7*cos(phi7(k))*dphi7(k)^2 +r8*cos(phi8(k))*dphi8(k)^2 + (r10a+r10b)*cos(phi10(k))*dphi10(k)^2;
%          -r3*sin(phi3(k))*dphi3(k)^2 + r6*sin(phi6(k))*dphi6(k)^2 +r7*sin(phi7(k))*dphi7(k)^2 +r8*sin(phi8(k))*dphi8(k)^2 + (r10a+r10b)*sin(phi10(k))*dphi10(k)^2;
%          r3*cos(phi3(k))*dphi3(k)^2 - r10a*cos(phi10(k))*dphi10(k)^2 - r2*sin(phi2(k))*ddphi2(k) -r2*cos(phi2(k))*dphi2(k)^2;
%          r3*sin(phi3(k))*dphi3(k)^2 - r10a*sin(phi10(k))*dphi10(k)^2 +r2*cos(phi2(k))*ddphi2(k) - r2*sin(phi2(k))*dphi2(k)^2;
%          r8*cos(phi8(k))*dphi8(k)^2 - r9a*cos(phi9(k))*dphi9(k)^2 - r2*sin(phi2(k))*ddphi2(k) -r2*cos(phi2(k))*dphi2(k)^2;
%          r8*sin(phi8(k))*dphi8(k)^2 - r9a*sin(phi9(k))*dphi9(k)^2 +r2*cos(phi2(k))*ddphi2(k) - r2*sin(phi2(k))*dphi2(k)^2;];
    x = A\B;
    % save results
    ddphi3(k) = x(1);
    ddphi4(k) = x(2);
    ddphi5(k) = x(3);
    ddphi6(k) = x(4);
    ddphi7(k) = x(5);
    ddphi8(k) = x(6);
    ddphi9(k) = x(7);
    ddphi10(k) = x(8);
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    phi5_init = phi5(k)+Ts*dphi5(k);
    phi6_init = phi6(k)+Ts*dphi6(k);
    phi7_init = phi7(k)+Ts*dphi7(k);
    phi8_init = phi8(k)+Ts*dphi8(k);
    phi9_init = phi9(k)+Ts*dphi9(k);
    phi10_init = phi10(k)+Ts*dphi10(k);

end % loop over positions


%disp(phi3);%% toont het verloop van phi 3 voor phi2 tussen 0 en omega*t_end






% *** create movie ***

% point P = fixed
A = 0;
% point S = fixed
C = r12*exp(j*phi12);
% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -2*r2;
y_bottom = -2*max(r2,r4);
x_right = r11+2*r4;
y_top = 2*max(r2,r4);

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    F = C + r2*exp(j*phi2(index));
    E1 = A + r3*exp(j*phi3(index));
    E2 = F + r10a * exp(j*phi10(index));
    
    loop1 = [A E1 E2 F C];
    
    I = E1 + r4*exp(j*phi4(index));
    G = I + r5*exp(j*phi5(index));
    loop2 = [E1 I G F ];
    
    B = A + r11 * exp(j*0);
    D = B + r8*exp(j*phi8(index));
    loop3 = [B D F C ];
    
    J = D + r7*exp(j*phi7(index));
    H = J + r6*exp(j*phi6(index));
    loop4 = [D J H F];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),real(loop2),imag(loop2),real(loop3),imag(loop3),real(loop4),imag(loop4),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
close(10)


% *** plot figures ***

if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    A = 0;
    C = r12*exp(j*phi12(index));
    F = C + r2*exp(j*phi2(index));
    E = F + r10a*exp(j*phi10(index));
    
%     P = 0;
%     S = r1*exp(j*phi1);
%     Q = P + r2 * exp(j*phi2(index));
%     R = Q + r3 * exp(j*phi3(index));
    
    figure
    assembly=[A, C, F, E];
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(311)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(312)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(313)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(312)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(313)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(311)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(312)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(313)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    xlabel('t [s]')
end
% 
% 
% 
