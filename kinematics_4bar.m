

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
% optim_options = optimset('Display','off');
optim_options = optimset('Display','off','TolFun',10^-12);
% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    
%     if(phi2(k) ~= 2*pi) 
%         disp(k*Ts);
%     end
    


    % *** position analysis ***
    
    % fsolve slves the non-linear set of equations
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
   

    %%AANGEPAST, FOUT OKE NU
    A=[-r3*sin(phi3(k)),-r4*sin(phi4(k)),-r5*sin(phi5(k)),0,0,r8*sin(phi8(k)),-(r9a+r9b)*sin(phi9(k)),0;
    r3*cos(phi3(k)),r4*cos(phi4(k)),r5*cos(phi5(k)),0,0,-r8*cos(phi8(k)),(r9a+r9b)*cos(phi9(k)),0;
    r3*sin(phi3(k)),0,0,-r6*sin(phi6(k)),-r7*sin(phi7(k)),-r8*sin(phi8(k)),0,-(r10a+r10b)*sin(phi10(k));
    -r3*cos(phi3(k)),0,0,r6*cos(phi6(k)),r7*cos(phi7(k)),r8*cos(phi8(k)),0,(r10a+r10b)*cos(phi10(k));
    -r3*sin(phi3(k)),0,0,0,0,0,0,r10a*sin(phi10(k));
    r3*cos(phi3(k)),0,0,0,0,0,0,-r10a*cos(phi10(k));
    0,0,0,0,0,-r8*sin(phi8(k)),r9a*sin(phi9(k)),0;
    0,0,0,0,0,r8*cos(phi8(k)),-r9a*cos(phi9(k)),0;];

    B=[r3*cos(phi3(k))*dphi3(k)^2+r4*cos(phi4(k))*dphi4(k)^2+r5*cos(phi5(k))*dphi5(k)^2-r8*cos(phi8(k))*dphi8(k)^2+(r9a+r9b)*cos(phi9(k))*dphi9(k)^2;
    r3*sin(phi3(k))*dphi3(k)^2+r4*sin(phi4(k))*dphi4(k)^2+r5*sin(phi5(k))*dphi5(k)^2-r8*sin(phi8(k))*dphi8(k)^2+(r9a+r9b)*sin(phi9(k))*dphi9(k)^2;
    -r3*cos(phi3(k))*dphi3(k)^2+r6*cos(phi6(k))*dphi6(k)^2+r7*cos(phi7(k))*dphi7(k)^2+r8*cos(phi8(k))*dphi8(k)^2+(r10a+r10b)*cos(phi10(k))*dphi10(k)^2;
    -r3*sin(phi3(k))*dphi3(k)^2+r6*sin(phi6(k))*dphi6(k)^2+r7*sin(phi7(k))*dphi7(k)^2+r8*sin(phi8(k))*dphi8(k)^2+(r10a+r10b)*sin(phi10(k))*dphi10(k)^2;
    -r2*sin(phi2(k))*ddphi2(k)-r2*cos(phi2(k))*dphi2(k)^2+r3*cos(phi3(k))*dphi3(k)^2-r10a*cos(phi10(k))*dphi10(k)^2;
    r2*cos(phi2(k))*ddphi2(k)-r2*sin(phi2(k))*dphi2(k)^2+r3*sin(phi3(k))*dphi3(k)^2-r10a*sin(phi10(k))*dphi10(k)^2;
    -r2*sin(phi2(k))*ddphi2(k)-r2*cos(phi2(k))*dphi2(k)^2+r8*cos(phi8(k))*dphi8(k)^2-r9a*cos(phi9(k))*dphi9(k)^2;
    +r2*cos(phi2(k))*ddphi2(k)-r2*sin(phi2(k))*dphi2(k)^2+r8*sin(phi8(k))*dphi8(k)^2-r9a*sin(phi9(k))*dphi9(k)^2;];

  
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
F_diff1=[];
F_diff2=[];
F_diff3=[];
F_diff4=[];
for m=1:length(index_vec)  
    index = index_vec(m);
    F1 = C + r2*exp(j*phi2(index));
    E1 = A + r3*exp(j*phi3(index));
    E2 = F1 + r10a * exp(j*phi10(index));
    
%     E_diff= [E_diff sqrt(real(E1-E2)^2 + imag(E1-E2)^2)];
%     V_E1= (r3*dphi3(index))-
    
    loop1 = [A E1 E2 F1 C];
    
    I = E1 + r4*exp(j*phi4(index));
    G = I + r5*exp(j*phi5(index));
    loop2 = [E1 I G F1 ];
    
    B = A + r11 * exp(j*0);
    D = B + r8*exp(j*phi8(index));
    loop3 = [B D F1 C ];
    
    J = D + r7*exp(j*phi7(index));
    H = J + r6*exp(j*phi6(index));
    loop4 = [D J H F1];
    
    C2 = B + r13*exp(j*phi13);
    loop5 = [A C C2 B];
    
    F2= E1 - r10a*exp(j*phi10(index));
    F3= G + r9b*exp(j*phi9(index));
    F4= D - r9a*exp(j*phi9(index));
    F5= H+ r10b*exp(j*phi10(index));
    
    F_diff1= [F_diff1 sqrt(real(F1-F2)^2 + imag(F1-F2)^2)];
    F_diff2= [F_diff2 sqrt(real(F1-F3)^2 + imag(F1-F3)^2)];
    F_diff3= [F_diff3 sqrt(real(F1-F4)^2 + imag(F1-F4)^2)];
    F_diff4= [F_diff4 sqrt(real(F1-F5)^2 + imag(F1-F5)^2)];
   
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),real(loop2),imag(loop2),real(loop3),imag(loop3),real(loop4),imag(loop4),real(loop5),imag(loop5),'-o')
    
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
    
    C2 = B + r13*exp(j*phi13);
    loop5 = [A C C2 B];
    
   
    
    
    
    
    
    
    figure
plot(real(loop1),imag(loop1),real(loop2),imag(loop2),real(loop3),imag(loop3),real(loop4),imag(loop4),real(loop5),imag(loop5),'-o')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(331)
    plot(t,F_diff1);
    xlabel('t [s]')
    ylabel('F_{diff,1} [rad]')
    subplot(332)
    plot(t,F_diff2);
    xlabel('t [s]')
    ylabel('F_{diff,2} [rad]')
    subplot(333)
    plot(t,F_diff3);
    xlabel('t [s]')
    ylabel('F_{diff,3} [rad]')
    
    figure
    plot(t,F_diff4);
    xlabel('t [s]')
    ylabel('F_{diff,4} [rad]')

    
    %%Position
    figure
    subplot(331)
    plot(t,phi2)
    xlabel('t [s]')
    ylabel('\phi_2 [rad]')
    subplot(332)
    plot(t,phi3)
    xlabel('t [s]')
    ylabel('\phi_3 [rad]')
    subplot(333)
    plot(t,phi4)
    xlabel('t [s]')
    ylabel('\phi_4 [rad]')
    subplot(334)
    plot(t,phi5)
    xlabel('t [s]')
    ylabel('\phi_5 [rad]')
    subplot(335)
    plot(t,phi6)
    xlabel('t [s]')
    ylabel('\phi_6 [rad]')
    subplot(336)
    plot(t,phi7)
    xlabel('t [s]')
    ylabel('\phi_7 [rad]')
    subplot(337)
    plot(t,phi8)
    xlabel('t [s]')
    ylabel('\phi_8 [rad]')
    subplot(338)
    plot(t,phi9)
    xlabel('t [s]')
    ylabel('\phi_9 [rad]')
    subplot(339)
    plot(t,phi10)
    ylabel('\phi_10 [rad]')
    xlabel('t [s]')
    
    
 
    %%Velocity
    figure
    subplot(331)
    plot(t,dphi2)
    xlabel('t [s]')
    ylabel('d\phi_2 [rad/s]')
    subplot(332)
    plot(t,dphi3)
    xlabel('t [s]')
    ylabel('d\phi_3 [rad/s]')
    subplot(333)
    plot(t,dphi4)
    xlabel('t [s]')
    ylabel('d\phi_4 [rad/s]')
    subplot(334)
    plot(t,dphi5)
    xlabel('t [s]')
    ylabel('d\phi_5 [rad/s]')
    subplot(335)
    plot(t,dphi6)
    xlabel('t [s]')
    ylabel('d\phi_6 [rad/s]')
    subplot(336)
    plot(t,dphi7)
    xlabel('t [s]')
    ylabel('d\phi_7 [rad/s]')
    subplot(337)
    plot(t,dphi8)
    xlabel('t [s]')
    ylabel('d\phi_8 [rad/s]')
    subplot(338)
    plot(t,dphi9)
    xlabel('t [s]')
    ylabel('d\phi_9 [rad/s]')
    subplot(339)
    plot(t,dphi10)
    ylabel('d\phi_10 [rad/s]')
    xlabel('t [s]')
    
    
    
    %%Acceleration
    figure
    subplot(331)
    plot(t,ddphi2)
    xlabel('t [s]')
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(332)
    plot(t,ddphi3)
    xlabel('t [s]')
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(333)
    plot(t,ddphi4)
    xlabel('t [s]')
    ylabel('dd\phi_4 [rad/s^2]')
    subplot(334)
    plot(t,ddphi5)
    xlabel('t [s]')
    ylabel('dd\phi_5 [rad/s^2]')
    subplot(335)
    plot(t,ddphi6)
    xlabel('t [s]')
    ylabel('dd\phi_6 [rad/s^2]')
    subplot(336)
    plot(t,ddphi7)
    xlabel('t [s]')
    ylabel('dd\phi_7 [rad/s^2]')
    subplot(337)
    plot(t,ddphi8)
    xlabel('t [s]')
    ylabel('dd\phi_8 [rad/s^2]')
    subplot(338)
    plot(t,ddphi9)
    xlabel('t [s]')
    ylabel('dd\phi_9 [rad/s^2]')
    subplot(339)
    plot(t,ddphi10)
    ylabel('dd\phi_10 [rad/s^2]')
    xlabel('t [s]')
end

