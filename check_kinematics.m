
%%% TO CALCULATE ERRORS

function [ diffphi3,diffphi4,diffphi5,diffphi6,diffphi7,diffphi8,diffphi9,diffphi10,...
           ddiffphi3,ddiffphi4,ddiffphi5,ddiffphi6,ddiffphi7,ddiffphi8,ddiffphi9,ddiffphi10 ] = ...
check_kinematics( phi3,phi4,phi5,phi6,phi7,phi8,phi9,phi10, ...
           dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,dphi9,dphi10, ...
           ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,ddphi9,ddphi10, ...
           Ts, fig_kin_check, t)
       
% kin_check verifies kinematica

% Central difference approximation of derivative:
% f'(x) = (f(x+1)-f(x-1))/(2*Ts)

diffphi3  = (phi3(3:size(phi3)) - phi3(1:size(phi3)-2)) / (2*Ts);
diffphi4  = (phi4(3:size(phi3)) - phi4(1:size(phi3)-2)) / (2*Ts);
diffphi5  = (phi5(3:size(phi3)) - phi5(1:size(phi3)-2)) / (2*Ts);
diffphi6  = (phi6(3:size(phi3)) - phi6(1:size(phi3)-2)) / (2*Ts);
diffphi7  = (phi7(3:size(phi3)) - phi7(1:size(phi3)-2)) / (2*Ts);
diffphi8  = (phi8(3:size(phi3)) - phi8(1:size(phi3)-2)) / (2*Ts);
diffphi9  = (phi9(3:size(phi3)) - phi9(1:size(phi3)-2)) / (2*Ts);
diffphi10  = (phi10(3:size(phi3)) - phi10(1:size(phi3)-2)) / (2*Ts);
% Central difference approximation of second derivative:
% f"(x) = (f(x+1)+f(x-1)-2*f(x))/(Ts^2)

ddiffphi3  = (phi3(3:size(phi3)) + phi3(1:size(phi3)-2) - 2*phi3(2:size(phi3)-1)) / (Ts^2);
ddiffphi4  = (phi4(3:size(phi3)) + phi4(1:size(phi3)-2) - 2*phi4(2:size(phi3)-1)) / (Ts^2);
ddiffphi5  = (phi5(3:size(phi3)) + phi5(1:size(phi3)-2) - 2*phi5(2:size(phi3)-1)) / (Ts^2);
ddiffphi6  = (phi6(3:size(phi3)) + phi6(1:size(phi3)-2) - 2*phi6(2:size(phi3)-1)) / (Ts^2);
ddiffphi7  = (phi7(3:size(phi3)) + phi7(1:size(phi3)-2) - 2*phi7(2:size(phi3)-1)) / (Ts^2);
ddiffphi8  = (phi8(3:size(phi3)) + phi8(1:size(phi3)-2) - 2*phi8(2:size(phi3)-1)) / (Ts^2);
ddiffphi9  = (phi9(3:size(phi3)) + phi9(1:size(phi3)-2) - 2*phi9(2:size(phi3)-1)) / (Ts^2);
ddiffphi10  = (phi10(3:size(phi3)) + phi10(1:size(phi3)-2) - 2*phi10(2:size(phi3)-1)) / (Ts^2);

% Calculate the errors

erdphi3  = diffphi3 - dphi3(2:size(phi3)-1,:);
erdphi4  = diffphi4 - dphi4(2:size(phi3)-1,:);
erdphi5  = diffphi5 - dphi5(2:size(phi3)-1,:);
erdphi6  = diffphi6 - dphi6(2:size(phi3)-1,:);
erdphi7  = diffphi7 - dphi7(2:size(phi3)-1,:);
erdphi8  = diffphi8 - dphi8(2:size(phi3)-1,:);
erdphi9  = diffphi9 - dphi9(2:size(phi3)-1,:);
erdphi10  = diffphi10 - dphi10(2:size(phi3)-1,:);

erddphi3  = ddiffphi3 - ddphi3(2:size(phi3)-1,:);
erddphi4  = ddiffphi4 - ddphi4(2:size(phi3)-1,:);
erddphi5  = ddiffphi5 - ddphi5(2:size(phi3)-1,:);
erddphi6  = ddiffphi6 - ddphi6(2:size(phi3)-1,:);
erddphi7  = ddiffphi7 - ddphi7(2:size(phi3)-1,:);
erddphi8  = ddiffphi8 - ddphi8(2:size(phi3)-1,:);
erddphi9  = ddiffphi9 - ddphi9(2:size(phi3)-1,:);
erddphi10  = ddiffphi10 - ddphi10(2:size(phi3)-1,:);


% Plots all possible functions when fig_kin_check = 1

if fig_kin_check
    
figure

subplot(331)
    plot(t(1:size(erdphi3),:),diffphi3)
    ylabel('d\phi_3 approximation [rad]')
    xlabel('t [s]')
subplot(332)
    plot(t,dphi3)
    ylabel('d\phi_3 exact [rad]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(erdphi3),:),erdphi3)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(erdphi4),:),diffphi4)
    ylabel('d\phi_4 approximation [rad]')
    xlabel('t [s]')
subplot(335)
    plot(t,dphi4)
    ylabel('d\phi_4 exact [rad]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(erdphi4),:),erdphi4)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(erdphi5),:),diffphi5)
    ylabel('d\phi_5 approximation [rad]')
    xlabel('t [s]')
subplot(338)
    plot(t,dphi5)
    ylabel('d\phi_5 exact [rad]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(erdphi5),:),erdphi5)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(erdphi6),:),diffphi6)
    ylabel('d\phi_6 approximation [rad]')
    xlabel('t [s]')
subplot(332)
    plot(t,dphi6)
    ylabel('d\phi_6 exact [rad]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(erdphi6),:),erdphi6)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(erdphi7),:),diffphi7)
    ylabel('d\phi_7 approximation [rad]')
    xlabel('t [s]')
subplot(335)
    plot(t,dphi7)
    ylabel('d\phi_7 exact [rad]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(erdphi7),:),erdphi7)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(erdphi8),:),diffphi8)
    ylabel('d\phi_8 approximation [rad]')
    xlabel('t [s]')
subplot(338)
    plot(t,dphi8)
    ylabel('d\phi_8 exact [rad]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(erdphi8),:),erdphi8)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
    
figure 

subplot(340)
    plot(t(1:size(erdphi9),:),diffphi9)
    ylabel('d\phi_9 approximation [rad]')
    xlabel('t [s]')
subplot(341)
    plot(t,dphi9)
    ylabel('d\phi_9 exact [rad]')
    xlabel('t [s]')
subplot(342)
    plot(t(1:size(erdphi9),:),erdphi9)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(343)
    plot(t(1:size(erdphi10),:),diffphi10)
    ylabel('d\phi_10 approximation [rad]')
    xlabel('t [s]')
subplot(344)
    plot(t,dphi10)
    ylabel('d\phi_10 exact [rad]')
    xlabel('t [s]')
subplot(345)
    plot(t(1:size(erdphi10),:),erdphi10)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(erddphi3),:),ddiffphi3)
    ylabel('dd\phi_3 approximation [rad]')
    xlabel('t [s]')
subplot(332)
    plot(t,ddphi3)
    ylabel('dd\phi_3 exact [rad]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(erddphi3),:),erddphi3)
    axis([0 10 -0.2 0.2])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(erddphi4),:),ddiffphi4)
    ylabel('dd\phi_4 approximation [rad]')
    xlabel('t [s]')
subplot(335)
    plot(t,ddphi4)
    ylabel('dd\phi_4 exact [rad]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(erddphi4),:),erddphi4)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(erddphi5),:),ddiffphi5)
    ylabel('dd\phi_5 approximation [rad]')
    xlabel('t [s]')
subplot(338)
    plot(t,ddphi5)
    ylabel('dd\phi_5 exact [rad]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(erddphi5),:),erddphi5)
    axis([0 10 -0.3 0.3])
    ylabel('error [rad]')
    xlabel('t [s]')
    
figure

subplot(331)
    plot(t(1:size(erddphi6),:),ddiffphi6)
    ylabel('dd\phi_6 approximation [rad]')
    xlabel('t [s]')
subplot(332)
    plot(t,ddphi6)
    ylabel('dd\phi_6 exact [rad]')
    xlabel('t [s]')
subplot(333)
    plot(t(1:size(erddphi6),:),erddphi6)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(334)
    plot(t(1:size(erddphi7),:),ddiffphi7)
    ylabel('dd\phi_7 approximation [rad]')
    xlabel('t [s]')
subplot(335)
    plot(t,ddphi7)
    ylabel('dd\phi_7 exact [rad]')
    xlabel('t [s]')
subplot(336)
    plot(t(1:size(erddphi7),:),erddphi7)
    axis([0 10 -0.3 0.3])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(337)
    plot(t(1:size(erddphi8),:),ddiffphi8)
    ylabel('dd\phi_8 approximation [rad]')
    xlabel('t [s]')
subplot(338)
    plot(t,ddphi8)
    ylabel('dd\phi_8 exact [rad]')
    xlabel('t [s]')
subplot(339)
    plot(t(1:size(erddphi8),:),erddphi8)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
    
figure
subplot(340)
    plot(t(1:size(erddphi9),:),ddiffphi9)
    ylabel('dd\phi_9 approximation [rad]')
    xlabel('t [s]')
subplot(341)
    plot(t,ddphi9)
    ylabel('dd\phi_9 exact [rad]')
    xlabel('t [s]')
subplot(342)
    plot(t(1:size(erddphi9),:),erddphi9)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')
subplot(343)
    plot(t(1:size(erddphi10),:),ddiffphi10)
    ylabel('dd\phi_10 approximation [rad]')
    xlabel('t [s]')
subplot(344)
    plot(t,ddphi10)
    ylabel('dd\phi_10 exact [rad]')
    xlabel('t [s]')
subplot(345)
    plot(t(1:size(erddphi10),:),erddphi10)
    axis([0 10 -0.15 0.15])
    ylabel('error [rad]')
    xlabel('t [s]')

end

