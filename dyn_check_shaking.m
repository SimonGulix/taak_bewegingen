
function dyn_check_shaking(vel_2x,vel_2y,vel_3x,vel_3y,vel_4x,vel_4y,vel_5x,vel_5y,vel_6x,vel_6y,vel_7x,vel_7y,vel_8x,vel_8y,vel_9x,vel_9y,vel_10x,vel_10y,...
                   acc_2x,acc_2y,acc_3x,acc_3y,acc_4x,acc_4y,acc_5x,acc_5y,acc_6x,acc_6y,acc_7x,acc_7y,acc_8x,acc_8y,acc_9x,acc_9y,acc_10x,acc_10y,...
                   M_C,F_A_x,F_A_y,F_C_x,F_C_y,F_B_x,F_B_y,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10,m2,m3,m4,m5,m6,m7,m8,m9,m10,...
                   J2,J3,J4,J5,J6,J7,J8,J9,J10,CF_vec,AE_vec,EI_vec,IG_vec,JH_vec,DJ_vec,BD_vec,GD_vec,HE_vec,CA_x,CA_y,CB_x,CB_y,t,fig_dyn_check_shaking,g)
 %%FORCES 
 F_shake_x= -(m2*acc_2x+m3*acc_3x+m4*acc_4x+m5*acc_5x+m6*acc_6x+m7*acc_7x+m8*acc_8x+m9*acc_9x+m10*acc_10x);   
 F_shake_y= -(m2*acc_2y+m3*acc_3y+m4*acc_4y+m5*acc_5y+m6*acc_6y+m7*acc_7y+m8*acc_8y+m9*acc_9y+m10*acc_10y);
             
 F_external_x= F_A_x+F_C_x+F_B_x;
 F_external_y= F_A_y+F_C_y+F_B_y;
 
 F_err_x=F_external_x+F_shake_x;
 F_err_y=F_external_y+F_shake_y;
 
 
 %%DRIVING MOMENT
 acc_2=[acc_2x acc_2y zeros(size(acc_2x))];
 acc_3=[acc_3x acc_3y zeros(size(acc_3x))];
 acc_4=[acc_4x acc_4y zeros(size(acc_4x))];
 acc_5=[acc_5x acc_5y zeros(size(acc_5x))];
 acc_6=[acc_6x acc_6y zeros(size(acc_6x))];
 acc_7=[acc_7x acc_7y zeros(size(acc_7x))];
 acc_8=[acc_8x acc_8y zeros(size(acc_8x))];
 acc_9=[acc_9x acc_9y zeros(size(acc_9x))];
 acc_10=[acc_10x acc_10y zeros(size(acc_10x))];
 
 CA_vec=[CA_x*ones(size(acc_2x)) CA_y*ones(size(acc_2x)) zeros(size(acc_2x))];
 CB_vec=[CB_x*ones(size(acc_2x)) CB_y*ones(size(acc_2x)) zeros(size(acc_2x))];
 
 C_vec_2=CF_vec/2;
 C_vec_3=CA_vec+AE_vec/2;
 C_vec_4=CA_vec+AE_vec+EI_vec/2;
 C_vec_5=CA_vec+AE_vec+EI_vec+IG_vec/2;
 C_vec_6=CB_vec+BD_vec+DJ_vec+JH_vec/2;
 C_vec_7=CB_vec+BD_vec+DJ_vec/2;
 C_vec_8=CB_vec+BD_vec/2;
 C_vec_9=CA_vec+AE_vec+EI_vec+IG_vec+GD_vec/2;
 C_vec_10=CB_vec+BD_vec+DJ_vec+JH_vec+HE_vec/2;

 
 M_shake=-(J2*alpha2(:,3)+J3*alpha3(:,3)+J4*alpha4(:,3)+J5*alpha5(:,3)+J6*alpha6(:,3)+J7*alpha7(:,3)+J8*alpha8(:,3)+J9*alpha9(:,3)+J10*alpha10(:,3)+...
            m2*cross(C_vec_2,acc_2)+m3*cross(C_vec_3,acc_3)+m4*cross(C_vec_4,acc_4)+m5*cross(C_vec_5,acc_5)+m6*cross(C_vec_6,acc_6)+m7*cross(C_vec_7,acc_7)+...
            m8*cross(C_vec_8,acc_8)+m9*cross(C_vec_9,acc_9)+m10*cross(C_vec_10,acc_10));
 
 M_external=M_C+ cross(CA_vec,[F_A_x F_A_y zeros(size(F_A_x))])+cross(CB_vec,[F_B_x F_B_y zeros(size(F_B_x))]);
 M_error=M_shake(:,3)+M_external(:,3);
 
    if fig_dyn_check_shaking

      figure
      subplot(311)
      plot(t,F_external_x),grid
      xlabel('t [s]')        
      ylabel('F_{external_x} [N]')
      subplot(312)
      plot(t,-F_shake_x),grid
      xlabel('t [s]')
      ylabel('F_{shake,x} [N]')
      subplot(313)
      plot(t,F_err_x),grid
      xlabel('t [s]')
      ylabel('F_{error,x} [N]')

      figure
      subplot(311)
      plot(t,F_external_y),grid
      xlabel('t [s]')        
      ylabel('F_{external_y} [N]')
      subplot(312)
      plot(t,-F_shake_y),grid
      xlabel('t [s]')
      ylabel('F_{shake,y} [N]')
      subplot(313)
      plot(t,F_err_y),grid
      xlabel('t [s]')
      ylabel('F_{error,y} [N]')
      
      
      figure
      subplot(311)
      plot(t,M_external(:,3)),grid
      xlabel('t [s]')        
      ylabel('M_{external} [Nm]')
      subplot(312)
      plot(t,-M_shake(:,3)),grid
      xlabel('t [s]')
      ylabel('M_{shake} [Nm]')
      subplot(313)
      plot(t,M_error),grid
      xlabel('t [s]')
      ylabel('M_{error} [Nm]')
    end               
end