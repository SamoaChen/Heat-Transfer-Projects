clc;
clear all;
close all;
TH=400;
TL=100;
R1=0.002/(17*0.3*0.5);
R2=0.003/(372*0.3*0.5);
R3=0.002/(17*0.3*0.5);
C1=400*0.002*0.3*0.5*8000;
C2=385*0.003*0.3*0.5*8954;
C3=400*0.002*0.3*0.5*8000;

s=tf('s');
% syms s;
% 
% t_final=7;
% dt=0.01;
% time=(0:dt:t_final).';
% 
% A=[(R1/2)+(1/(s*C1)),-1/(s*C1),0,0;
%     -1/(s*C1),(1/(s*C1))+(R1/2)+(R2/2)+(1/(s*C2)),-1/(s*C2),0;
%     0,-1/(s*C2),(1/(s*C2))+R2/2+R3/2+(1/(s*C3)),-1/(s*C3);
%     0,0,-1/(s*C3),(1/(s*C3))+R3/2];
% B=[0,0,0,(TH-TL)/s].';
% x_s=A\B;
% 
% [x,t] = impulse(x_s,time);
% figure(2)
% plot(time,x);
% %%
% t_final=7;
% dt=0.01;
% t=(0:dt:t_final).';
% 
% T1=(TH/s)-x_s(1)*R1/2;
% T1_t=impulse(T1,t);
% T2=T1-x_s(2)*R1/2;
% T2_t=impulse(T2,t);
% T3=T2-x_s(2)*R2/2;
% T3_t=impulse(T3,t);
% T4=T3-x_s(3)*R2/2;
% T4_t=impulse(T4,t);
% T5=T4-x_s(3)*R3/2;
% T5_t=impulse(T5,t);
% 
% 
figure(5)
plot(t,T1_t,'LineWidth',1.1)
hold on
plot(t,T2_t,'LineWidth',1.1)
hold on
plot(t,T3_t,'LineWidth',1.1)
hold on
plot(t,T4_t,'LineWidth',1.1)
hold on
plot(t,T5_t,'LineWidth',1.1)
title('Heat Transfer Inside the Wall','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Heat Transfer (J/s)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'q1','q2','q3','q4','q5'},'Location','southeast')
%     
% 
% 
C2=385*0.002*0.3*0.5*8954;

% % part1 q2
A1=[1/(R1/2)+1/(R1/2)+(s*C1),0,0;
    0,4/R2+(s*C2),0;
    0,0,4/R3+(s*C3)];
B1=[(TH/s)/(R1/2)+(TH/s)*(s*C1)+(TL/s)/(R1/2),(TH/s)/(R2/2)+(TH/s)*(s*C2)+(TL/s)/(R2/2),(TH/s)/(R3/2)+(TH/s)*(s*C3)+(TL/s)/(R3/2)].';
x_s1=A1\B1;

t_final=1;
dt=0.01;
t=(0:dt:t_final).';
x_t1=impulse(x_s1,t);

figure(2)
plot(t,x_t1,'LineWidth',1.1)
title('Heat Transfer Inside the Wall','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Heat Transfer (J/s)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'q1','q2','q3'},'Location','southeast')   

% %% part 2
% R=1;
% C=0.45*200;
% 
% B2=((pi*5/12)/(s^2+(pi*1/60)^2))+(25/s)-(300/s);
% A2=R+(1/(s*C));
% x_s3=A2\B2;
% % x_t3=ilaplace(x_s3);
% 
% t_final=1000;
% dt=0.01;
% t=(0:dt:t_final).';
% x_t3=impulse(x_s3,t);
% 
% % Verifying the calculated laplace domain driven function
% Df=impulse((pi*5/12)/(s^2+(pi*1/60)^2)+(25/s),t);
% figure(5)
% plot(t,Df,'LineWidth',1.5)
% title('Driven Temperature','FontSize',14)
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Temperature (K)','FontSize',12)
%     set(gca,'FontSize',10)
%     grid on
%     set(gca,'GridAlpha',0.3)
    
% % Ex#2 part b)   
% figure(1)
% plot(t,x_t3,'LineWidth',1.5)
% title('Heat Transfer of Iron Block','FontSize',14)
% xlabel('Time (s)','FontSize',12)
% ylabel('Heat Transfer Rate (J/s))','FontSize',12)
% set(gca,'FontSize',10)
% grid on
% set(gca,'GridAlpha',0.3)
% 
% T=((pi*5/12)/(s^2+(pi*1/60)^2))+(25/s)-R*x_s3+(273.15/s);
% 
% figure(2)
% T_t3=impulse(T,t);
% plot(t,T_t3,'LineWidth',1.5)
% title('Temperature of Iron Block','FontSize',14)
% xlabel('Time (s)','FontSize',12)
% ylabel('Temperature (K)','FontSize',12)
% set(gca,'FontSize',10)
% grid on
% set(gca,'GridAlpha',0.3)
% 
% %% Ex#2 part c)
% for R=[0,1,2.5,5,10]
%     B2=((pi*5/12)/(s^2+(pi*1/60)^2))+(25/s)-(300/s);
%     A2=R+(1/(s*C));
%     x_s3=A2\B2;
% 
%     t_final=1000;
%     dt=0.01;
%     t=(0:dt:t_final).';
%     x_t3=impulse(x_s3,t);
%     
%     figure(3)
%     plot(t,x_t3,'LineWidth',1.1)
%     title('Heat Transfer of Iron Block with Varying R','FontSize',14)
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Heat Transfer rate (J/s)','FontSize',12)
%     set(gca,'FontSize',10)
%     grid on
%     set(gca,'GridAlpha',0.3)
%     legend({'R=0','R=1','R=2.5','R=5','R=10'},'Location','southeast')
%     hold on
% end 
% 
% for R=[0,1,2.5,5,10]
%     B2=((pi*5/12)/(s^2+(pi*1/60)^2))+(25/s)-(300/s);
%     A2=R+(1/(s*C));
%     x_s3=A2\B2;
% 
%     t_final=2500;
%     dt=0.01;
%     t=(0:dt:t_final).';
%     
%     T=((pi*5/12)/(s^2+(pi*1/60)^2))+(25/s)-R*x_s3+(273.15/s);
% 
%     figure(4)
%     T_t3=impulse(T,t);
%     plot(t,T_t3,'LineWidth',1.1)
%     title('Temperature of Iron Block','FontSize',14)
%     xlabel('Time (s)','FontSize',12)
%     ylabel('Temperature (K)','FontSize',12)
%     set(gca,'FontSize',10)
%     grid on
%     set(gca,'GridAlpha',0.3)
%     legend({'R=0','R=1','R=2.5','R=5','R=10'},'Location','northeast')
%     hold on
% end 
% 




