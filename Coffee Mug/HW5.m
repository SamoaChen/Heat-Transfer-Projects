clc;
clear all;
close all;
%% QUESTION 2
t_final=1e6;
dt=4;
t=(0:dt:t_final).';
TH=400+273.15; %c
TL=100+273.15; %c 
T0=400+273.15; % c
k_w=0.16; %W/m*K
Rw=0.2/(k_w*0.3*0.5);%K/w
e=0.89;
sig=5.669E-8; %w/m^2*K^4
hrH=sig*e*4*TH^3;
hrL=sig*e*4*TL^3;
RrH=1/(hrH*0.3*0.5); %K/w
RrL=1/(hrL*0.3*0.5); %K/w

hair=6.8; %w/m^2*k
Rh=1/(hair*0.3*0.5); %K/w
C=2380*0.2*0.3*0.5*740; %J/K
%part a)
s=tf('s');
A=[Rh+Rw/2+1/(s*C),-1/(s*C);
    -1/(s*C),1/(s*C)+Rw/2+Rh];
B=[(TH-T0)/s,(T0-TL)/s].';
q_s=A\B;
T_left=TH/s-q_s(1)*Rh;
T_left_t=impulse(T_left,t)-273.15;
T_wall=TH/s-q_s(1)*(Rh+Rw/2);
T_wall_t=impulse(T_wall,t)-273.15;
T_right=TL/s+q_s(2)*Rh;
T_right_t=impulse(T_right,t)-273.15;

figure(1)
plot(t,T_left_t,'LineWidth',1.1)
hold on
plot(t,T_wall_t,'LineWidth',1.1)
hold on
plot(t,T_right_t,'LineWidth',1.1)
title('Temperature of the Wall and Its Surfaces','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Temperature (degree celsius)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'T LeftWall','T Wall','T RightWall'},'Location','bestoutside')

%% part b
RhH=(RrH*Rh)/(RrH+Rh);
RhL=(RrL*Rh)/(RrL*Rh);

A=[RhH+Rw/2+1/(s*C),-1/(s*C);
    -1/(s*C),1/(s*C)+Rw/2+RhL];
B=[(TH-T0)/s,(T0-TL)/s].';
q_s=A\B;
T_left=TH/s-q_s(1)*RhH;
T_left_t=impulse(T_left,t)-273.15;
T_wall=TH/s-q_s(1)*(RhH+Rw/2);
T_wall_t=impulse(T_wall,t)-273.15;
T_right=TL/s+q_s(2)*RhL;
T_right_t=impulse(T_right,t)-273.15;

figure(2)
plot(t,T_left_t,'LineWidth',1.1)
hold on
plot(t,T_wall_t,'LineWidth',1.1)
hold on
plot(t,T_wall_t,'LineWidth',1.1)
hold on
plot(t,T_right_t,'LineWidth',1.1)
title('Temperature of the Wall and Its Surfaces with Radiation','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Temperature (degree celsius)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'T LeftWall','T Wall','T RightWall'},'Location','bestoutside')


%% part c
T_left_tem=T_left_t(250000)+273.15;
T_left_sur=400+273.15;
T_right_tem=T_right_t(250000)+273.15;
T_right_sur=100+273.15
while (T_left_sur-T_left_tem) >0.001
    T_left_sur=T_left_sur-0.6*(T_left_sur-T_left_tem);
    while (T_right_tem-T_right_sur) >0.001;
        T_right_sur=T_right_sur+0.6*(T_right_tem-T_right_sur);
        hrH=sig*e*4*T_left_sur^3;
        hrL=sig*e*4*T_right_sur^3;
        RrH=1/(hrH*0.3*0.5)
        RrL=1/(hrL*0.3*0.5)
        
        RhH=(RrH*Rh)/(RrH+Rh);
        RhL=(RrL*Rh)/(RrL*Rh);

        A=[RhH+Rw/2+1/(s*C),-1/(s*C);
            -1/(s*C),1/(s*C)+Rw/2+RhL];
        B=[(TH-T0)/s,(T0-TL)/s].';
        q_s=A\B;
        T_left=TH/s-q_s(1)*RhH;
        T_right=TL/s+q_s(2)*RhL;
        T_left_t=impulse(T_left,t);
        T_left_tem=T_left_t(250000);
        T_right_t=impulse(T_right,t);
        T_right_tem=T_right_t(250000);
    end
end

hrH=sig*e*4*T_left_sur^3;
hrL=sig*e*4*T_right_sur^3;
RrH=1/(hrH*0.3*0.5); 
RrL=1/(hrL*0.3*0.5);  

A=[RhH+Rw/2+1/(s*C),-1/(s*C);
    -1/(s*C),1/(s*C)+Rw/2+RhL];
B=[(TH-T0)/s,(T0-TL)/s].';
q_s=A\B;
T_left=TH/s-q_s(1)*RhH;
T_left_t=impulse(T_left,t)-273.15;
T_wall=TH/s-q_s(1)*(RhH+Rw/2);
T_wall_t=impulse(T_wall,t)-273.15;
T_right=TL/s+q_s(2)*RhL;
T_right_t=impulse(T_right,t)-273.15;

figure(3)
plot(t,T_left_t,'LineWidth',1.1)
hold on
plot(t,T_wall_t,'LineWidth',1.1)
hold on
plot(t,T_wall_t,'LineWidth',1.1)
hold on
plot(t,T_right_t,'LineWidth',1.1)
title('Temperature of the Wall and Its Surfaces Optimized','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Temperature (degree celsius)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'T LeftWall','T Wall','T RightWall'},'Location','bestoutside')


