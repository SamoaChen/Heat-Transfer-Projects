clc; clear all; close all;
s=tf('s');
n=30;
%%PARAMETERS
T0=25+273.15;
Tinf=25+273.15;
Tin_s=(650/s)+(300*s/(s^2+64*pi^2)); %internal temp in laplace domain
k=14.4; %conductivity
T=0.01; %m
dT=T/n;
W=0.06*2*pi; %m
H=0.15; %m
c=500;
rho=7700;%density
C=c*dT*W*H*rho; %capacity
R=dT/(k*W*H*2); %resistance

%%IMPEDENCE MATRIX

%first term
A1=[(1/R)+(1/(2*R))+s*C,-1/(2*R)];
for num =3:n
    A1=[A1,0];
end
A=A1;

%middle term
for elem=2:n-1
    Amid=[];
    for zeroL=1:elem-2
        Amid=[Amid,0];
    end
    Amid=[Amid,-1/(2*R),(1/R)+s*C,-1/(2*R)];
    for zeroR=elem:n-2
        Amid=[Amid,0];
    end
    A=cat(1,A,Amid);
end

%end term
An=[];
for num=1:n-2
    An=[An,0];
end
An=[An,-1/(2*R),(1/(2*R))+(1/R)+s*C];

%combine terms
A=cat(1,A,An);

%% FORCING MATRIX
B=[((Tinf/s)/R)+(T0/s)*s*C];
for step=2:n-1
    B=[B,(T0/s)*s*C];
end
B=[B,((Tin_s)/R)+(T0/s)*s*C].';

%% SOL FOR MULTIPLE NODES
Twall_s=A\B;

%% PLOT TEMP RESPONSE
t_final=5;
dt=0.01;
t=(0:dt:t_final).';
figure(1)
Twall_t=impulse(Twall_s,t);
% plot(t,Twall_t)
%impulse(Tin_s)
%plot for sub slabs
% twall=Twall_t(1000,1)
%create x array
x_axis=[dT/2];
for inter=1:1:n-1
    x_axis=[x_axis,inter*dT+dT/2];
end
y1_axis=[];
for inter=1:1:n
    y1_axis=[y1_axis,Twall_t(100,inter)];
end
y2_axis=[];
for inter=1:1:n
    y2_axis=[y2_axis,Twall_t(100,inter)];
end
y3_axis=[];
for inter=1:1:n
    y3_axis=[y3_axis,Twall_t(300,inter)];
end
y4_axis=[];
for inter=1:1:n
    y4_axis=[y4_axis,Twall_t(400,inter)];
end
y5_axis=[];
for inter=1:1:n
    y5_axis=[y5_axis,Twall_t(500,inter)];
end
plot(x_axis,y1_axis,'LineWidth',1.1)
hold on
plot(x_axis,y2_axis,'LineWidth',1.1)
hold on
plot(x_axis,y3_axis,'LineWidth',1.1)
hold on
plot(x_axis,y4_axis,'LineWidth',1.1)
hold on
plot(x_axis,y5_axis,'LineWidth',1.1)
title('Temperature Responses Within Cylinder wall of Thickness 0.001m','FontSize',14)
xlabel('Distance (m)','FontSize',12)
ylabel('Temperature (Kevlin)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'T at 1s','T at 2s','T at 3s','T at 4s','T at 5s'},'Location','Northwest')