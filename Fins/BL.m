clc;
clear all;
close all;

%% PARAMETERS
k_air=0.0257; %W/m/k
k_water=0.6; %W/m/k
rou_air=1.205; %kg/m^3
rou_water=1000; %kg/m^3
c_p_air=1009; %J/kg/K
c_p_water=4200; %J/kg/k
mu_air=1.82e-5; %kg/s/m
mu_water=1e-3; %kg/s/m
beta_air=3.41e-3; %1/k
beta_water=2.07e-4; %1/k
D_f=0.02; %m
D_a=0.08; %m
T_s=32; %degree celsius
%% FORCED CONVECTION COEFFICENT
hf_falw=cal_h(rou_air,0.447,D_f,mu_air,c_p_air,k_air) %finger at light wind

hf_fahw=cal_h(rou_air,4.47,D_f,mu_air,c_p_air,k_air) %finger at a windy day

ha_falw=cal_h(rou_air,4.47,D_a,mu_air,c_p_air,k_air) %forearm at a windy day

hf_fwrw=cal_h(rou_water,1,D_f,mu_water,c_p_water,k_water) %finger at running water

ha_fwrw=cal_h(rou_water,1,D_a,mu_water,c_p_water,k_water) %forearm at running water

%% NATURAL CONVECTION
hf_nac=cal_nh(beta_air,rou_air,0,T_s,D_f,mu_air,c_p_air,k_air) %finger at T_inf=0

hf_narc=cal_nh(beta_air,rou_air,-15,T_s,D_f,mu_air,c_p_air,k_air) %finger at T_inf=-15

ha_nac=cal_nh(beta_air,rou_air,0,T_s,D_a,mu_air,c_p_air,k_air) %forearm at 0

hf_nwc=cal_nh(beta_water,rou_water,0,T_s,D_f,mu_water,c_p_water,k_water) %finger in water at 0

ha_nwc=cal_nh(beta_water,rou_water,0,T_s,D_a,mu_water,c_p_water,k_water) %forearm in water at 0


%% FUNCTIONS
function [h]=cal_h(rou,V,D,mu,c_p,k)
    C=[0.989,0.911,0.683,0.193,0.0266];
    Re_list=[0.4,4,35,4083,40045,400000];
    n=[0.33,0.385,0.446,0.618,0.805];
    Re=rou*V*D/mu;
    Re_loc=sort([Re_list,Re]);
    for i=1:1:6
        if Re_loc(i)==Re
            num=i-1;
        end
    end
    Pr=(mu/rou)/(k/(rou*c_p));
    Nud=C(num)*(Re^n(num))*Pr^(1/3);
    h=Nud*k/D;
end

function [h]=cal_nh(beta,rou,T_inf,T_s,D,mu,c_p,k)
    C=[0.53,0.13];
    N=[0.25,0.3333];
    Ray_list=[1e4,2.12e7];
    Pr=(mu/rou)/(k/(rou*c_p));
    Ray=Pr*9.8*beta*(D^3)*sqrt((T_s-T_inf)^2)/((mu/rou)^2);
    Ray_loc=sort([Ray_list,Ray]);
    for i=1:1:3
        if Ray_loc(i)==Ray
            num=i-1;
        end
    end
    Nud=C(num)*Ray^N(num);
    h=Nud*k/D;
end
