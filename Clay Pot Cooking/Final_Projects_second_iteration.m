clear all;
close all;
s=tf('s');
step=20; %time step for execusion
laps=100; %total time=laps*step
%% PROPERTIES USED FOR ACCOUNTING EVAPORATING BEHAVIOR
heat_eva=2256.4e3; %J/kg
c_water_boil=heat_eva/3; %J/kg.K
% %% THERMAL PROPERTIES FOR DIFFERENT MATERIALS BEFORE BOILING
% k_clay=0.15; %W/m/K
% k_rice=(7442910468903242361*T0_rice)/4611686018427387904000 - (1090438763915875231147*T0_rice^2)/188894659314785808547840000 + 1804227/3906250; %W/m/K
% k_pork=(8083797673922461199*T0_meat)/4611686018427387904000 - (117603812228207589099*T0_meat^2)/18889465931478580854784000 + 26458201/50000000; %W/m/K
% k_carroT0_veg=(8083797673922461199*T0_veg)/4611686018427387904000 - (117603812228207589099*T0_veg^2)/18889465931478580854784000 + 26458201/50000000; %W/m/K
% rou_clay=1089; %kg/m^3
% rou_carrot=57534959/50000 - (1955971*T0_veg)/25000000; %kg/m^3
% rou_rice=121458761/100000 - (27927117*T0_rice)/250000000; %kg/m^3
% rou_pork=57534959/50000 - (1955971*T0_meat)/25000000; %kg/m^3
% c_clay=1381; %J/kg/K
% c_pork=(1536343*T0_meat)/5000000 - (403044410064435971*T0_meat^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
% c_carrot=(1536343*T0_veg)/5000000 - (403044410064435971*T0_veg^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
% c_rice=(1928763*T0_rice)/3125000 - (20298878805552515853*T0_rice^2)/11529215046068469760000 + 5916197053716758783/1717986918400000; %J/kg/K
h_air=6.8; %W/m^2/K

% %% HEAT CAPACITANCE FOR DIFFERENT MATERIALS AFTER BOILING POINT
% c_rice_boil =(1928763*T0_rice)/3125000 - (20298878805552515853*T0_rice^2)/11529215046068469760000 + 3492382347590893943/6710886400000; %capacitance of rice after boiling
%  
% c_carrot_boil =(1536343*T0_veg)/5000000 - (403044410064435971*T0_veg^2)/576460752303423488000 + 1739251903851659389/2684354560000; %capacitance of carrot after boiling
%  
% c_meat_boil =(1536343*T0_meat)/5000000 - (403044410064435971*T0_meat^2)/576460752303423488000 + 1739251903851659389/2684354560000; %capacitance of pork after boiling
 

%% DRIVING PARAMETERS (INITIAL TEMPERATURE AND CURRENT SOURCE)
q=800; %W The heat transfer from the stove
T0_bottom=25;
T0_cap=25;%initial temperature
T0_capwall=25;
T0_rice=25;
T0_veg=25;
T0_meat=25;
T0_bwall1=25;
T0_bwall2=25;
T0_bwall3=25;
T_inf=25; %initial environment temperature
T_wall=50; %temperature of the wall due to convection
%% INSTRUMENT PARAMETERS
R_pot=0.1; %m radius of pot
t_pot=0.005; %m thickness of pot
hc_pot=0.05; %m height of cap 
hb_pot=0.1; %m height of base wall
V_cap=pi*R_pot^2*t_pot; %m^3 volume of cap plate
V_base=V_cap; %m^3 volume of base plate
V_twall=(R_pot^2-(R_pot-t_pot)^2)*pi*(hc_pot-t_pot); %volume of top cap wall
V_bwall=(R_pot^2-(R_pot-t_pot)^2)*pi*(hb_pot-t_pot); %volume of bottom wall
V_veg=(R_pot-t_pot)^2*pi*(hb_pot-t_pot)/3; %volume of vegetable
V_meat=V_veg; %volume of meat
V_rice=V_veg; %volume of rice
A_in=(R_pot-t_pot)^2*pi;%inner area of cap plate

% %% CALCULATING CONDUCTIVE RESISTANCES
% Rk_vegy=((hb_pot-t_pot)/3)/(k_carrot*A_in);
% Rk_meaty=((hb_pot-t_pot)/3)/(k_pork*A_in);
% Rk_ricey=((hb_pot-t_pot)/3)/(k_rice*A_in);
% Rk_poty=t_pot/(k_carrot*R_pot^2*pi);
% Rk_vegx=((R_pot-t_pot)/2)/(k_carrot*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
% Rk_meatx=((R_pot-t_pot)/2)/(k_pork*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
% Rk_ricex=((R_pot-t_pot)/2)/(k_rice*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
% Rk_potxt=(t_pot/2)/(k_clay*(R_pot-t_pot/2)*2*pi*(hc_pot-t_pot));
% Rk_potxb=(t_pot/2)/(k_clay*(R_pot-t_pot/2)*2*pi*((hb_pot-t_pot)/3));
% 
%% CALCULATING CONVECTIVE RESISTANCE
Rh_capy=1/(h_air*A_in);
Rh_capx=1/(h_air*(R_pot-t_pot)*2*pi*(hc_pot-t_pot));
Rh_up=1/(h_air*R_pot^2*pi);
Rh_cap=1/(h_air*R_pot*2*pi*hc_pot);
Rh_base=1/(h_air*((R_pot*2*pi)*((hb_pot-t_pot)/3)));
% 
% %% CALCULATING CAPACITANCE
% C_top=V_cap*rou_clay*c_clay;
% C_cap=V_twall*rou_clay*c_clay;
% C_veg=V_veg*rou_carrot*c_carrot;
% C_meat=V_meat*rou_pork*c_pork;
% C_rice=V_rice*rou_rice*c_rice;
% C_wall=V_bwall*rou_clay*c_clay/3;

% %% IMPEDENCE MATRIX
% s=tf('s');
% 
% A=[ (2/(Rk_poty+Rk_ricey)+s*C_top)  (-2/(Rk_poty+Rk_ricey)) 0 0 0 0 0 0 0 0;
% 
%   (-2/(Rk_poty+Rk_ricey))  (2/(Rk_poty+Rk_ricey)+2/(Rk_ricey+Rk_meaty)+2/(Rk_ricex+Rk_potxb)+s*C_rice)  (-2/(Rk_ricey+Rk_meaty))  0 0 0 0  (-2/(Rk_ricex+Rk_potxb))  0 0;
% 
%   0  (-2/(Rk_ricey+Rk_meaty))  (2/(Rk_ricey+Rk_meaty)+2/(Rk_meaty+Rk_vegy)+2/(Rk_meatx+Rk_potxb)+s*C_meat)  (-2/(Rk_meaty+Rk_vegy))  0 0 0 0  (-2/(Rk_meatx+Rk_potxb))  0;
% 
%   0 0  (-2/(Rk_meaty+Rk_vegy))  (2/(Rk_meaty+Rk_vegy)+2/(Rk_vegy+2*Rh_capy)+2/(Rk_vegx+Rk_potxb)+s*C_veg)  (-2/(Rk_vegy+2*Rh_capy))  0 0 0 0  (-2/(Rk_vegx+Rk_potxb));
% 
%   0 0 0  (-2/(Rk_vegy+2*Rh_capy))  (2/(Rk_vegy+2*Rh_capy)+2/(Rk_poty)+2/(Rk_potxt+2*Rh_capx))  (-2/(Rk_poty))  (-2/(Rk_potxt+2*Rh_capx))  0 0 0;
% 
%   0 0 0 0  (-2/(Rk_poty))  (2/(Rk_poty+2*Rh_up)+2/(Rk_poty)+s*C_top) 0 0 0 0;
% 
%   0  (-2/(2*Rk_ricex+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_ricex+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall)  0 0;
% 
%   0 0  (-2/(2*Rk_meatx+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_meatx+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall)  0;
% 
%   0 0 0  (-2/(2*Rk_vegx+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_vegx+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall);
% 
%   0 0 0 0  (-2/(Rk_poty))  0  (2/(2*Rh_capx+Rk_potxt)+2/(2*Rh_cap+Rk_potxt)+s*C_cap)  0 0 0];
% 
% 
% F=[ q+s*C_top*(T0_bottom/s);
%   s*C_rice*(T0_rice/s);
%   s*C_meat*(T0_meat/s);
%   s*C_veg*(T0_veg/s);
%   0;
%   2/(Rk_poty+2*Rh_up)*(T_inf/s)+s*C_top*(T0_cap/s);
%   2/(Rk_potxb+2*Rh_base)*(T_inf/s)+s*C_wall*(T0_bwall1/s);
%   2/(Rk_potxb+2*Rh_base)*(T_inf/s)+s*C_wall*(T0_bwall2/s);
%   2/(Rk_potxb+2*Rh_base)*(T_inf/s)+s*C_wall*(T0_bwall3/s);
%   2/(Rk_potxt+2*Rh_up)*(T_inf/s)+s*C_cap*(T0_capwall/s) ];



%% ITERATIVE PROCESS OF PLOTTING THE TEMPERATURE RESPONSE
for n=1:1:laps
    time=[0:step/10:step].';
    %%CALCULATING FOOD COMPONENTS PROPERTIES
    k_clay=0.15; %W/m/K
    k_rice=(7442910468903242361*T0_rice)/4611686018427387904000 - (1090438763915875231147*T0_rice^2)/188894659314785808547840000 + 1804227/3906250; %W/m/K
    k_pork=(8083797673922461199*T0_meat)/4611686018427387904000 - (117603812228207589099*T0_meat^2)/18889465931478580854784000 + 26458201/50000000; %W/m/K
    k_carrot=(8083797673922461199*T0_veg)/4611686018427387904000 - (117603812228207589099*T0_veg^2)/18889465931478580854784000 + 26458201/50000000; %W/m/K
    rou_clay=1089; %kg/m^3
    rou_carrot=57534959/50000 - (1955971*T0_veg)/25000000; %kg/m^3
    rou_rice=121458761/100000 - (27927117*T0_rice)/250000000; %kg/m^3
    rou_pork=57534959/50000 - (1955971*T0_meat)/25000000; %kg/m^3
    c_clay=1381; %J/kg/K
%     c_pork=(1536343*T0_meat)/5000000 - (403044410064435971*T0_meat^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
%     c_carrot=(1536343*T0_veg)/5000000 - (403044410064435971*T0_veg^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
%     c_rice=(1928763*T0_rice)/3125000 - (20298878805552515853*T0_rice^2)/11529215046068469760000 + 5916197053716758783/1717986918400000; %J/kg/K
    h_air=6.8; %W/m^2/K
    %%CALCULATING THE CONDUCTIVE RESISTANCE AGAIN
    Rk_vegy=((hb_pot-t_pot)/3)/(k_carrot*A_in);
    Rk_meaty=((hb_pot-t_pot)/3)/(k_pork*A_in);
    Rk_ricey=((hb_pot-t_pot)/3)/(k_rice*A_in);
    Rk_poty=t_pot/(k_carrot*R_pot^2*pi);
    Rk_vegx=((R_pot-t_pot)/2)/(k_carrot*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
    Rk_meatx=((R_pot-t_pot)/2)/(k_pork*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
    Rk_ricex=((R_pot-t_pot)/2)/(k_rice*((R_pot-t_pot)/2)*2*pi*((hb_pot-t_pot)/3));
    Rk_potxt=(t_pot/2)/(k_clay*(R_pot-t_pot/2)*2*pi*(hc_pot-t_pot));
    Rk_potxb=(t_pot/2)/(k_clay*(R_pot-t_pot/2)*2*pi*((hb_pot-t_pot)/3));
    if T0_rice<100
        %%CALCULATING CAPACITANCE AGAIN
        c_rice=(1928763*T0_rice)/3125000 - (20298878805552515853*T0_rice^2)/11529215046068469760000 + 5916197053716758783/1717986918400000; %J/kg/K
        C_top=V_cap*rou_clay*c_clay;
        C_cap=V_twall*rou_clay*c_clay;
        C_rice=V_rice*rou_rice*c_rice;
        C_wall=V_bwall*rou_clay*c_clay/3;
    else
        %%HEAT CAPACITANCE FOR DIFFERENT MATERIALS AFTER BOILING POINT 
        c_rice_boil =(1928763*T0_rice)/3125000 - (20298878805552515853*T0_rice^2)/11529215046068469760000 + 3492382347590893943/6710886400000; %capacitance of rice after boiling
        C_top=V_cap*rou_clay*c_clay;
        C_cap=V_twall*rou_clay*c_clay;
        C_rice=V_rice*rou_rice*c_rice_boil;
        C_wall=V_bwall*rou_clay*c_clay/3;
    end
    if T0_meat<100
        %%CALCULATING CAPACITANCE AGAIN
        c_pork=(1536343*T0_meat)/5000000 - (403044410064435971*T0_meat^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
        C_meat=V_meat*rou_pork*c_pork;
    else
        %%HEAT CAPACITANCE FOR DIFFERENT MATERIALS AFTER BOILING POINT 
        c_meat_boil =(1536343*T0_meat)/5000000 - (403044410064435971*T0_meat^2)/576460752303423488000 + 1739251903851659389/2684354560000; %capacitance of pork after boiling
        C_meat=V_meat*rou_pork*c_meat_boil;
    end
    if T0_veg<100
        %%CALCULATING CAPACITANCE AGAIN
        c_carrot=(1536343*T0_veg)/5000000 - (403044410064435971*T0_veg^2)/576460752303423488000 + 2723544178070912959/687194767360000; %J/kg/K
        C_veg=V_veg*rou_carrot*c_carrot;
    else
        %%HEAT CAPACITANCE FOR DIFFERENT MATERIALS AFTER BOILING POINT 
        c_carrot_boil =(1536343*T0_veg)/5000000 - (403044410064435971*T0_veg^2)/576460752303423488000 + 1739251903851659389/2684354560000; %capacitance of carrot after boiling
        C_veg=V_veg*rou_carrot*c_carrot_boil;
    end
    %%CALCULATING THE IMPEDENCE AND FORCING MATRIX AGAIN
    A=[ (2/(Rk_poty+Rk_ricey)+s*C_top)  (-2/(Rk_poty+Rk_ricey)) 0 0 0 0 0 0 0 0;

      (-2/(Rk_poty+Rk_ricey))  (2/(Rk_poty+Rk_ricey)+2/(Rk_ricey+Rk_meaty)+2/(Rk_ricex+Rk_potxb)+s*C_rice)  (-2/(Rk_ricey+Rk_meaty))  0 0 0 0  (-2/(Rk_ricex+Rk_potxb))  0 0;

      0  (-2/(Rk_ricey+Rk_meaty))  (2/(Rk_ricey+Rk_meaty)+2/(Rk_meaty+Rk_vegy)+2/(Rk_meatx+Rk_potxb)+s*C_meat)  (-2/(Rk_meaty+Rk_vegy))  0 0 0 0  (-2/(Rk_meatx+Rk_potxb))  0;

      0 0  (-2/(Rk_meaty+Rk_vegy))  (2/(Rk_meaty+Rk_vegy)+2/(Rk_vegy+2*Rh_capy)+2/(Rk_vegx+Rk_potxb)+s*C_veg)  (-2/(Rk_vegy+2*Rh_capy))  0 0 0 0  (-2/(Rk_vegx+Rk_potxb));

      0 0 0  (-2/(Rk_vegy+2*Rh_capy))  (2/(Rk_vegy+2*Rh_capy)+2/(Rk_poty)+2/(Rk_potxt+2*Rh_capx))  (-2/(Rk_poty))  (-2/(Rk_potxt+2*Rh_capx))  0 0 0;

      0 0 0 0  (-2/(Rk_poty))  (2/(Rk_poty+2*Rh_up)+2/(Rk_poty)+s*C_top) 0 0 0 0;

      0  (-2/(2*Rk_ricex+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_ricex+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall)  0 0;

      0 0  (-2/(2*Rk_meatx+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_meatx+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall)  0;

      0 0 0  (-2/(2*Rk_vegx+Rk_potxb))  0 0 0 0 0  (2/(2*Rk_vegx+Rk_potxb)+2/(2*Rh_base+Rk_potxb)+s*C_wall);

      0 0 0 0  (-2/(Rk_poty))  0  (2/(2*Rh_capx+Rk_potxt)+2/(2*Rh_cap+Rk_potxt)+s*C_cap)  0 0 0];
    F=[ q/s+s*C_top*(T0_bottom/s);
      s*C_rice*(T0_rice/s);
      s*C_meat*(T0_meat/s);
      s*C_veg*(T0_veg/s);
      0;
      2/(Rk_poty+2*Rh_up)*(T_inf/s)+s*C_top*(T0_cap/s);
      2/(Rk_potxb+2*Rh_base)*(T_wall/s)+s*C_wall*(T0_bwall1/s);
      2/(Rk_potxb+2*Rh_base)*(T_wall/s)+s*C_wall*(T0_bwall2/s);
      2/(Rk_potxb+2*Rh_base)*(T_wall/s)+s*C_wall*(T0_bwall3/s);
      2/(Rk_potxt+2*Rh_up)*(T_wall/s)+s*C_cap*(T0_capwall/s) ];
  
    X=A\F; %X is composed of temp matrix Tb;T1;T2;T3;T4;Ttop;Twtop;Twb1;Twb2;Twb3;
    
    %%TEMPERATURE IN S DOMAIN
    Tb_s=X(1);
    T1_s=X(2);
    T2_s=X(3);
    T3_s=X(4);
    T4_s=X(5);
    Ttop_s=X(6);
    Twtop_s=X(7);
    Twb1_s=X(8);
    Twb2_s=X(9);
    Twb3_s=X(10);
    %%IN t DOMAIN
    Tb_t=impulse(Tb_s,time);
    T1_t=impulse(T1_s,time);
    T2_t=impulse(T2_s,time);
    T3_t=impulse(T3_s,time);
    T4_t=impulse(T4_s,time);
    Ttop_t=impulse(Ttop_s,time);
    Twtop_t=impulse(Twtop_s,time);
    Twb1_t=impulse(Twb1_s,time);
    Twb2_t=impulse(Twb2_s,time);
    Twb3_t=impulse(Twb3_s,time);
    %%PLOT
    figure(1)
%     plot(time+step*n,Tb_t,'r','LineWidth',1.1)
%     hold on
    plot(time+step*n,T1_t,'y','LineWidth',1.1)
    hold on
    plot(time+step*n,T2_t,'--','LineWidth',1.1)
    hold on
    plot(time+step*n,T3_t,'--','LineWidth',1.1)
    hold on
    plot(time+step*n,T4_t,'--','LineWidth',1.1)
    hold on
    plot(time+step*n,Ttop_t,'g','LineWidth',1.1)
    hold on
    plot(time+step*n,Twtop_t,'b','LineWidth',1.1)
    hold on
    plot(time+step*n,Twb1_t,'k','LineWidth',1.1)
    hold on
    plot(time+step*n,Twb2_t,'m','LineWidth',1.1)
    hold on
    plot(time+step*n,Twb3_t,'c','LineWidth',1.1)
    hold on
    %%CHANGING DRIVING TEMPERATURES
    T0_bottom=Tb_t(11);%initial temperature
    T0_rice=T1_t(11);
    T0_meat=T2_t(11);
    T0_veg=T3_t(11);
    T0_cap=Ttop_t(11);
    T0_capwall=Twtop_t(11);
    T0_bwall1=Twb1_t(11);
    T0_bwall2=Twb2_t(11);
    T0_bwall3=Twb3_t(11);
   
end

% %% TEMPERATURE RESPONSES 
% %%IN S DOMAIN
% Tb_s=X(1);
% T1_s=X(2);
% T2_s=X(3);
% T3_s=X(4);
% T4_s=X(5);
% Ttop_s=X(6);
% Twtop_s=X(7);
% Twb1_s=X(8);
% Twb2_s=X(9);
% Twb3_s=X(10);
% 
% 
% %%IN t DOMAIN
% % time=(0:0.1:100).';
% Tb_t=impulse(Tb_s,time);
% T1_t=impulse(T1_s,time);
% T2_t=impulse(T2_s,time);
% T3_t=impulse(T3_s,time);
% T4_t=impulse(T4_s,time);
% Ttop_t=impulse(Ttop_s,time);
% Twtop_t=impulse(Twtop_s,time);
% Twb1_t=impulse(Twb1_s,time);
% Twb2_t=impulse(Twb2_s,time);
% Twb3_t=impulse(Twb3_s,time);
% 
% figure(1)
% plot(t,Tb_t,'LineWidth',1.1)
% hold on
% plot(time,T1_t,'LineWidth',1.1)
% hold on
% plot(time,T2_t,'LineWidth',1.1)
% hold on
% plot(time,T3_t,'LineWidth',1.1)
% hold on
% plot(time,T4_t,'LineWidth',1.1)
% hold on
% plot(time,Ttop_t,'LineWidth',1.1)
% hold on
% plot(time,Twtop_t,'LineWidth',1.1)
% hold on
% plot(time,Twb1_t,'LineWidth',1.1)
% hold on
% plot(time,Twb2_t,'LineWidth',1.1)
% hold on
% plot(time,Twb3_t,'LineWidth',1.1)
title('Temperature Responses of the Base','FontSize',14)
xlabel('Time (s)','FontSize',12)
ylabel('Temperature (Degree Celsius)','FontSize',12)
set(gca,'FontSize',10)
grid on
set(gca,'GridAlpha',0.3)
legend({'T-rice','T-meat','T-veg','T4','Ttop','Twtop','Twb1','Twb2','Twb3'},'Location','bestoutside')



