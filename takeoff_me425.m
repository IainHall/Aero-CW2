clear all;
%% Known Variables
Bmat = 1:0.1:38;
matv = zeros(length(Bmat),1);
y = 1.4;
yg = 1.3;
Ec = 0.9;
Cp = 1005;
Cpg = 1244;
Ma_crit = 1;
R = 287;
U = 231.3;

%Design
T_04 = 1500;         %K
r_fan = 1.5;         %-, fan pressure ratio = P3/P1
r_boo = 1.667;      %-, booster pressure ratio = P23/P2
r_hpc = 18;          %-, HPC pressure ratio = P3/P2
F_T = 38900;        %N, thrust of each engine

%Cruise
M_2 = 0.78;        %-, Mach at cruise
H_cr = 35000;       %ft, cruise altitude, values for density, pressure and temperature have been calculated separately
rho_2 = 0.3794;     %kg/m3, density of free stream
P_2 = 23.8418;      %kPa, pressure of free stream
T_2 = 218.9076;     %K, temperature of free stream, calculated using standard day assumption


%% Simple shit
%             ----------------------------------------
%             ||(13)
%             ||------------------------------------
%             ||      |  \     |     \
%         (2) ||(13)  |   |(23)|      |(3)  -----[    ]
%             ||      |  /     |     /
%             ||------------------------------------
%             ||
%             -----------------------------------------
V_2 = M_2*sqrt(y*R*T_2);
T_02 = T_2*(1+((y-1)/2)*M_2^2);         %K, stag temp before fan
P_02 = P_2*(T_02/T_2)^(y/(y-1));        %kPa, stag pressure before fan

P_013 = r_fan*P_02;                     %kPa, stag pressure after fan
t_fan = r_fan^((y-1)/y);
T_013 = t_fan*T_02;

P_023 = r_boo*P_013;
t_boo = r_boo^((y-1)/(y*Ec));
T_023 = t_boo*T_013;

P_03 = r_hpc*P_023;
t_hpc = r_hpc.^((y-1)/(y*Ec));
T_03 = t_hpc*T_023;

P_04 = P_03;                            %kPa, stag pressure before HPC, no temperature drop across comb chamber
%T_04 is known

T_045 = T_04 - (Cp/Cpg)*(T_03-T_023);
% t_hpt = (T_04-T_045)/T_04;
t_hpt = T_045/T_04;
r_hpt = t_hpt^((yg)/(Ec*(yg-1)));
P_045 = r_hpt*P_04;

%19
ratioC=(1+((y-1)./2)).^(y./(y-1));    %P_013/P*
ratioD=P_013./P_2;                    %P_013/P2

if ratioD > ratioC
    disp('choked at 19');
else
    disp('not choked at 19')
end 


P_19 = P_013.*(1./(ratioC));
%P_19 = P_019/((1+0.5*(y-1)*Ma_crit^2)^((y)/(y-1)));
T_19 = T_013.*(P_19./P_013).^((y-1)./(y));
V_19 = (y.*R.*T_19).^(0.5);
rho_19 = P_19/(0.287*T_19);

foundB=0;
% for matnum = 1:length(matv)
%     B = Bmat(matnum)
% 
%     T_05 = T_045-(Cp/Cpg)*((1+B)*(T_013-T_02)+(T_023-T_013));
% %     T_05 = T_04 - (Cp./(Cpg)).*(((1+B).*(T_013-T_02))+(T_023-T_013)+(T_03-T_023));
%     P_05 = P_045*(T_05/T_045)^((yg)/(Ec*(yg-1)));
% 
%     
%  
% 
%     ratioA=(P_05./P_2);                             %P_05/P2
%     ratioB=(1+((y-1)./2).*Ma_crit).^(y./(y-1));     %P_05/P*
% 
%     if ratioA > ratioB
%     disp('choked flow')
%     else 
%     disp('not choked flow')
%     end 
% 
%     P_9=P_05.*(1./ratioB);
%     T_9=T_05.*(P_9./P_05).^((yg-1)./yg);
%     V_9=(yg.*R.*T_9).^0.5;
%     rho_9 = P_9/(R*T_9);
%     matv(matnum,1) = V_9-V_19;
% %     if abs(V_9-V_19)<5 && foundB==0
% %         matv(matnum,1) = V_9-V_19;
% %         foundB=1;
% %         rightB=B;
% %     end
% end 

B = 14.35

T_05 = T_045-(Cp/Cpg)*((1+B)*(T_013-T_02)+(T_023-T_013));
%     T_05 = T_04 - (Cp./(Cpg)).*(((1+B).*(T_013-T_02))+(T_023-T_013)+(T_03-T_023));
P_05 = P_045*(T_05/T_045)^((yg)/(Ec*(yg-1)));
r_lpt = P_05/P_045;

ratioA=(P_05./P_2);                             %P_05/P2
ratioB=(1+((yg-1)./2).*Ma_crit).^(yg./(yg-1));     %P_05/P*

if ratioA > ratioB
    disp('choked flow')
    P_9=P_05.*(1./ratioB);
    T_9=T_05.*(P_9./P_05).^((yg-1)./yg);
    V_9=(yg.*R.*T_9).^0.5;
    rho_9c = P_9/(0.287*T_9);
    diff = V_9-V_19
else 
    disp('not choked flow')
    P_9=P_2;
    T_9=T_05.*(P_9./P_05).^((yg-1)./yg);
    
    M_9 = ((2./(yg-1)).*((P_05./P_9).^((yg-1)./yg)-1)).^0.5;
    %T_9 = ((V_19./M_9).^2).*yg.*R;
    V_9 = sqrt(Cpg*(T_05-T_9))
     %V_9=M_9.*(yg.*R.*T_9).^0.5;
     
    rho_9 = P_9./(0.287.*T_9);
    diff = V_9-V_19
end 

%Q9 = (1+B)*(Cp/Cpg)*((T_013-T_02)/(T_045-T_05));


X = (V_9-V_2)+(P_9-P_2)/(rho_9*V_9);
Y = (V_19-V_2)+(P_19-P_2)/(rho_19*V_19);
Mc = F_T/(X+B*Y);
Mb = B*Mc;
Mtot = (1+B)*Mc;

Q9A = (Mb./Mc)*(Cp/Cpg)*((T_013-T_02)/(T_045-T_05));

%% Areas
% Core Nozzle Area
A_9 = Mc./(V_9.*rho_9);

% Bypass Nozzle Area
A_19 = Mb./(V_19.*rho_19);

% Turbine Area
T_4 = T_04./(1+((yg-1)./2)*Ma_crit);
V_4 = (yg.*R.*T_4);
P_4 = P_04./(1+((yg-1)./2)).^(yg./(yg-1));
rho_4 = P_4./(R.*T_4);
A_4 = Mc./(V_4.*rho_4);






%% Take-Off

%Sea-Level Conditions
T_2t = 288.15;
P_2t = 101.325;
rho_2t = 1.225;

V_t_horoz = 90;
V_t_vert = 3;
V_t = (((V_t_horoz).^2)+((V_t_vert).^2)).^(0.5);


Ma_2t = V_t./((y.*R.*T_2t).^0.5);

T_02t = T_2t.*(1+((y-1)./2).*(Ma_2t).^2);
P_02t = P_2t.*(1+((y-1)./2).*(Ma_2t).^2).^(y./(y-1));


%%GUESS
FPR2 = 3.35085;                         %2.7385046
pi_c = 18;

P_013t = FPR2.*P_02t;                     %kPa, stag pressure after fan
t_fan2 = FPR2.^((y-1)./y);
T_013t = t_fan2.*T_02t;

ratioE=(1+((y-1)./2)).^(y./(y-1));      %P_013/P*
ratioF=P_013t./P_2t;                    %P_013/P2

P_19t = P_013t.*(1./(ratioE));
T_19t = T_013t.*(P_19t./P_013t).^((y-1)./(y));
V_19t = (y.*R.*T_19t).^(0.5);
rho_19t = P_19t./(0.287.*T_19t);


T_023t = (T_013t-T_02t)*t_boo +T_02t;
P_023t = (P_013t).*(T_023t./T_013t).^((y.*Ec)./(y-1));

%
%T_03t = t_hpc.*T_023t;
%P_03t = P_023t.*(T_03t./T_023t).^((y.*Ec)./(y-1));

%P_03t = P_023t.*r_hpc;
%T_03t = t_hpc.*T_023t;

P_03t = pi_c * P_023t
T_03t = T_023t*(P_03t/P_023t)^((y-1)/(y*Ec))

T_04t = (T_02t*T_04/T_02);
P_4t=P_013t.*(1./ratioE);
%P_04t = P_03t;
T_4t = T_04t./(1+((yg-1)./2));
P_04t = P_4t.*(1+((yg-1)./2)).^(yg./(yg-1));
V_4t = (yg.*R.*T_4t).^(0.5);
rho_4t = P_4t./(R.*T_4t);
Mc4t = rho_4t.*A_4.*V_4t;

T_045t = T_04t.*(t_hpt+1);
P_045t = P_04t*(T_045t/T_04t)^(yg/(Ec*(y-1)));


%T_03t = (Cpg./Cp).*(T_04t-T045t) + T_023t;


T_05t = T_045t-((Cp/Cpg).*((1+B)*(T_013t-T_02t)+(T_023t-T_013t)));
P_05t = P_045t.*(T_05t./T_045t).^((yg)/(Ec*(yg-1)));
r_lpt_t = P_05t./P_045t;

ratioG=(P_05t./P_2t);                             %P_05/P2
ratioH=(1+((yg-1)./2).*Ma_crit).^(yg./(yg-1));     %P_05/P*

if ratioG > ratioH
    disp('choked flow')
    P_9t=P_05t.*(1./ratioH);
    T_9t=T_05t.*(P_9t./P_05t).^((yg-1)./yg);
    V_9t=(yg.*R.*T_9t).^0.5;
    rho_9t = P_9t/(0.287*T_9t);
    
else 
    disp('not choked flow')
    P_9t=P_2t;
    T_9t=T_05t.*(P_9t./P_05t).^((yg-1)./yg);
    
    M_9t = ((2./(yg-1)).*((P_05t./P_9t).^((yg-1)./yg)-1)).^0.5;
    %T_9 = ((V_19./M_9).^2).*yg.*R;
    V_9t = sqrt(Cpg*(T_05t-T_9t))
     %V_9=M_9.*(yg.*R.*T_9).^0.5;
     
    rho_9t = P_9t./(0.287.*T_9t);
    diff = V_9t-V_19t
end 

%Q9 = (1+B)*(Cp/Cpg)*((T_013-T_02)/(T_045-T_05));


Xt = (V_9t-V_t)+(P_9t-P_2t)./(rho_9t.*V_9t);
Yt = (V_19t-V_t)+(P_19t-P_2t)./(rho_19t.*V_19t);
%F_T = 
Mc9t = rho_9t.*A_9.*V_9t;
%Mc9t = F_T./(Xt+(B.*Yt));
Mb9t = B.*Mc9t;
Mtot_t = (1+B).*Mc9t;

diff = Mc4t - Mc9t



