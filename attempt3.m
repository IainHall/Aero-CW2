


%% attempt 3


%%

%ambient conditions
Pa = 23.842;
rhoa = 0.3799;
Ta = 218.9;
Ma = 0.78;
Ua = 231.3;

%other properties
gamma_comp = 1.4;  %gamma air
gamma_turb = 1.3;  %gamma exh
Cp_comp = 1005;
Cp_turb = 1244;
Rideal = 287;
ec = 0.9;
Fengine = 38900;

% guess
B = 13.74;

%compressor ratios
R_fan= 1.5; %fan
R_boo = 1.667; %booster
R_hpc = 18; %bypass

%inlet connies
Po_2 = Pa*(1+((gamma_comp-1)/2)*Ma^2)^(gamma_comp/(gamma_comp-1));
To_2 = Ta + (Ua^2)/(2*Cp_comp);

% fan outlet conditions
Po_13 = R_fan*Po_2;
To_13 = To_2 * R_fan^((gamma_comp-1)/(gamma_comp*ec));

% booster outlet conditions

Po_23 = R_boo*Po_13;
To_23 = To_13 * R_boo^((gamma_comp-1)/(gamma_comp*ec));


% hpc outlet conditions

Po_3 = R_hpc*Po_23;
To_3 = To_23 * R_hpc^((gamma_comp-1)/(gamma_comp*ec));


% combustion chamber outlet conditions

Po_4 = Po_3;
To_4 = 1500;

% HP turbine outlet conditions: 45
To_45 = To_4 - (Cp_comp/Cp_turb)* (To_3-To_23);
HPTtempratio = To_45/To_4;
Po_45 = Po_4*(HPTtempratio)^((gamma_turb)/(ec*(gamma_turb-1)));
R_hpt= Po_45/Po_4;

% LP turbine outlet conditions: 5

To_5 = To_45 - (Cp_comp/Cp_turb)*(To_23- To_13 + (B+1)*(To_13-To_2));
LPTtempratio = To_5/To_45;
Po_5 = Po_45*(LPTtempratio)^((gamma_turb)/(ec*(gamma_turb-1)));
R_lpt = Po_5/Po_45;



ChokedRatioCore = ((gamma_turb+1)/2)^(gamma_turb/(gamma_turb-1));
CoreexitRatio = Po_5/Pa;

if CoreexitRatio>=ChokedRatioCore

    disp('core is choked')
    % 8 refers to section 9 of the engine
    T_9 = T_9*(1/ChokedRatioCore)^((gamma_turb-1)/gamma_turb);
    U_9 = sqrt(gamma_turb*Rideal*T(8))
else
    disp('core is not choked')
    P_9 = Pa;
    T_9 = To_5*(P_9/Po_5)^((gamma_turb-1)/gamma_turb);
    %Mexitcore = (CoreexitRatio^((gamma_turb-1)/gamma_turb)-1)*(2/(gamma_turb-1));
    %Uexitcore = Mexitcore*sqrt(gamma_turb*Rideal*T(8))
    U9 = sqrt(2 * Cp_turb *  (To_5-T_9  ) )
end


ChokedRatioBypass = ((gamma_comp+1)/2)^(gamma_comp/(gamma_comp-1));
BypassexitRatio = Po_13/Pa;

if BypassexitRatio>=ChokedRatioBypass
    P_19 = Po_13/ChokedRatioBypass;
    disp('Bypass is choked')
    % 8 refers to section 9 of the engine
    T_19= To_13*(1/ChokedRatioBypass)^((gamma_comp-1)/gamma_comp);
    U19= sqrt(gamma_comp*Rideal*T_19)
else
    disp('Bypass is not choked in cruise')
end


% for quiz
Worklpt= Cp_turb*(To_45-To_5);
Worklptonbypass = Cp_comp* B*(To_13-To_2);
proportionWorklpt_onBypass = Worklptonbypass/Worklpt;


FspecCore = (1/(1+B))*(U9-Ua);
FspecBypass = (B/(1+B))*(U19-Ua)+ (B/(1+B))* ((Rideal*T_19)/(1000*P_19*U19) ) * (1000*P_19- 1000*Pa);
Fspec = FspecBypass+FspecCore;
mtotal = Fengine/Fspec;
mbypass = (B/(B+1))*mtotal;
mcore = (1/(B+1))*mtotal;


rho_19 =  1000*P_19/(Rideal*T_19);
A19 = mbypass /(rho_19*U19);


rho_9 = 1000*P_9/(Rideal*T_9);
A9 = mcore /(rho_19*U9);


rho_o_4 = 1000*Po_4/(Rideal*To_4);
rho_4 = rho_o_4/(1+(gamma_turb -1)/2  ) ^(1/(gamma_turb-1)); % as Ma = 1
T_4 = To_4/(1+(gamma_turb -1)/2  );     % as Ma = 1
U4 = sqrt(gamma_turb*Rideal*T_4); % as Ma = 1
A4 = mcore/(rho_4*U4);



%% take-off

% guess
R_fan_to = 1.2;

%ambient conditions
Pa_to = 101.325;
rhoa_to = 1.225;
Ta_to = 288.15;
%Vtakeoff = [90 0 3] ; % i j k
Ua_to = 90;
ca_to= sqrt(gamma_comp*Rideal*Ua);
Ma_to = Ua_to/ca_to;


%inlet connies
Po_to_2 = Pa_to*(1+((gamma_comp-1)/2)*Ma_to^2)^(gamma_comp/(gamma_comp-1));
To_to_2 = Ta_to + (Ua_to^2)/(2*Cp_comp);


% fan outlet conditions
Po_to_13 = R_fan_to*Po_to_2;
To_to_13 = To_to_2 * R_fan_to^((gamma_comp-1)/(gamma_comp*ec));

% booster outlet conditions
BTR_relator = (To_23 - To_2)/(To_13-To_2);
To_to_23 = BTR_relator*(To_to_13-To_to_2)+To_to_2;
Bootempratio = To_to_23/To_to_13;
Po_to_23 = Po_to_13*Bootempratio^(gamma_comp*ec/(gamma_comp-1));


% combustion chamber temperature conditions
inlet_comb_temp_ratio = To_4/To_2;
To_to_4 = To_to_2 * inlet_comb_temp_ratio;


% HP turbine temperature conditions
HPT_relator = (To_4 - To_45)/(To_4);
To_to_45 = To_to_4-HPT_relator*(To_to_4);
HPTtempratio = To_to_45/To_to_4;

% HPC temp and pressure calc

To_to_3 = (Cp_turb/Cp_comp)*(To_to_4-To_to_45)+To_to_23;
HPCtempratio = To_to_3/To_to_23;
Po_to_3 = Po_to_23*HPCtempratio^(gamma_comp*ec/(gamma_comp-1));

% combustion and high pressure turbine pressures at outlets
Po_to_4 = Po_to_3;

Po_to_45 = Po_to_4* (To_to_45/To_to_4)^(gamma_comp/((gamma_comp-1)*ec));

% calculating mass flow in core

rho_to_o_4 = 1000*Po_to_4/(Rideal*To_to_4);
rho_to_4 = rho_to_o_4/(1+(gamma_turb -1)/2  ) ^(1/(gamma_turb-1)); % as Ma = 1
T_to_4 = To_to_4/(1+(gamma_turb -1)/2  );     % as Ma = 1
U4 = sqrt(gamma_turb*Rideal*T_to_4);
mcore_to_4 = rho_to_4 * A4 * T_to_4;
















