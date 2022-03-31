


%%

%ambient conditions
Pa = 23.842;
rhoa = 0.3799;
Ta = 218.9;
Ma = 0.78;
Ua = 231.3;

%other properties
gamma_comp = 1.4;
gamma_turb = 1.3;
Cp_comp = 1005;
Cp_turb = 1244;
Rideal = 287;
ec = 0.9;
Fengine = 38900;
%compressor ratios
R(2) = 1.5; %fan
R(3) = 1.667; %booster
R(4) = 18; %bypass

% Inlet Conditions
P(1) = Pa*(1+((gamma_comp-1)/2)*Ma^2)^(gamma_comp/(gamma_comp-1));
T(1) = Ta + (Ua^2)/(2*Cp_comp);

%compressor stages: 13, 23, 3
 i = 2;
P(i) = R(i)*P(i-1);
T(i) = T(i-1)*R(i)^((gamma_comp-1)/(gamma_comp*ec));


i = 3;
P(i) = R(i)*P(i-1);
T(i) = T(i-1)*R(i)^((gamma_comp-1)/(gamma_comp*ec));
%T(3) = R(i)^((gamma_comp-1)/(gamma_comp*ec))*(T(2)-T(1))+T(1);



i = 4;
P(i) = R(i)*P(i-1);
T(i) = T(i-1)*R(i)^((gamma_comp-1)/(gamma_comp*ec));
%combustion chamber exit conditions: 4
P(5)= P(4);
T(5)= 1500;

% HP turbine outlet conditions: 45
T(6) = T(5) - (Cp_comp/Cp_turb)* (T(4)-T(3));
HPTtempratio = (T(5)-T(6))/T(5);
P(6) = P(5)*(T(6)/T(5))^((gamma_turb)/(ec*(gamma_turb-1)));

B = 13.74;

%LP turbine outlet conditions: 5
T(7) = T(6)-(Cp_comp/Cp_turb)*( T(3)-T(2) + (B+1)*(T(2)-T(1)));
P(7) = P(6)*(T(7)/T(6))^((gamma_turb)/(ec*(gamma_turb-1)));


ChokedRatioCore = ((gamma_turb+1)/2)^(gamma_turb/(gamma_turb-1));
CoreexitRatio = P(7)/Pa;

if CoreexitRatio>=ChokedRatioCore

    disp('core is choked')
    % 8 refers to section 9 of the engine
    T(8) = T(7)*(1/ChokedRatioCore)^((gamma_turb-1)/gamma_turb);
    Uexitcore = sqrt(gamma_turb*Rideal*T(8))
else
    disp('core is not choked')
    P(8) = Pa;
    T(8) = T(7)*(P(8)/P(7))^((gamma_turb-1)/gamma_turb);
    %Mexitcore = (CoreexitRatio^((gamma_turb-1)/gamma_turb)-1)*(2/(gamma_turb-1));
    %Uexitcore = Mexitcore*sqrt(gamma_turb*Rideal*T(8))
    Uexitcore = sqrt(2 * Cp_turb *  (T(7)-T(8))   )
end

ChokedRatioBypass = ((gamma_comp+1)/2)^(gamma_comp/(gamma_comp-1));
BypassexitRatio = P(2)/Pa;
P(9) = P(2)/ChokedRatioBypass;
if BypassexitRatio>=ChokedRatioBypass

    disp('Bypass is choked')
    % 8 refers to section 9 of the engine
    T(9) = T(2)*(1/ChokedRatioBypass)^((gamma_comp-1)/gamma_comp);
    Uexitbypass = sqrt(gamma_comp*Rideal*T(9))
end



%results

Rhpt = P(6)/P(5);
MaBypassexit = 1;
Rlpt = P(7)/P(6);

Worklpt= Cp_turb*(T(6)-T(7));
Worklptonbypass = Cp_comp* B*(T(2)-T(1));
proportionWorklpt_onBypass = Worklptonbypass/Worklpt;

FspecCore = (1/(1+B))*(Uexitcore-Ua);
FspecBypass = (B/(1+B))*(Uexitcore-Ua)+(B/(1+B))* ( (Rideal*T(9))/(P(9)*Uexitbypass) ) * ( P(9)- Pa);
Fspec = FspecBypass+FspecCore;
mtotal = Fengine/Fspec;

mcore = mtotal*(1/(1+B));
mbypass = mtotal*(B/(1+B));

Fbypass = mtotal*FspecBypass;
AbypassExit = (Fbypass-mbypass*(Uexitbypass-Ua))/((P(9)-Pa)*1000);


for i = 1:9
rho(i) = (P(i)*1000)/(Rideal*T(i));

end

%rhostatic(7) = rho(7)/( (1+(gamma_turb-1)/2)^(1/(gamma_turb-1))   )


Acore= mcore/(rho(8)*Uexitcore);
V4 = sqrt(gamma_turb*Rideal*T(5));
A4 = mcore/(rho(5)*V4);



%% take-off
Pa2 = 101.325;
rhoa2 = 1.225;
Ta2 = 288.15;
Vtakeoff = [90 0 3] ; % i j k
Ua2 = norm(Vtakeoff);
speedofsound= sqrt(gamma_comp*Rideal*Ua);
Ma2 = Ua2/speedofsound;

%inlet conditions
Pto(1)= Pa*(1+((gamma_comp-1)/2)*Ma2^2)^(gamma_comp/(gamma_comp-1));
Tto(1) = Ta2 + (Ua2^2)/(2*Cp_comp);

Rto(2) = 1.6;

%compressor stages: 13, 23, 3
 i = 2;
Pto(i) = Rto(i)*Pto(i-1);
Tto(i) = Tto(i-1)*Rto(i)^((gamma_comp-1)/(gamma_comp*ec));

T4toT2ratio = T(5)/T(1);
Tto(5) = Tto(1)*T4toT2ratio;
tauHPT = (T(5)-T(6))/T(5);
Tto(6) = tauHPT*Tto(5);


tauBooster = (T(3)-T(1))/(T(2)-T(1));
Tto(3) =tauBooster*Tto(2);
Pto(3) = Pto(2)*tauBooster^(gamma_comp*ec/(gamma_comp-1));

Tto(4)= Tto(3)+ (Cp_turb/Cp_comp)*(Tto(6)-Tto(5));

% index 5 is equalt to hpt inlet
rhotostag5 = 1000*P(5)/(Rideal*Tto(5))
rho5to = rhotostag5/(1+(gamma_turb -1)/2  ) ^(); % assuming hpt inlet is choked