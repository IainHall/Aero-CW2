


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

B = 13.44;

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

if BypassexitRatio>=ChokedRatioBypass

    disp('Bypass is choked')
    % 8 refers to section 9 of the engine
    P(9) = Pa;
    T(9) = T(2)*(P(9)/P(2))^((gamma_comp-1)/gamma_comp);
    %Uexitbypass = sqrt(gamma_comp*Rideal*T(9))
    speedofsoundbypassexit = sqrt(gamma_comp*Rideal*T(9));
    Uexitbypass = sqrt(2 * Cp_comp *  (T(2)-T(9))   )


end


%results

Rhpt = P(6)/P(5);
MaBypassexit = Uexitbypass/speedofsoundbypassexit;
Rlpt = P(7)/P(6);

Worklpt= Cp_turb*(T(6)-T(7));
Worklptonbypass = Cp_comp* B*(T(2)-T(1));
proportionWorklpt_onBypass = Worklptonbypass/Worklpt;




