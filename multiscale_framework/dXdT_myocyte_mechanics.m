%%%% Differential equations for myocyte model

% Adapted by Abby Teitgen based on code from: 

% Tewari, S. G., Bugenhagen, S. M., Palmer, B. M, Beard,
% D. A. (2016). Dynamics of corss-bridge cycling, ATP hydrolysis, force
% generation, and deformation in cardiac muscle. J. Mol. Cell Cardiol. 96:
% 11-25.

% Tewari, S. G., Bugenhagen, S. M., Vinnakota, K. C., Rice, J. J., Janssen,
% M. L., Beard, D. A. (2016). Influence of metabolic dysfuntion on cardiac
% mechanics in decompensated hypertrophy and heart failure. J. Mol.
% Cell Cardiol. 94: 162-175.

% Lopez, R. Marzban, B., Gao, X. Lauinger, E., Van den Bergh, F.,
% Whitesall, S. E., Converso-Baran, K., Burant, C. F., Michele, D. E.,
% Beard, D. A. (2020). Impaired myocardial energetics causes mechanical dysfunction
% in decompensated failing hearts. Function, 1(2): zqaa018

% Marzban, B., Lopez, R., Beard, D. A. (2020). Computational modeling of coupled
% energetics and mechanics in the rat ventricular myocardium. Physiome.

function [dYdT, dSL, F_XB, F_passive] = dXdT_myocyte_mechanics(t, y, para, Ca_fraction, flag, Ca_flag)
%% Model parameters
% Metabolite concentrations based on mean sham rat (Lopez et al. 2020)
MgATP = para(1); % mM
MgADP = para(2); % mM
Pi = para(3); % mM

% Crossbridge parameters
kstiff1 = para(4); % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = para(5); % Stiffness constant due to working stroke of XBs (kPa/um)
k_passive = para(6); % Passive stiffness constant (kPa/um)
Kse = para(7); % Series element stiffness (mmHg/um) 
k_coop = para(8); % Strength of thin filament cooperativity
k_on = para(9); % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = para(10); % Rate constant of Ca2+ unbinding from troponin C (s^-1)
km = para(11); % OFF to ON transition rate (s^-1)
krecruit = para(12); % Force dependence of transition to OFF state (N^-1 m^-1)
k_m = para(13); % ON to OFF transition rate (s^-1)
kf = para(14); % Myosin actin associaiton rate (P to A1) (s^-1)
k_f = para(15); % Myosin actin dissociation rate (A1 to P) (s^-1) 
kw = para(16); % A1 to A2 transition rate (s^-1)
k_w = para(17); % A2 to A1 transition rate (s^-1) 
kp = para(18); % A2 to A3 transition rate (s^-1)
k_p = para(19); % A3 to A2 transition rate (s^-1) 
kg = para(20); % A3 to P transition rate (s^-1)
visc = para(21); % Viscosity (mmHg*s/um)
stim_period = para(22); % Hz

% Calcium transient
a_Ca = para(23);
b_Ca = para(24);
c_Ca = para(25);

L_rest_pas = 1.51; % Length at which passive force = 0 (um)
alpha1 = 10; % Stretch sensing parameter for kw and k_w (1/um)
alpha2 = 9.1; % Stretch sensing parameter for kp and k_p
alpha3 = 5.93; % Stretch sensing parameter for kg
s3 = 9.9e-3; % Strain at which kg is minimum (um)
K_Pi = 4.00; % Pi dissociation constant (mM)
K_T = 0.4897 ; % MgATP dissociation constant (mM)
K_D = 0.194;  % MgADP dissociation constant (mM)
Lthin = 1.200; % Length of thin filament (nm)
Lthick = 1.670; % Length of thick filament (nm)
Lbare = 0.100; % Bare length of thin filament (nm)
dr = 0.01; % Power-stroke size (um)

% Defining the metabolite dependent coeficient, rapid equilibrium of the
% cross bridge sub-states
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T);
g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi);
f2 = 1/(1 + Pi/K_Pi);

% Adjust for metabolite concentrations
k_f = k_f*f1; 
kw = kw*f2; 
k_p = k_p*g1; 
kg = kg*g2;

if flag == 0 % force pCa
SLset = 2.25; 
Ca_i = Ca_fraction;

else if flag == 1 % Twitch
SLset = 1.84;
        if Ca_flag == 0 % ATP
            % Korte
            %a = 0.106;
            %b = 0.5635;
            %c = 1.8017;
            %Ca0 = 0;
            
            % Average
            %a = 0.1208;
            %b = 0.6651;
            %c = 1.7374;
            a = a_Ca;
            b = b_Ca;
            c = c_Ca;
            Ca0 = 0;
            
        else % dATP
            % Korte
            %a = 0.0534;
            %b = 0.5484;
            %c = 2.4732;
            %Ca0 = 0;
            
            % KM 
            %a = 0.1114;    
            %b = 0.8391;   
            %c = 1.8965;
            %Ca0 = 0;
            
            % Average
            %a = 0.0864;
            %b = 0.5971;
            %c = 2.0301;
            a = a_Ca;
            b = b_Ca;
            c = c_Ca;
            Ca0 = 0;
            
        end
        phi = mod(t+0.001, stim_period)/stim_period;
        Ca_i = ((a/phi)*exp(-b*(log(phi)+c)^2) + Ca0);
end
end


%% State variables
% Moments of state variables
p1_0 = y(1);
p1_1 = y(2);
p1_2 = y(3);
p2_0 = y(4);
p2_1 = y(5);
p2_2 = y(6);
p3_0 = y(7);
p3_1 = y(8);
p3_2 = y(9);
N = y(10); % Non-permissible state 
P = 1.0 - N - p1_0 - p2_0 - p3_0; % Permssible state 
SL = y(11); % Sarcomere length (um)
U_NR = y(12); % ON state
U_SR = 1.0 - U_NR; % OFF state

%% Stretch-sensitive rates
f_alpha1o = (p1_0 - alpha1*p1_1 + 0.5*(alpha1*alpha1)*p1_2);
f_alpha1i = (p1_1 - alpha1*p1_2);

f_alpha0o = (p2_0 + alpha1*p2_1 + 0.5*alpha1*alpha1*p2_2);
f_alpha0i = (p2_1 + alpha1*p2_2);

f_alpha2o = (p2_0 - alpha2*p2_1 + 0.5*(alpha2*alpha2)*p2_2);
f_alpha2i = (p2_1 - alpha2*p2_2);

f_alpha3o = (p3_0 + alpha3*(s3*s3*p3_0 + 2.0*s3*p3_1 + p3_2)); 
f_alpha3i = (p3_1 + alpha3*(s3*s3*p3_1 + 2.0*s3*p3_2));

%% Compute active and passive force
% Overlap function 
OV_Zaxis = min(Lthick/2, SL/2); % Overlap region closest to Z-axis (nm)
OV_Mline = max(SL/2 - (SL - Lthin), Lbare/2); % Overal region closest to M-line (nm)
LOV = OV_Zaxis - OV_Mline; % Length of overlap (nm)
N_overlap_thick = LOV*2/(Lthick - Lbare); % Fraction of thick filament overlap

% Active Force
B_process = kstiff2*dr*p3_0;   % Force due to XB ratcheting
C_process = kstiff1*(p2_1 + p3_1 );% Force due to stretching of XBs
F_XB = N_overlap_thick*(B_process + C_process); % Active force

% Passive force
gamma = 8;
F_passive = k_passive*(SL - L_rest_pas)^gamma; 

% Total force
Ftotal = F_XB + F_passive;

%% Simulate various experimental protocols
if flag == 0 % Force pCa
%dSL = 0; % Fully isometric
intf = (- Ftotal  + Kse*(SLset - SL));
dSL = intf/visc; % Muscle isometric but not sarcomere isometric

else if flag == 1 % Twitch
intf = (-Ftotal + Kse*(SLset - SL));
%intf = -Ftotal; % Fully unloaded shortening
dSL = intf/visc;
%dSL = 0; % Fully isometric
end
end

dp1_0 = kf*P*N_overlap_thick*U_NR - k_f*p1_0 - kw*f_alpha1o + k_w*f_alpha0o;
dp1_1 = 1*dSL*p1_0 - k_f*p1_1 - kw*f_alpha1i + k_w*f_alpha0i;
dp1_2 = 2*dSL*p1_1 - k_f*p1_2 - kw*p1_2 + k_w*p2_2;
dp2_0 = - k_w*f_alpha0o - kp*f_alpha2o + k_p*p3_0 + kw*f_alpha1o;
dp2_1 = 1*dSL*p2_0 - k_w*f_alpha0i - kp*f_alpha2i + k_p*p3_1 + kw*f_alpha1i;
dp2_2 = 2*dSL*p2_1 - k_w*p2_2 - kp*p2_2 + k_p*p3_2 + kw*p1_2;
dp3_0 = + kp*f_alpha2o - k_p*p3_0 - kg*f_alpha3o;
dp3_1 = 1*dSL*p3_0 + kp*f_alpha2i - k_p*p3_1 - kg*f_alpha3i;
dp3_2 = 2*dSL*p3_1 + kp*p2_2 - k_p*p3_2 - kg*p3_2;


%% Campbell Ca activation model (Campbell et al. 2018)
Jon = k_on*Ca_i*N*(1 + k_coop*(1 - N));
Joff = k_off*P*(1 + k_coop*N);
dNp = -Jon + Joff;   

%% Transitions between OFF and ON states
dU_NR = km * (1 + krecruit * F_XB) * U_SR - k_m * U_NR; % ON state
dYdT = [dp1_0; dp1_1; dp1_2; dp2_0; dp2_1; dp2_2; dp3_0; dp3_1; dp3_2; dNp; dSL; dU_NR];

end