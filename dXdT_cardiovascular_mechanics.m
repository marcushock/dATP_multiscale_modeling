%%%% Differential equations for ventricular model

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

function dXdT = dXdT_cardiovascular_mechanics(t,x,input)

%% Parameters
% Heart 
stim_period = input(1); % Stimulation period

C_Ao = input(2);  % Proximal aortic compliance (mL/mmHg)
C_SA = input(3); % Systemic arterial compliance (mL/mmHg)
C_SV = input(4); % Systemic venous compliance (mL/mmHg) 
C_PV = input(5); % Pulmonary venous compliance (mL/mmHg)
C_PA = input(6); % Pulmonary arterial compliance (mL/mmHg)
R_Ao = input(7); % Aortic resistance (mmHg*sec/mL)
R_SA = input(8); % Systemic vasculature resistance (mmHg*sec/mL)
R_PA = input(9); % Pulmonary arterial resistance (mmHg*sec/mL)
R_SV = input(10); % Systemic vascular resistance (mmHg*sec/mL)
R_PV = input(11); % Pulmonary venous resistance (mmHg*sec/mL)
R_vlv = input(12); % Valve resistance (mmHg*sec/mL)
R_AV = input(13); % Resistance across aortic valve (mmHg*sec/mL)
R_tAo = input(14); % Transmural aortic resistance (mmHg*sec/mL)
R_tSA = input(15); % Transmural systemic arterial resistance (mmHg*sec/mL)

% Geometry
Amref_LV = input(16); % LV midwall reference surface area (cm^2)
Amref_SEP = input(17); % SEP midwall reference surface area (cm^2)
Amref_RV = input(18); % RV midwall reference surface area (cm^2)

Vw_LV = input(19); % LV wall volume (mL)
Vw_SEP = input(20); % SEP wall volume (mL)
Vw_RV = input(21); % RV wall volume (mL)

% Energetics 
MgATP_LV = input(22); % Cytosolic MgATP, LV (M/L cytosol water)
MgADP_LV = input(23); % Cytosolic MgADP, LV (M/L cytosol water)
Pi_LV = input(24); % Cytosolic Pi, LV (M/L cytosol water)
MgATP_SEP = input(25); % Cytosolic MgATP, SEP (M/L cytosol water)
MgADP_SEP = input(26); % Cytosolic MgADP, SEP (M/L cytosol water)
Pi_SEP = input(27); % Cytosolic Pi, SEP (M/L cytosol water)
MgATP_RV = input(28); % Cytosolic MgATP, RV (M/L cytosol water)
MgADP_RV = input(29); % Cytosolic MgADP, RV (M/L cytosol water)
Pi_RV = input(30); % Cytosolic Pi, RV (M/L cytosol water)

% Crossbridge
Lsref = input(31); % Resting sarcomere length (um)
L_rest_pas = input(32); % Length at which passive force = 0 (um)
gamma = input(33); % For calculating passive force

K_Pi = input(34); 
K_T = input(35); 
K_D = input(36);

alpha1 = input(37); % Stretch sensing parameter for kw and k_1, 1/um
alpha2 = input(38); % Stretch sensing parameter for kp and k_2, 1/um
alpha3 =  input(39); % Stretch sensing parameter for k_3, 1/um
s3 = input(40);  % Strain at which kg is minimum, um

Kse = input(41); % Series stiffness (mmHg/um)
k_passive_LV = input(42); % Passive stiffness (mmHg/um)
k_passive_SEP = input(43);
k_passive_RV = input(44);
eta = input(45); % Viscosity (mmHg*s/um)

L_thick = input(46); % Length of thick filament (um)
L_hbare = input(47); % Length of bare region of thick filament (um)
L_thin  = input(48); % Length of thin filament (um)
deltaR  = input(49); % um

kstiff1 = input(50); % Stiffness constant due to myosin actin interaction (kPa/um)
kstiff2 = input(51); % Stiffness constant due to working stroke of XBs (kPa/um)

kf = input(52); % P to A1 transition rate (1/s)
k_f = input(53);% A1 to P transition rate (1/s)
kw = input(54); % A1 to A2 transition rate (1/s)
k_w = input(55); % A2 to A1 transition rate (1/s)
kp = input(56); % A2 to A3 transition rate (1/s)
k_p = input(57); % A3 to A2 transition rate (1/s)
kg = input(58); % A3 to P transition rate (1/s)
k_coop = input(59); % Cooperative constant
k_on = input(60); % N to P transition rate (1/Ms)
k_off = input(61); % P to N transition rate (1/s)
km = input(62); % OFF to on transition rate transition rate (1/s)
krecruit = input(63);  % Force dependent constant (1/Nm^s)
k_m = input(64); % ON to OFF transition rate (1/s)

% Correcting rate constants for metabolite levels 
k_f_LV  = k_f*(Pi_LV/K_Pi)/(1.0 + Pi_LV/K_Pi);
kw_LV  = kw/(1.0 + Pi_LV/K_Pi);
k_p_LV = k_p*(MgADP_LV/K_D)/(1.0 + MgADP_LV/K_D + MgATP_LV/K_T);
kg_LV  = kg*(MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);

k_f_SEP  = k_f*(Pi_SEP/K_Pi)/(1.0 + Pi_SEP/K_Pi);
kw_SEP  = kw/(1.0 + Pi_SEP/K_Pi);
k_p_SEP = k_p*(MgADP_SEP/K_D)/(1.0 + MgADP_SEP/K_D + MgATP_SEP/K_T);
kg_SEP  = kg*(MgATP_SEP/K_T)/(1.0 + MgATP_SEP/K_T + MgADP_SEP/K_D);

k_f_RV  = k_f*(Pi_RV/K_Pi)/(1.0 + Pi_RV/K_Pi);
kw_RV  = kw/(1.0 + Pi_RV/K_Pi);
k_p_RV = k_p*(MgADP_RV/K_D)/(1.0 + MgADP_RV/K_D + MgATP_RV/K_T);
kg_RV  = kg*(MgATP_RV/K_T)/(1.0 + MgATP_RV/K_T + MgADP_RV/K_D);

dATP = input(69);

Ca_flag = input(70);

% Calcium
if Ca_flag == 0 % Always use ATP Ca transient
    % Orig
    a = input(65);
    b = input(66);
    c = input(67);
    Ca0 = input(68);
    
else if Ca_flag == 1 % Use dATP Ca transient if dATP > 0
        if dATP == 0
            % Orig
            a = 0.2968;
            b = 1.5900;
            c = 1.7171;
            Ca0 = 0.1622;
        else % dATP > 0
            % Korte
            %a = 0.1447;
            %b = 2.1011;
            %c = 2.4775;
            %Ca0 = 0.1622;

            % KM
            %a = 0.2864;
            %b = 1.9679;
            %c = 1.7832;
            %Ca0 = 0.1622;
            
            % Average
            %a = 0.2218;
            %b = 1.4980;
            %c = 1.9989;
            %Ca0 = 0.1622;
            a = input(65);
            b = input(66);
            c = input(67);
            Ca0 = input(68);
        end
    end
end

%% State variables
xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm
SL_LV  = x(5); % sarcomere length, LV, micron
SL_SEP = x(6); % sarcomere length, septum, micron
SL_RV  = x(7); % sarcomere length, RV, micron
V_LV   = x(8); % volume LV, mL
V_RV   = x(9); % volume RV, mL

P1_0_LV = x(10); % 0th moment state A1, LV
P1_1_LV = x(11); % 1st moment state A1, LV
P1_2_LV = x(12); % 2nd moment state A1, LV
P2_0_LV = x(13); % 0th moment state A2, LV
P2_1_LV = x(14); % 1st moment state A2, LV
P2_2_LV = x(15); % 2nd moment state A2, LV
P3_0_LV = x(16); % 0th moment state A3, LV
P3_1_LV = x(17); % 1st moment state A3, LV
P3_2_LV = x(18); % 2nd moment state A3, LV
N_LV    = x(19); % Non permissible state, LV
U_NR_LV = x(20); % ON state, LV

P1_0_SEP = x(21); % 0th moment state A1, SEP
P1_1_SEP = x(22); % 1st moment state A1, SEP
P1_2_SEP = x(23); % 2nd moment state A1, SEP
P2_0_SEP = x(24); % 0th moment state A2, SEP
P2_1_SEP = x(25); % 1st moment state A2, SEP
P2_2_SEP = x(26); % 2nd moment state A2, SEP
P3_0_SEP = x(27); % 0th moment state A3, SEP
P3_1_SEP = x(28); % 1st moment state A3, SEP
P3_2_SEP = x(29); % 2nd moment state A3, SEP
N_SEP    = x(30); % Non permissible state, SEP
U_NR_SEP = x(31); % ON state, SEP

P1_0_RV = x(32); % 0th moment state A1, RV
P1_1_RV = x(33); % 1st moment state A1, RV
P1_2_RV = x(34); % 2nd moment state A1, RV
P2_0_RV = x(35); % 0th moment state A2, RV
P2_1_RV = x(36); % 1st moment state A2, RV
P2_2_RV = x(37); % 2nd moment state A2, RV
P3_0_RV = x(38); % 0th moment state A3, RV
P3_1_RV = x(39); % 1st moment state A3, RV
P3_2_RV = x(40); % 2nd moment state A3, RV
N_RV    = x(41); % Non permissible state, RV
U_NR_RV = x(42); % ON state, RV

V_SV   = x(43); % Systemic veins volume (mL)
V_PV   = x(44); % Pulmonary veins volume (mL)
V_SA   = x(45); % Systemic arterys volume (mL)
V_PA   = x(46); % Pulmonary arterys volume (mL)
V_Ao   = x(47); % Proximal aorta volume (mL)

%% Calcium
phi = mod(t+0.0001,stim_period)/stim_period;
Ca_i = (a/phi)*exp(-b*(log(phi)+c)^2) + Ca0;

%% Heart and Sarcomere Model
% Geometry
Vm_LV  = (pi/6)*xm_LV*(xm_LV^2 + 3*ym^2); % LV midwall volume (mL)
Vm_SEP = (pi/6)*xm_SEP*(xm_SEP^2 + 3*ym^2); % SEP midwall volume (mL)
Vm_RV  = (pi/6)*xm_RV*(xm_RV^2 + 3*ym^2); % RV midwall volume (mL)
Am_LV  = pi*(xm_LV^2 + ym^2); % LV midwall surface area (cm^2)
Am_SEP = pi*(xm_SEP^2 + ym^2); % SEP midwall surface area (cm^2)
Am_RV  = pi*(xm_RV^2 + ym^2); % RV midwall surface area (cm^2)
Cm_LV  = 2*xm_LV/(xm_LV^2 + ym^2); % LV midwall surface curvature
Cm_SEP = 2*xm_SEP/(xm_SEP^2 + ym^2); % SEP midwall surface curvature
Cm_RV  = 2*xm_RV/(xm_RV^2 + ym^2); % RV midwall surface curvature
z_LV   = 3*Cm_LV*Vw_LV/(2*Am_LV);
z_SEP  = 3*Cm_SEP*Vw_SEP/(2*Am_SEP);
z_RV   = 3*Cm_RV*Vw_RV/(2*Am_RV);

% Fiber strain
epsf_LV = 0.5*log(Am_LV/Amref_LV) - 0.083333*z_LV^2 - 0.019*z_LV^4;
epsf_SEP = 0.5*log(Am_SEP/Amref_SEP) - 0.083333*z_SEP^2 - 0.019*z_SEP^4;
epsf_RV = 0.5*log(Am_RV/Amref_RV) - 0.083333*z_RV^2 - 0.019*z_RV^4;
SLo_LV = Lsref*exp(epsf_LV); 
SLo_SEP = Lsref*exp(epsf_SEP); 
SLo_RV = Lsref*exp(epsf_RV);

% Passive force
sigmapas_LV  = k_passive_LV*(SLo_LV - L_rest_pas)^(gamma);
sigmapas_SEP  = k_passive_SEP*(SLo_SEP - L_rest_pas)^(gamma);
sigmapas_RV  = k_passive_RV*(SLo_RV - L_rest_pas)^(gamma);

% Active force
sovr_ze = min(L_thick*0.5, SL_LV*0.5); % Overlap region closest to Z-axis (nm)
sovr_cle = max(SL_LV*0.5 - (SL_LV-L_thin),L_hbare*0.5); % Overlap region closes to M-line (nm)
L_sovr = sovr_ze - sovr_cle; % Length of overlap (nm)
N_overlap_LV = L_sovr*2/(L_thick - L_hbare); % Fraction of thick filament overlap
 
sovr_ze = min(L_thick*0.5, SL_SEP*0.5); % Overlap region closest to Z-axis (nm)
sovr_cle = max(SL_SEP*0.5 - (SL_SEP-L_thin),L_hbare*0.5); % Overlap region closes to M-line (nm)
L_sovr = sovr_ze - sovr_cle; % Length of overlap (nm)
N_overlap_SEP = L_sovr*2/(L_thick - L_hbare); % Fraction of thick filament overlap

sovr_ze = min(L_thick*0.5, SL_RV*0.5); % Overlap region closest to Z-axis (nm)
sovr_cle = max(SL_RV*0.5 - (SL_RV-L_thin),L_hbare*0.5); % Overlap region closes to M-line (nm)
L_sovr = sovr_ze - sovr_cle; % Length of overlap (nm)
N_overlap_RV = L_sovr*2/(L_thick - L_hbare); % Fraction of thick filament overlap

sigmaact_LV  = N_overlap_LV*(kstiff2*deltaR*(P3_0_LV) + kstiff1*(P2_1_LV + P3_1_LV)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_SEP = N_overlap_SEP*(kstiff2*deltaR*(P3_0_SEP) + kstiff1*(P2_1_SEP + P3_1_SEP)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_RV  = N_overlap_RV*(kstiff2*deltaR*(P3_0_RV) + kstiff1*(P2_1_RV + P3_1_RV)); % Includes force due to XB ratcheting and stretching of XBs

% Total force
sigmaf_LV = -Kse*(SL_LV - SLo_LV);
sigmaf_SEP = -Kse*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse*(SL_RV - SLo_RV);

% Tension
Tm_LV  = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5); % LV midwall tension
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5); % SEP midwall tension
Tm_RV  = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5); % RV midwall tension

Tx_LV  = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2); % LV axial tension
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2); % SEP axial tension
Tx_RV  = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2); % RV axial tension

Ty_LV  = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2); % LV radial tension
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2); % SEP radial tension
Ty_RV  = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2); % RV radial tension

% Ventricular pressures
ptrans1 = 2*Tx_LV/ym;
ptrans3 = 2*Tx_RV/ym;
P_LV = -ptrans1; % LV transmural pressure (mmHg)
P_RV = ptrans3; % RV transmural pressure (mmHg)

% Pulmonary pressures
P_SV = V_SV/C_SV; % Systemic venous pressure (mmHg)
P_PV = V_PV/C_PV; % Pulmonary venous pressure (mmHg)
P_PA = V_PA/C_PA; % Pulmonary arterial pressure (mmHg)

% Lumped circulatory model
% Ao valve closed equations
QOUT_LV = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV)*(V_LV>0) 

% Ao valve open equations 
  P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_AV*V_SA + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_AV + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
  QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
  Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_AV*V_SA - C_SA*R_SA*R_AV*V_Ao - C_SA*R_tSA*R_AV*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
end

QIN_LV  = max((P_PV - P_LV)/R_PV,0); 
QIN_RV  = max((P_SV - P_RV)/R_SV,0); 
QOUT_RV = max((P_RV - P_PA)/R_vlv,0)*(V_RV>0); 

% TriSeg
dXdT(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV;
dXdT(2) = (Tx_LV + Tx_SEP + Tx_RV); 
dXdT(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; 
dXdT(4) = (Ty_LV + Ty_SEP + Ty_RV);  

% Myofiber Mechanics: SL_LV
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_LV - alpha1*P1_1_LV + 0.5*(alpha1*alpha1)*P1_2_LV);
f_alpha1i = (P1_1_LV - alpha1*P1_2_LV);

f_alpha0o = (P2_0_LV + alpha1*P2_1_LV + 0.5*alpha1*alpha1*P2_2_LV);
f_alpha0i = (P2_1_LV + alpha1*P2_2_LV);

f_alpha2o = (P2_0_LV - alpha2*P2_1_LV + 0.5*(alpha2*alpha2)*P2_2_LV);
f_alpha2i = (P2_1_LV - alpha2*P2_2_LV);

f_alpha3o = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV));
f_alpha3i = (P3_1_LV + alpha3*(s3*s3*P3_1_LV + 2.0*s3*P3_2_LV));

dSL_LV = (sigmaf_LV - sigmapas_LV - sigmaact_LV)/eta;

P0_LV = 1.0 - N_LV - P1_0_LV - P2_0_LV - P3_0_LV; 
dXdT(10) = kf*P0_LV*U_NR_LV*N_overlap_LV - k_f_LV*P1_0_LV - kw_LV*f_alpha1o + k_w*f_alpha0o;
dXdT(11) = dSL_LV*P1_0_LV - k_f_LV*P1_1_LV - kw_LV*f_alpha1i + k_w*f_alpha0i;
dXdT(12) = 2*dSL_LV*P1_1_LV - k_f_LV*P1_2_LV - kw_LV*P1_2_LV + k_w*P2_2_LV;

dXdT(13) = -k_w*f_alpha0o - kp*f_alpha2o + k_p_LV*P3_0_LV + kw_LV*f_alpha1o;
dXdT(14) = dSL_LV*P2_0_LV - kw_LV*f_alpha0i - kp*f_alpha2i + k_p_LV*P3_1_LV + kw_LV*f_alpha1i;
dXdT(15) = 2*dSL_LV*P2_1_LV - k_w*P2_2_LV - kp*P2_2_LV + k_p_LV*P3_2_LV + kw_LV*P1_2_LV;

dXdT(16) = +kp*f_alpha2o - k_p_LV*P3_0_LV - kg_LV*f_alpha3o;
dXdT(17) = dSL_LV*P3_0_LV + kp*f_alpha2i - k_p_LV*P3_1_LV - kg_LV*f_alpha3i;
dXdT(18) = 2*dSL_LV*P3_1_LV + kp*P2_2_LV - k_p_LV*P3_2_LV - kg_LV*P3_2_LV;

U_SR_LV = 1 - U_NR_LV;
Jon = k_on*Ca_i*N_LV*(1 + k_coop*(1 - N_LV));
Joff = k_off*P0_LV*(1 + k_coop*N_LV);
dXdT(19) = - Jon + Joff; 
dXdT(20) = km * (1 + krecruit * sigmaact_LV) * U_SR_LV - k_m*U_NR_LV ; 

% Myofiber Mechanics: SL_SEP
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_SEP - alpha1*P1_1_SEP + 0.5*(alpha1*alpha1)*P1_2_SEP);
f_alpha1i = (P1_1_SEP - alpha1*P1_2_SEP);

f_alpha0o = (P2_0_SEP + alpha1*P2_1_SEP + 0.5*alpha1*alpha1*P2_2_SEP);
f_alpha0i = (P2_1_SEP + alpha1*P2_2_SEP);

f_alpha2o = (P2_0_SEP - alpha2*P2_1_SEP + 0.5*(alpha2*alpha2)*P2_2_SEP);
f_alpha2i = (P2_1_SEP - alpha2*P2_2_SEP);

f_alpha3o = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP));
f_alpha3i = (P3_1_SEP + alpha3*(s3*s3*P3_1_SEP + 2.0*s3*P3_2_SEP));

dSL_SEP = (sigmaf_SEP - sigmapas_SEP - sigmaact_SEP)/eta;
 
P0_SEP = 1 - N_SEP - P1_0_SEP - P2_0_SEP - P3_0_SEP;  
dXdT(21) = kf*P0_SEP*U_NR_SEP*N_overlap_SEP - k_f_SEP*P1_0_SEP - kw_SEP*f_alpha1o + k_w*f_alpha0o; 
dXdT(22) = dSL_SEP*P1_0_SEP - k_f_SEP*P1_1_SEP - kw_SEP*f_alpha1i + k_w*f_alpha0i;
dXdT(23) = 2*dSL_SEP*P1_1_SEP - k_f_SEP*P1_2_SEP - kw_SEP*P1_2_SEP + k_w*P2_2_SEP;

dXdT(24) = -k_w*f_alpha0o - kp*f_alpha2o + k_p_SEP*P3_0_SEP + kw_SEP*f_alpha1o;
dXdT(25) = dSL_SEP*P2_0_SEP - kw_SEP*f_alpha0i - kp*f_alpha2i + k_p_SEP*P3_1_SEP + kw_SEP*f_alpha1i;
dXdT(26) = 2*dSL_SEP*P2_1_SEP - k_w*P2_2_SEP       - kp*P2_2_SEP + k_p_SEP*P3_2_SEP + kw_SEP*P1_2_SEP;

dXdT(27) = +kp*f_alpha2o - k_p_SEP*P3_0_SEP - kg_SEP*f_alpha3o;
dXdT(28) = dSL_SEP*P3_0_SEP + kp*f_alpha2i - k_p_SEP*P3_1_SEP - kg_SEP*f_alpha3i;
dXdT(29) = 2*dSL_SEP*P3_1_SEP + kp*P2_2_SEP       - k_p_SEP*P3_2_SEP - kg_SEP*P3_2_SEP;

U_SR_SEP = 1 - U_NR_SEP;
Jon = k_on*Ca_i*N_SEP*(1 + k_coop*(1 - N_SEP));
Joff = k_off*P0_SEP*(1 + k_coop*N_SEP);
dXdT(30) = - Jon + Joff; 
dXdT(31) = km*(1 + krecruit * sigmaact_SEP) * U_SR_SEP - k_m*U_NR_SEP ; 

% Myofiber Mechanics: SL_RV
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_RV - alpha1*P1_1_RV + 0.5*(alpha1*alpha1)*P1_2_RV);
f_alpha1i = (P1_1_RV - alpha1*P1_2_RV);

f_alpha0o = (P2_0_RV + alpha1*P2_1_RV + 0.5*alpha1*alpha1*P2_2_RV);
f_alpha0i = (P2_1_RV + alpha1*P2_2_RV);

f_alpha2o = (P2_0_RV - alpha2*P2_1_RV + 0.5*(alpha2*alpha2)*P2_2_RV);
f_alpha2i = (P2_1_RV - alpha2*P2_2_RV);

f_alpha3o = (P3_0_RV + alpha3*(s3*s3*P3_0_RV + 2.0*s3*P3_1_RV + P3_2_RV));
f_alpha3i = (P3_1_RV + alpha3*(s3*s3*P3_1_RV + 2.0*s3*P3_2_RV));

dSL_RV = (sigmaf_RV - sigmapas_RV - sigmaact_RV)/eta;

P0_RV = 1.0 - N_RV - P1_0_RV - P2_0_RV - P3_0_RV; 
dXdT(32) = kf*P0_RV*U_NR_RV*N_overlap_RV - k_f_RV*P1_0_RV - kw_RV*f_alpha1o + k_w*f_alpha0o; 
dXdT(33) = dSL_RV*P1_0_RV - k_f_RV*P1_1_RV - kw_RV*f_alpha1i + k_w*f_alpha0i;
dXdT(34) = 2*dSL_RV*P1_1_RV - k_f_RV*P1_2_RV - kw_RV*P1_2_RV + k_w*P2_2_RV;

dXdT(35) = -k_w*f_alpha0o - kp*f_alpha2o + k_p_RV*P3_0_RV + kw_RV*f_alpha1o;
dXdT(36) = dSL_RV*P2_0_RV - kw_RV*f_alpha0i - kp*f_alpha2i + k_p_RV*P3_1_RV + kw_RV*f_alpha1i;
dXdT(37) = 2*dSL_RV*P2_1_RV - k_w*P2_2_RV       - kp*P2_2_RV + k_p_RV*P3_2_RV + kw_RV*P1_2_RV;

dXdT(38) = +kp*f_alpha2o - k_p_RV*P3_0_RV - kg_RV*f_alpha3o;
dXdT(39) = dSL_RV*P3_0_RV + kp*f_alpha2i - k_p_RV*P3_1_RV - kg_RV*f_alpha3i;
dXdT(40) = 2*dSL_RV*P3_1_RV + kp*P2_2_RV       - k_p_RV*P3_2_RV - kg_RV*P3_2_RV;

U_SR_RV = 1.0 - U_NR_RV;
Jon = k_on*Ca_i*N_RV*(1 + k_coop*(1 - N_RV));
Joff = k_off*P0_RV*(1 + k_coop*N_RV);
dXdT(41) = - Jon + Joff; 
dXdT(42) = km*(1 + krecruit * sigmaact_RV) * U_SR_RV - k_m * U_NR_RV ; 

% Sarcomere length
dXdT(5) = dSL_LV;
dXdT(6) = dSL_SEP;
dXdT(7) = dSL_RV;

% Lumped parameter circulation model variables
dXdT(8) = QIN_LV - QOUT_LV; % V_LV
dXdT(9) = QIN_RV - QOUT_RV; % V_RV
dXdT(43) = (P_SA - P_SV)/R_SA - QIN_RV;  % V_SV
dXdT(44) = (P_PA - P_PV)/R_PA - QIN_LV;  % V_PV
dXdT(45) = Q_Ao - (P_SA - P_SV)/R_SA; % V_SA 
dXdT(46) = QOUT_RV - (P_PA - P_PV)/R_PA; % V_PA 
dXdT(47) = QOUT_LV - Q_Ao; % V_Ao

dXdT = dXdT(:);

end