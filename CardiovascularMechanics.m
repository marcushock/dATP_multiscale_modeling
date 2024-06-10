%%%% Main code to run ventricular model
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


function [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F] = CardiovascularMechanics(Ca, HF, dATP_percent, kf_i, k_f_i, kw_i, k_w_i, kp_i, k_p_i, kg_i, k_recruit_i, k_on_i, k_off_i, k_coop_i, a_i, b_i, c_i)
%% Flags
% Calcium
Ca_flag = Ca; % 0 = ATP, 1 = dATP

% Protocol
HF_flag = HF; % 0 = healthy, 1 = HF

% Percent dATP
dATP = dATP_percent; 
ATP = 100-dATP_percent;

%% Parameters
% Heart
BW = 614; % Body weight (g)
LVW = 1.1127e+03; % LV weight (mg)
RVW = 335.8125; % RV weight (mg)
LW = 2.3414e+03; % Lung weight (mg)
HR = 344.75; % Heart rate (beats/min)

C_Ao = 0.0015; % Proximal aortic compliance (mL/mmHg)
C_SA = 0.0077157; % Systemic arterial compliance (mL/mmHg)
C_SV = 2.5; %2.5; % Systemic venous compliance (mL/mmHg) 
C_PV = 0.25; % Pulmonary venous compliance (mL/mmHg)
C_PA = 0.013778; % Pulmonary arterial compliance (mL/mmHg)
R_Ao = 2.5; % Aortic resistance (mmHg*sec/mL)
R_SA = 56.4565; % Systemic vasculature resistance (mmHg*sec/mL)
R_PA = 7.6986; % Pulmonary arterial resistance (mmHg*sec/mL)
R_SV = 0.25; % Systemic venous resistance (mmHg*sec/mL)
R_PV = 0.25; % Pulmonary venous resistance (mmHg*sec/mL)
R_vlv = 0.05; % Valve resistance (mmHg*sec/mL)
R_AV = R_vlv; % Resistance across aortic valve (mmHg*sec/mL)
R_tAo = 0.5; % Transmural aortic resistance (mmHg*sec/mL)
R_tSA = 4; % Transmural systemic arterial resistance (mmHg*sec/mL)

% Geometry
Amref_LV = 2.0776*1.28; % LV midwall reference surface area (cm^2)
Amref_SEP = 1.2475*1.28; % SEP midwall reference surface area (cm^2)
Amref_RV = 3.1751*1.28; % RV midwall reference surface area (cm^2)

Vw_LV = 0.7065; % LV wall volume (mL)
Vw_SEP = 0.35335; % SEP wall volume (mL)
Vw_RV = 0.31985; % RV wall volume (mL)

% Energetics
TAN = 0.0076; % Total adenine nucleotide pool (M/L cell)
CRtot = 0.0303; % Total creatine pool (M/L cell)
Ox_capacity = 1; % Oxidative capacity
Ox_capacity_sham = 1; % Oxidative capacity for mean sham
Ox_capacity_TAC = 0.7482; %Oxidative capacity for mean TAC
TAN_sham = 0.0076; % Total adenine nucleotide pool for mean sham (M/L cell)
CRtot_sham = 0.0303; % Total creatine pool for mean sham (M/L cell)
TAN_TAC = 0.0070; % Total adenine nucleotide pool for mean TAC (M/L cell)
CRtot_TAC = 0.0230; % Total creatine pool for mean TAC (M/L cell)

Ao = 10.260e-3; % (M/L cell)
Co = 43.007e-3; % (M/L cell)
Po = 35.446e-3; % (M/L cell)

TEP = Po - (0.283e-3)*(Ao-TAN)/(0.082e-3); % Total exchangeable phosphate (M/L cell)
TEP_sham = Po - (0.283e-3)*(Ao-TAN_sham)/(0.082e-3); % Total exchangeable phosphate for mean sham (M/L cell)
TEP_TAC = Po - (0.283e-3)*(Ao-TAN_TAC)/(0.082e-3); % Total exchangeable phosphate pool (failure), mol/L cell

tune_ATPase_LV = 0.0019; % ATP hydrolysis rate (M/s/L cytosol)
rate_of_XB_turnover_mean_sham = 5.0147; % Crossbridge turnover rate (1/s)

if HF_flag == 0 % Healthy
    TAN = TAN_sham;
    CRtot = CRtot_sham;
    TEP = TEP_sham;
    Ox_capacity = Ox_capacity_sham;
else % HF
    TAN = TAN_TAC;
    CRtot = CRtot_TAC;
    TEP = TEP_TAC;
    Ox_capacity = Ox_capacity_TAC;
end

% Crossbridge 
Lsref = 1.9; % Resting sarcomere length (um)
L_rest_pas = 1.51; % Length at which passive force = 0 (um)
gamma = 8; % For calculating passive force

K_Pi = 4.00;
K_T = 0.4897; 
K_D = 0.194;

alpha1 = 10.0; % Stretch sensing parameter for kw and k_w (1/um)
alpha2 = 9.1; % Stretch sensing parameter for kp and k_p (1/um)
alpha3 = 0.1*59.3; % Stretch sensing parameter for k_3 (1/um)
s3 = 9.9e-3;  % Strain at which kg is minimum (um)

Kse = 10000; % Series stiffness (mmHg/um)
k_passive_LV = 58000; % Passive stiffness (mmHg/um)
k_passive_SEP = 58000; 
k_passive_RV = 58000;
eta = 0.0001;  % Viscosity (mmHg*s/um)

L_thick = 1.67; % Length of thick filament (um)
L_hbare = 0.10; % Length of bare region of thick filament (um)
L_thin  = 1.20; % Length of thin filament (um)
deltaR  = 0.010; % um

kstiff1 = 1.4*1.7351e+03*7.5; % Stiffness constant due to myosin actin interaction (kPa/um)
kstiff2 = 1.4*45545*7.5; % Stiffness constant due to working stroke of XBs (kPa/um)

kf = 250*(ATP/100) + kf_i*(dATP/100); % P to A1 transition rate (1/s)
k_f = 304.6708*(ATP/100) + k_f_i*(dATP/100); % A1 to P transition rate (1/s)
kw = 112.3727*(ATP/100) + kw_i*(dATP/100); % A1 to A2 transition rate (1/s)
km1 = 21.296*(ATP/100) + k_w_i*(dATP/100); % A2 to A1 transition rate (1/s)
kp = 811.72*(ATP/100) + kp_i*(dATP/100); % A2 to A3 transition rate (1/s)
km2 = 43.25*(ATP/100) + k_p_i*(dATP/100); % A3 to A2 transition rate (1/s)
kg = 144.5586*(ATP/100) + kg_i*(dATP/100); % A3 to P transition rate (1/s)
K_coop = k_coop_i; % Cooperative constant
k_on = k_on_i;  % N to P transition rate (1/Ms)
k_off = k_off_i;  % P to N transition rate (1/s)
km = 15.4691375; % OFF to ON transition rate (N^-1 m^-1)
k_recruit = 0.2069*(ATP/100) + k_recruit_i*(dATP/100); % Force dependent constant (1/Nm^s)
k_m = 50.032; % ON to OFF transition rate (1/s)

% Calcium
% para_fitted_Ca = [2	3	4	5	6	7	8	9	10;
% 0.0838	0.1306	0.1802	0.2557	0.3099	0.3613	0.408	0.4539	0.4602;
% 0.7513	0.8756	1.0274	1.4988	1.6107	1.6741	1.7902	2.1398	1.9832;
% 2.237	2.0486	1.948	1.852	1.6932	1.6773	1.5988	1.4952	1.4524;
% 0.1865	0.1815	0.1709	0.1693	0.161	0.1661	0.1425	0.1458	0.1222];
% freq_all = para_fitted_Ca(1,:);
% A_HR_pchip = pchip(freq_all,para_fitted_Ca(2,:));
% A_HR = ppval(A_HR_pchip,HR/60);
% B_HR_pchip = pchip(freq_all,para_fitted_Ca(3,:));
% B_HR = ppval(B_HR_pchip,HR/60);
% C_HR_pchip = pchip(freq_all,para_fitted_Ca(4,:));
% C_HR = ppval(C_HR_pchip,HR/60);
% Ca0_HR_pchip = pchip(freq_all,para_fitted_Ca(5,:));
% Ca0_HR = ppval(Ca0_HR_pchip,HR/60);

A_HR = a_i;
B_HR = b_i;
C_HR = c_i;
Ca0_HR = 0.1622;

%% Initial conditions
% Heart
V_LV = 0.6167+0.1; % LV volume (mL)
V_RV = 0.6167+0.1; % RV volume (mL)
V_SA = 3.3043; % Systemic arteries volume (mL)
V_SV = 5.2868; % Systemic veins volume (mL)
V_PA = 0.5507; % Pulmonary arteries volume (mL)
V_PV = 1.1014; % Pulmonary veins volume (mL)
V_Ao = 1.1014; % Aorta volume (mL)

% Geometry
xm_LV = -0.60; % Maximal axial distance from LV midwall surface to origin (cm)
xm_SEP = 0.40; % Maximal axial distance from SEP midwall surface to origin (cm)
xm_RV = 1.0; % Maximal axial distance from RV midwall surface to origin (cm)
ym = 0.50; % Radius of midwall junction circle (cm)

% Energetics
energtics_output  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV); % Get metabolite concentrations

MgATP_cyto = energtics_output(1); % Cytosolic MgATP concentration (M/(L cytosol water)^-1)
MgADP_cyto = energtics_output(2); % Cytosolic MgADP concentration (M/(L cytosol water)^-1)
fPi_cytoplasm = energtics_output(3); % Cytosolic unchelated Pi concentration (M/(L cytosol water)^-1)
MVO2_tissue = energtics_output(5); % Oxygen consumption rate (uM/min^-1(g tissue)^-1)
dGrATPase = energtics_output(6); % ATP hydrolysis free energy, kJ/mol

PCrATP =  energtics_output(7); % Creatine phoosphate ATP ratio, unitless
ATP_cyto = energtics_output(8); % Cytosolic total ATP concentration (M/(L cytosol water)^-1)
ADP_cyto = energtics_output(9); % Cytosolic total ADP concentration (M/(L cytosol water)^-1)
Pi_cyto = energtics_output(10)*1000; % Cytosolic total Pi concentration (M/(L cytosol water)^-1)

ATP_store(1) = MgATP_cyto;
ADP_store(1) = MgADP_cyto*1000;
Pi_store(1) = Pi_cyto;
MVO2_tissue_store(1) = MVO2_tissue;
PCrATP_store(1) = PCrATP;

% Crossbridge
SL_LV = 2.2; % LV sarcomere length (um)
SL_SEP = 2.2; % SEP sarcomere length (um)
SL_RV = 2.2; % RV sarcomere length (um)

P1_0_LV = 0; % 0th moment state A1, LV
P1_1_LV = 0; % 1st moment state A1, LV
P1_2_LV = 0; % 2nd moment state A1, LV
P2_0_LV = 0; % 0th moment state A2, LV
P2_1_LV = 0; % 1st moment state A2, LV
P2_2_LV = 0; % 2nd moment state A2, LV
P3_0_LV = 0; % 0th moment state A3, LV
P3_1_LV = 0; % 1st moment state A3, LV
P3_2_LV = 0; % 2nd moment state A3, LV
N_LV = 1; % Non permissible state, LV
U_NR_LV = 0; % ON state, LV
P1_0_SEP = 0; % 0th moment state A1, SEP
P1_1_SEP = 0; % 1st moment state A1, SEP
P1_2_SEP = 0; % 2nd moment state A1, SEP
P2_0_SEP = 0; % 0th moment state A2, SEP
P2_1_SEP= 0; % 1st moment state A2, SEP
P2_2_SEP = 0; % 2nd moment state A2, SEP
P3_0_SEP = 0; % 0th moment state A3, SEP
P3_1_SEP = 0; % 1st moment state A3, SEP
P3_2_SEP = 0; % 2nd moment state A3, SEP
N_SEP = 1; % Non permissible state, SEP
U_NR_SEP = 0; % ON state, SEP
P1_0_RV = 0; % 0th moment state A1, RV
P1_1_RV = 0; % 1st moment state A1, RV
P1_2_RV = 0; % 2nd moment state A1, RV
P2_0_RV = 0; % 0th moment state A2, RV
P2_1_RV = 0; % 1st moment state A2, RV
P2_2_RV = 0; % 2nd moment state A2, RV
P3_0_RV = 0; % 0th moment state A3, RV
P3_1_RV = 0; % 1st moment state A3, RV
P3_2_RV = 0; % 2nd moment state A3, RV
N_RV = 1; % Non permissible state, RV
U_NR_RV = 0; % ON state, RV

init = [xm_LV ,xm_SEP ,xm_RV ,ym , SL_LV, SL_SEP, SL_RV, V_LV, V_RV, ...
       P1_0_LV, P1_1_LV, P1_2_LV ,P2_0_LV, P2_1_LV, P2_2_LV, P3_0_LV, P3_1_LV, P3_2_LV, N_LV, U_NR_LV,...
       P1_0_SEP,P1_1_SEP,P1_2_SEP,P2_0_SEP,P2_1_SEP,P2_2_SEP,P3_0_SEP,P3_1_SEP,P3_2_SEP,N_SEP,U_NR_SEP,...
       P1_0_RV, P1_1_RV, P1_2_RV, P2_0_RV, P2_1_RV, P2_2_RV, P3_0_RV, P3_1_RV, P3_2_RV, N_RV, U_NR_RV,...
       V_SV, V_PV ,V_SA ,V_PA, V_Ao]';

% Get consistent initial conditions for solving DAE
opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
x_triseg = fsolve(@TrisegEquations,init(1:4), opts, Vw_LV, Vw_SEP, Vw_RV, SL_LV, SL_SEP, SL_RV, V_LV, V_RV, Amref_LV, Amref_SEP, Amref_RV);
init(1:4) = x_triseg;
   
%% Run initially to steady state without coupled energetics model
stim_period = 1/(HR/60);
M = speye(47);
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 

input = [stim_period, C_Ao, C_SA, C_SV, C_PV, C_PA, R_Ao, R_SA, R_PA, R_SV, R_PV, R_vlv, R_AV, ...
    R_tAo, R_tSA, Amref_LV, Amref_SEP, Amref_RV, Vw_LV, Vw_SEP, Vw_RV, MgATP_cyto, ...
    MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, ...
    Pi_cyto, Lsref, L_rest_pas, gamma, K_Pi, K_T, K_D, alpha1, alpha2, alpha3, s3, Kse, ...
    k_passive_LV, k_passive_SEP, k_passive_RV, eta, L_thick, L_hbare, L_thin, deltaR, kstiff1, kstiff2, kf, k_f, kw, ...
    km1, kp, km2, kg, K_coop, k_on, k_off, km, k_recruit, k_m, A_HR, B_HR, C_HR, Ca0_HR, dATP, Ca_flag];
options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-7,'AbsTol',1e-7,'MaxStep',stim_period/200);

[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 120*stim_period], init, options, input);
t_store{1} = t;
Y_store{1} = Y;

% Calculate
% State variables
xm_LV  = Y(:,1); % Maximal axial distance from LV midwall surface to origin (cm)
xm_SEP = Y(:,2); % Maximal axial distance from SEP midwall surface to origin (cm)
xm_RV  = Y(:,3); % Maximal axial distance from RV midwall surface to origin (cm)
ym     = Y(:,4); % Radius of midwall junction circle (cm)
SL_LV  = Y(:,5); % LV sarcomere length (um)
SL_SEP = Y(:,6); % SEP sarcomere length (um)
SL_RV  = Y(:,7); % RV sarcomere length (um)
V_LV   = Y(:,8); % LV volume (mL)
V_RV   = Y(:,9); % RV volume (mL)

% Store
xm_LV_store{1} = xm_LV;
xm_SEP_store{1} = xm_SEP;
xm_RV_store{1} = xm_RV;
ym_store{1} = ym;
SL_LV_store{1} = SL_LV; 
SL_SEP_store{1} = SL_SEP; 
SL_RV_store{1} = SL_RV; 
V_LV_store{1} = V_LV;
V_RV_store{1} = V_RV;

P1_0_LV = Y(:,10); % 0th moment state A1, LV
P1_1_LV = Y(:,11); % 1st moment state A1, LV
P1_2_LV = Y(:,12); % 2nd moment state A1, LV
P2_0_LV = Y(:,13); % 0th moment state A2, LV
P2_1_LV = Y(:,14); % 1st moment state A2, LV
P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
N_LV    = Y(:,19); % Non permissible state, LV
U_NR_LV = Y(:,20); % ON, LV

% Store
A1_LV_store{1} = P1_0_LV;
A2_LV_store{1} = P2_0_LV;
A3_LV_store{1} = P3_0_LV;
N_LV_store{1} = N_LV;

P1_0_SEP = Y(:,21); % 0th moment state A1, SEP
P1_1_SEP = Y(:,22); % 1st moment state A1, SEP
P1_2_SEP = Y(:,23); % 2nd moment state A1, SEP
P2_0_SEP = Y(:,24); % 0th moment state A2, SEP
P2_1_SEP = Y(:,25); % 1st moment state A2, SEP
P2_2_SEP = Y(:,26); % 2nd moment state A2, SEP
P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP
N_SEP    = Y(:,30); % Non permissible state, SEP
U_NR_SEP = Y(:,31); % ON state, SEP

% Store
A1_SEP_store{1} = P1_0_SEP;
A2_SEP_store{1} = P2_0_SEP;
A3_SEP_store{1} = P3_0_SEP;
N_SEP_store{1} = N_SEP;

P1_0_RV = Y(:,32); % 0th moment state A1, RV
P1_1_RV = Y(:,33); % 1st moment state A1, RV
P1_2_RV = Y(:,34); % 2nd moment state A1, RV
P2_0_RV = Y(:,35); % 0th moment state A2, RV
P2_1_RV = Y(:,36); % 1st moment state A2, RV
P2_2_RV = Y(:,37); % 2nd moment state A2, RV
P3_0_RV = Y(:,38); % 0th moment state A3, RV
P3_1_RV = Y(:,39); % 1st moment state A3, RV
P3_2_RV = Y(:,40); % 2nd moment state A3, RV
N_RV    = Y(:,41); % Non permissible state, RV
U_NR_RV = Y(:,42); % ON state, RV 

% Store
A1_RV_store{1} = P1_0_RV;
A2_RV_store{1} = P2_0_RV;
A3_RV_store{1} = P3_0_RV;
N_RV_store{1} = N_RV;

V_SV   = Y(:,43); % Systemic veins volume (mL)
V_PV   = Y(:,44); % Pulmonary veins volume (mL)
V_SA   = Y(:,45); % Systemic arteries volume (mL)
V_PA   = Y(:,46); % Pulmonary arteries volume (mL)
V_Ao   = Y(:,47); % Aorta volume (mL)
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;

% Store
V_SV_store{1} = V_SV;
V_PV_store{1} = V_PV;
V_SA_store{1} = V_SA;
V_PA_store{1} = V_PA;
V_Ao_store{1} = V_Ao;

% Geometry
Vm_LV  = (pi/6).*xm_LV.*(xm_LV.^2 + 3.*ym.^2); % LV midwall volume (mL)
Vm_SEP = (pi/6).*xm_SEP.*(xm_SEP.^2 + 3.*ym.^2); % SEP midwall volume (mL)
Vm_RV  = (pi/6).*xm_RV.*(xm_RV.^2 + 3.*ym.^2); % RV midwall volume (mL)
Am_LV  = pi.*(xm_LV.^2 + ym.^2); % LV midwall surface area (cm^2)
Am_SEP = pi.*(xm_SEP.^2 + ym.^2); % SEP midwall surface area (cm^2)
Am_RV  = pi.*(xm_RV.^2 + ym.^2); % RV midwall surface area (cm^2)
Cm_LV  = 2.*xm_LV./(xm_LV.^2 + ym.^2); % LV midwall surface curvature
Cm_SEP = 2.*xm_SEP./(xm_SEP.^2 + ym.^2); % SEP midwall surface curvature
Cm_RV  = 2.*xm_RV./(xm_RV.^2 + ym.^2); % RV midwall surface curvature
z_LV   = 3.*Cm_LV.*Vw_LV./(2.*Am_LV); 
z_SEP  = 3.*Cm_SEP.*Vw_SEP./(2.*Am_SEP);
z_RV   = 3.*Cm_RV.*Vw_RV./(2.*Am_RV);

% Store
Vm_LV_store{1} = Vm_LV;
Vm_SEP_store{1} = Vm_SEP;
Vm_RV_store{1} = Vm_RV;
Am_LV_store{1} = Am_LV;
Am_SEP_store{1} = Am_SEP;
Am_RV_store{1} = Am_RV;
Cm_LV_store{1} = Cm_LV;
Cm_SEP_store{1} = Cm_SEP;
Cm_RV_store{1} = Cm_RV;
z_LV_store{1} = z_LV;
z_SEP_store{1} = z_SEP;
z_RV_store{1} = z_RV;

% Fiber strain
epsf_LV = 0.5.*log(Am_LV./Amref_LV) - 0.083333.*z_LV.^2 - 0.019.*z_LV.^4;
epsf_SEP = 0.5.*log(Am_SEP./Amref_SEP) - 0.083333.*z_SEP.^2 - 0.019.*z_SEP.^4;
epsf_RV = 0.5.*log(Am_RV./Amref_RV) - 0.083333.*z_RV.^2 - 0.019.*z_RV.^4;

% Sarcomere length
SLo_LV = Lsref.*exp(epsf_LV); 
SLo_SEP = Lsref.*exp(epsf_SEP); 
SLo_RV = Lsref.*exp(epsf_RV);

% Store
epsf_LV_store{1} = epsf_LV;
epsf_SEP_store{1} = epsf_SEP;
epsf_RV_store{1} = epsf_RV;
SLo_LV_store{1} = SLo_LV;
SLo_SEP_store{1} = SLo_SEP;
SLo_RV_store{1} = SLo_RV;

% Passive force
sigmapas_LV  = k_passive_LV.*(SLo_LV - L_rest_pas).^(gamma);
sigmapas_SEP  = k_passive_SEP.*(SLo_SEP - L_rest_pas).^(gamma);
sigmapas_RV  = k_passive_RV.*(SLo_RV - L_rest_pas).^(gamma);

% Store
sigmapas_LV_store{1} = sigmapas_LV;
sigmapas_SEP_store{1} = sigmapas_SEP;
sigmapas_RV_store{1} = sigmapas_RV;

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

sigmaact_LV  = N_overlap_LV.*(kstiff2.*deltaR.*(P3_0_LV) + kstiff1.*(P2_1_LV + P3_1_LV)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_SEP = N_overlap_SEP.*(kstiff2.*deltaR.*(P3_0_SEP) + kstiff1.*(P2_1_SEP + P3_1_SEP)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_RV  = N_overlap_RV.*(kstiff2.*deltaR.*(P3_0_RV) + kstiff1.*(P2_1_RV + P3_1_RV)); % Includes force due to XB ratcheting and stretching of XBs

% Store
sigmaact_LV_store{1} = sigmaact_LV;
sigmaact_SEP_store{1} = sigmaact_SEP;
sigmaact_RV_store{1} = sigmaact_RV;

% Total force
sigmaf_LV = -Kse.*(SL_LV - SLo_LV);
sigmaf_SEP = -Kse.*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse.*(SL_RV - SLo_RV);

% Store
sigmaf_LV_store{1} = sigmaf_LV;
sigmaf_SEP_store{1} = sigmaf_SEP;
sigmaf_RV_store{1} = sigmaf_RV;

% Tension
Tm_LV  = (Vw_LV.*sigmaf_LV./(2.*Am_LV)).*(1 + (z_LV.^2)/3 + (z_LV.^4)/5); % LV midwall tension
Tm_SEP = (Vw_SEP.*sigmaf_SEP./(2.*Am_SEP)).*(1 + (z_SEP.^2)/3 + (z_SEP.^4)/5); % SEP midwall tension
Tm_RV  = (Vw_RV.*sigmaf_RV./(2.*Am_RV)).*(1 + (z_RV.^2)/3 + (z_RV.^4)/5); % RV midwall tension

Tx_LV  = Tm_LV.*2.*xm_LV.*ym./(xm_LV.^2 + ym.^2); % LV axial tension
Tx_SEP = Tm_SEP.*2.*xm_SEP.*ym./(xm_SEP.^2 + ym.^2); % SEP axial tension
Tx_RV  = Tm_RV.*2.*xm_RV.*ym./(xm_RV.^2 + ym.^2); % RV axial tension

Ty_LV  = Tm_LV.*(-xm_LV.^2 + ym.^2)./(xm_LV.^2 + ym.^2); % LV radial tension
Ty_SEP = Tm_SEP.*(-xm_SEP.^2 + ym.^2)./(xm_SEP.^2 + ym.^2); % SEP radial tension
Ty_RV  = Tm_RV.*(-xm_RV.^2 + ym.^2)./(xm_RV.^2 + ym.^2); % RV radial tension

% Store
Tm_LV_store{1} = Tm_LV;
Tm_SEP_store{1} = Tm_SEP;
Tm_RV_store{1} = Tm_RV;

Tx_LV_store{1} = Tx_LV;
Tx_SEP_store{1} = Tx_SEP;
Tx_RV_store{1} = Tx_RV;

Ty_LV_store{1} = Ty_LV;
Ty_SEP_store{1} = Ty_SEP;
Ty_RV_store{1} = Ty_RV;

% Ventricular pressures
ptrans1 = 2.*Tx_LV./ym;
ptrans3 = 2.*Tx_RV./ym;
P_LV = -ptrans1; % LV transmural pressure (mmHg)
P_RV = ptrans3; % RV transmural pressure (mmHg)

% Store
P_LV_store{1} = P_LV;
P_RV_store{1} = P_RV;

% Pulmonary pressures
P_PV = V_PV./C_PV; % Pulmonary venous pressure (mmHg)
P_SV = V_SV./C_SV; % Systemic venous pressure (mmHg)
P_PA = V_PA./C_PA; % Pulmonary arterial pressure (mmHg)
P_SA = V_SA./C_SA; % Systemic arterial pressure (mmHg)

% Store
P_PV_store{1} = P_PV;
P_SV_store{1} = P_SV;
P_PA_store{1} = P_PA;
P_SA_store{1} = P_SA;

% Lumped circulatory model
% Ao valve closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));

% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + P_Ao_closed.*(P_LV<=P_Ao_open);
P_Ao_store{1} = P_Ao;

SV_LV_sim = max(1e3*V_LV) - min(1e3*V_LV); % Stroke volume
EF_LV_sim = SV_LV_sim/max(1e3*V_LV) * 100; % Ejection fraction

edLV_sim =  max(1e3*V_LV);
esLV_sim =  min(1e3*V_LV);
edRV_sim =  max(1e3*V_RV);
esRV_sim =  min(1e3*V_RV);

% Defining the metabolite dependent coeficient, rapid equilibrium of the
% cross bridge sub-states
g2_LV = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
kg_LV = kg*g2_LV;
g2_SEP = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
kg_SEP = kg*g2_SEP;
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 

% Detachment rates
ti = 0:0.00001:stim_period;
MAP = mean(interp1(t,P_Ao,ti)); % Mean arterial pressure (mmHg)

r_LV  = interp1(t,kg_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,kg_SEP*f_alpha3o_SEP,ti);

% Crossbridge turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W);

ATP_ase_mechanics_Averge_LV_SEP = (1.3/rate_of_XB_turnover_mean_sham)*rate_of_XB_turnover_ave; %  1.31 Kstiff - ATP hydrolized (mM/s/L cell) per X-bridge turnover rate in LV
tune_ATPase_LV =  ATP_ase_mechanics_Averge_LV_SEP * (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate (M/s/L cytosol)
ATPase_store(1) = tune_ATPase_LV;


%% Run coupled cardiovascular mechanics/energetics model
p = 2;
beat = 1;
jspan = 122:420; 

for j = jspan % Run for 300 beats
tspan_beat = [(stim_period)*(j-1) (stim_period)*j];

% Run energetics model
if rem(j,3) == 0 % Run every 3 beats (for model stability)
energtics_output  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV); % Get metabolite concentrations

MgATP_cyto = energtics_output(1); % Cytosolic MgATP concentration (M/(L cytosol water)^-1)
MgADP_cyto = energtics_output(2); % Cytosolic MgADP concentration (M/(L cytosol water)^-1)
fPi_cytoplasm = energtics_output(3); % Cytosolic unchelated Pi concentration (M/(L cytosol water)^-1)
MVO2_tissue = energtics_output(5); % Oxygen consumption rate (uM/min^-1(g tissue)^-1)
dGrATPase = energtics_output(6); % ATP hydrolysis free energy, kJ/mol

PCrATP =  energtics_output(7); % Creatine phoosphate ATP ratio, unitless
ATP_cyto = energtics_output(8); % Cytosolic total ATP concentration (M/(L cytosol water)^-1)
ADP_cyto = energtics_output(9); % Cytosolic total ADP concentration (M/(L cytosol water)^-1)
Pi_cyto = energtics_output(10)*1000; % Cytosolic total Pi concentration (M/(L cytosol water)^-1)
end

% Store
ATP_store(p) = MgATP_cyto;
ADP_store(p) = MgADP_cyto*1000;
Pi_store(p) = Pi_cyto;
MVO2_tissue_store(p) = MVO2_tissue;
PCrATP_store(p) = PCrATP;

% Run mechanics model
init = Y(end,:);
input = [stim_period, C_Ao, C_SA, C_SV, C_PV, C_PA, R_Ao, R_SA, R_PA, R_SV, R_PV, R_vlv, R_AV, ...
    R_tAo, R_tSA, Amref_LV, Amref_SEP, Amref_RV, Vw_LV, Vw_SEP, Vw_RV, MgATP_cyto, ...
    MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, ...
    Pi_cyto, Lsref, L_rest_pas, gamma, K_Pi, K_T, K_D, alpha1, alpha2, alpha3, s3, Kse, ...
    k_passive_LV, k_passive_SEP, k_passive_RV, eta, L_thick, L_hbare, L_thin, deltaR, kstiff1, kstiff2, kf, k_f, kw, ...
    km1, kp, km2, kg, K_coop, k_on, k_off, km, k_recruit, k_m, A_HR, B_HR, C_HR, Ca0_HR, dATP, Ca_flag];
options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-7,'AbsTol',1e-7,'MaxStep',stim_period/200);

[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,tspan_beat, init, options, input);
t_store{p} = t;
Y_store{p} = Y;

% Calculate
% State variables
xm_LV  = Y(:,1); % Maximal axial distance from LV midwall surface to origin (cm)
xm_SEP = Y(:,2); % Maximal axial distance from SEP midwall surface to origin (cm)
xm_RV  = Y(:,3); % Maximal axial distance from RV midwall surface to origin (cm)
ym     = Y(:,4); % Radius of midwall junction circle (cm)
SL_LV  = Y(:,5); % LV sarcomere length (um)
SL_SEP = Y(:,6); % SEP sarcomere length (um)
SL_RV  = Y(:,7); % RV sarcomere length (um)
V_LV   = Y(:,8); % LV volume (mL)
V_RV   = Y(:,9); % RV volume (mL)

% Store
xm_LV_store{p} = xm_LV;
xm_SEP_store{p} = xm_SEP;
xm_RV_store{p} = xm_RV;
ym_store{p} = ym;
SL_LV_store{p} = SL_LV; 
SL_SEP_store{p} = SL_SEP; 
SL_RV_store{p} = SL_RV; 
V_LV_store{p} = V_LV;
V_RV_store{p} = V_RV;

P1_0_LV = Y(:,10); % 0th moment state A1, LV
P1_1_LV = Y(:,11); % 1st moment state A1, LV
P1_2_LV = Y(:,12); % 2nd moment state A1, LV
P2_0_LV = Y(:,13); % 0th moment state A2, LV
P2_1_LV = Y(:,14); % 1st moment state A2, LV
P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
N_LV    = Y(:,19); % Non permissible state, LV
U_NR_LV = Y(:,20); % ON state, LV

% Store
A1_LV_store{p} = P1_0_LV;
A2_LV_store{p} = P2_0_LV;
A3_LV_store{p} = P3_0_LV;
N_LV_store{p} = N_LV;

P1_0_SEP = Y(:,21); % 0th moment state A1, SEP
P1_1_SEP = Y(:,22); % 1st moment state A1, SEP
P1_2_SEP = Y(:,23); % 2nd moment state A1, SEP
P2_0_SEP = Y(:,24); % 0th moment state A2, SEP
P2_1_SEP = Y(:,25); % 1st moment state A2, SEP
P2_2_SEP = Y(:,26); % 2nd moment state A2, SEP
P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP
N_SEP    = Y(:,30); % Non permissible state, SEP
U_NR_SEP = Y(:,31); % ON state, SEP

% Store
A1_SEP_store{p} = P1_0_SEP;
A2_SEP_store{p} = P2_0_SEP;
A3_SEP_store{p} = P3_0_SEP;
N_SEP_store{p} = N_SEP;

P1_0_RV = Y(:,32); % 0th moment state A1, RV
P1_1_RV = Y(:,33); % 1st moment state A1, RV
P1_2_RV = Y(:,34); % 2nd moment state A1, RV
P2_0_RV = Y(:,35); % 0th moment state A2, RV
P2_1_RV = Y(:,36); % 1st moment state A2, RV
P2_2_RV = Y(:,37); % 2nd moment state A2, RV
P3_0_RV = Y(:,38); % 0th moment state A3, RV
P3_1_RV = Y(:,39); % 1st moment state A3, RV
P3_2_RV = Y(:,40); % 2nd moment state A3, RV
N_RV    = Y(:,41); % Non permissible state, RV
U_NR_RV = Y(:,42); % ON state, RV

% Store
A1_RV_store{p} = P1_0_RV;
A2_RV_store{p} = P2_0_RV;
A3_RV_store{p} = P3_0_RV;
N_RV_store{p} = N_RV;

V_SV   = Y(:,43); % Systemic veins volume (mL)
V_PV   = Y(:,44); % Pulmonary veins volume (mL)
V_SA   = Y(:,45); % Systemic arteries volume (mL)
V_PA   = Y(:,46); % Pulmonary arteries volume (mL)
V_Ao   = Y(:,47); % Aorta volume (mL)
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;

% Store
V_SV_store{p} = V_SV;
V_PV_store{p} = V_PV;
V_SA_store{p} = V_SA;
V_PA_store{p} = V_PA;
V_Ao_store{p} = V_Ao;

% Geometry
Vm_LV  = (pi/6).*xm_LV.*(xm_LV.^2 + 3.*ym.^2); % LV midwall volume (mL)
Vm_SEP = (pi/6).*xm_SEP.*(xm_SEP.^2 + 3.*ym.^2); % SEP midwall volume (mL)
Vm_RV  = (pi/6).*xm_RV.*(xm_RV.^2 + 3.*ym.^2); % RV midwall volume (mL)
Am_LV  = pi.*(xm_LV.^2 + ym.^2); % LV midwall surface area (cm^2)
Am_SEP = pi.*(xm_SEP.^2 + ym.^2); % SEP midwall surface area (cm^2)
Am_RV  = pi.*(xm_RV.^2 + ym.^2); % RV midwall surface area (cm^2)
Cm_LV  = 2.*xm_LV./(xm_LV.^2 + ym.^2); % LV midwall surface curvature
Cm_SEP = 2.*xm_SEP./(xm_SEP.^2 + ym.^2); % LV midwall surface curvature
Cm_RV  = 2.*xm_RV./(xm_RV.^2 + ym.^2); % LV midwall surface curvature
z_LV   = 3.*Cm_LV.*Vw_LV./(2.*Am_LV);
z_SEP  = 3.*Cm_SEP.*Vw_SEP./(2.*Am_SEP);
z_RV   = 3.*Cm_RV.*Vw_RV./(2.*Am_RV);

% Store
Vm_LV_store{p} = Vm_LV;
Vm_SEP_store{p} = Vm_SEP;
Vm_RV_store{p} = Vm_RV;
Am_LV_store{p} = Am_LV;
Am_SEP_store{p} = Am_SEP;
Am_RV_store{p} = Am_RV;
Cm_LV_store{p} = Cm_LV;
Cm_SEP_store{p} = Cm_SEP;
Cm_RV_store{p} = Cm_RV;
z_LV_store{p} = z_LV;
z_SEP_store{p} = z_SEP;
z_RV_store{p} = z_RV;

% Fiber strain
epsf_LV = 0.5.*log(Am_LV./Amref_LV) - 0.083333.*z_LV.^2 - 0.019.*z_LV.^4;
epsf_SEP = 0.5.*log(Am_SEP./Amref_SEP) - 0.083333.*z_SEP.^2 - 0.019.*z_SEP.^4;
epsf_RV = 0.5.*log(Am_RV./Amref_RV) - 0.083333.*z_RV.^2 - 0.019.*z_RV.^4;
SLo_LV = Lsref.*exp(epsf_LV); 
SLo_SEP = Lsref.*exp(epsf_SEP); 
SLo_RV = Lsref.*exp(epsf_RV);

% Store
epsf_LV_store{p} = epsf_LV;
epsf_SEP_store{p} = epsf_SEP;
epsf_RV_store{p} = epsf_RV;
SLo_LV_store{p} = SLo_LV;
SLo_SEP_store{p} = SLo_SEP;
SLo_RV_store{p} = SLo_RV;

% Passive force
sigmapas_LV  = k_passive_LV.*(SLo_LV - L_rest_pas).^(gamma);
sigmapas_SEP  = k_passive_SEP.*(SLo_SEP - L_rest_pas).^(gamma);
sigmapas_RV  = k_passive_RV.*(SLo_RV - L_rest_pas).^(gamma);

% Store
sigmapas_LV_store{p} = sigmapas_LV;
sigmapas_SEP_store{p} = sigmapas_SEP;
sigmapas_RV_store{p} = sigmapas_RV;

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

sigmaact_LV  = N_overlap_LV.*(kstiff2.*deltaR.*(P3_0_LV) + kstiff1.*(P2_1_LV + P3_1_LV)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_SEP = N_overlap_SEP.*(kstiff2.*deltaR.*(P3_0_SEP) + kstiff1.*(P2_1_SEP + P3_1_SEP)); % Includes force due to XB ratcheting and stretching of XBs
sigmaact_RV  = N_overlap_RV.*(kstiff2.*deltaR.*(P3_0_RV) + kstiff1.*(P2_1_RV + P3_1_RV)); % Includes force due to XB ratcheting and stretching of XBs

% Store
sigmaact_LV_store{p} = sigmaact_LV;
sigmaact_SEP_store{p} = sigmaact_SEP;
sigmaact_RV_store{p} = sigmaact_RV;

% Total force
sigmaf_LV = -Kse.*(SL_LV - SLo_LV);
sigmaf_SEP = -Kse.*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse.*(SL_RV - SLo_RV);

% Store
sigmaf_LV_store{p} = sigmaf_LV;
sigmaf_SEP_store{p} = sigmaf_SEP;
sigmaf_RV_store{p} = sigmaf_RV;

% Tension
Tm_LV  = (Vw_LV.*sigmaf_LV./(2.*Am_LV)).*(1 + (z_LV.^2)/3 + (z_LV.^4)/5); % LV midwall tension
Tm_SEP = (Vw_SEP.*sigmaf_SEP./(2.*Am_SEP)).*(1 + (z_SEP.^2)/3 + (z_SEP.^4)/5); % SEP midwall tension
Tm_RV  = (Vw_RV.*sigmaf_RV./(2.*Am_RV)).*(1 + (z_RV.^2)/3 + (z_RV.^4)/5); % RV midwall tension

Tx_LV  = Tm_LV.*2.*xm_LV.*ym./(xm_LV.^2 + ym.^2); % LV axial tension
Tx_SEP = Tm_SEP.*2.*xm_SEP.*ym./(xm_SEP.^2 + ym.^2); % SEP axial tension
Tx_RV  = Tm_RV.*2.*xm_RV.*ym./(xm_RV.^2 + ym.^2); % RV axial tension

Ty_LV  = Tm_LV.*(-xm_LV.^2 + ym.^2)./(xm_LV.^2 + ym.^2); % LV radial tension
Ty_SEP = Tm_SEP.*(-xm_SEP.^2 + ym.^2)./(xm_SEP.^2 + ym.^2); % SEP radial tension
Ty_RV  = Tm_RV.*(-xm_RV.^2 + ym.^2)./(xm_RV.^2 + ym.^2); % RV radial tension

% Store
Tm_LV_store{p} = Tm_LV;
Tm_SEP_store{p} = Tm_SEP;
Tm_RV_store{p} = Tm_RV;

Tx_LV_store{p} = Tx_LV;
Tx_SEP_store{p} = Tx_SEP;
Tx_RV_store{p} = Tx_RV;

Ty_LV_store{p} = Ty_LV;
Ty_SEP_store{p} = Ty_SEP;
Ty_RV_store{p} = Ty_RV;

% Ventricular pressures
ptrans1 = 2.*Tx_LV./ym;
ptrans3 = 2.*Tx_RV./ym;
P_LV = -ptrans1; % LV transmural pressure (mmHg)
P_RV = ptrans3; % RV transmural pressure (mmHg)

% Store
P_LV_store{p} = P_LV;
P_RV_store{p} = P_RV;

% Pulmonary pressures
P_PV = V_PV./C_PV; % Pulmonary venous pressure (mmHg)
P_SV = V_SV./C_SV; % Systemic venous pressure (mmHg)
P_PA = V_PA./C_PA; % Pulmonary arterial pressure (mmHg)
P_SA = V_SA./C_SA; % Systemic arterial pressure (mmHg)

% Store
P_PV_store{p} = P_PV;
P_SV_store{p} = P_SV;
P_PA_store{p} = P_PA;
P_SA_store{p} = P_SA;

% Lumped circulatory model
% Ao valve closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));

% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + P_Ao_closed.*(P_LV<=P_Ao_open);
P_Ao_store{p} = P_Ao;

SV_LV_sim = max(1e3*V_LV) - min(1e3*V_LV); % Stroke volume
EF_LV_sim = SV_LV_sim/max(1e3*V_LV) * 100; % Ejection fraction

edLV_sim =  max(1e3*V_LV);
esLV_sim =  min(1e3*V_LV);
edRV_sim =  max(1e3*V_RV);
esRV_sim =  min(1e3*V_RV);

g2_LV = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
kg_LV = kg*g2_LV;
g2_SEP = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
kg_SEP = kg*g2_SEP;
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 

% Detachment rates
ti = stim_period*(j-1):0.00001:stim_period*j;
MAP = mean(interp1(t,P_Ao,ti)); % Mean arterial pressure (mmHg)

r_LV  = interp1(t,kg_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,kg_SEP*f_alpha3o_SEP,ti);

% Crossbridge turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W);

ATP_ase_mechanics_Averge_LV_SEP = (1.3/rate_of_XB_turnover_mean_sham)*rate_of_XB_turnover_ave; %  1.31 Kstiff - ATP hydrolized (mM/s/L cell) per XB turnover rate in LV
tune_ATPase_LV =  ATP_ase_mechanics_Averge_LV_SEP * (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate (M/s/L cytosol)
ATPase_store(p) = tune_ATPase_LV;

beat = p
p = p+1;

end


%% Finalize output
% Initialize
k = 1;
t_f = [];
xm_LV_f = [];
xm_SEP_f = [];
xm_RV_f = [];
ym_f = [];
SL_LV_f = [];
SL_SEP_f = [];
SL_RV_f = [];
V_LV_f = [];
V_RV_f = [];
A1_LV_f = [];
A2_LV_f = [];
A3_LV_f = [];
N_LV_f = [];
A1_SEP_f = [];
A2_SEP_f = [];
A3_SEP_f = [];
N_SEP_f = [];
A1_RV_f = [];
A2_RV_f = [];
A3_RV_f = [];
N_RV_f = [];
V_SV_f = [];
V_PV_f = [];
V_SA_f = [];
V_PA_f = [];
V_Ao_f = [];
Vm_LV_f = [];
Vm_SEP_f = [];
Vm_RV_f = [];
Am_LV_f = [];
Am_SEP_f = [];
Am_RV_f = [];
Cm_LV_f = [];
Cm_SEP_f = [];
Cm_RV_f = [];
z_LV_f = [];
z_SEP_f = [];
z_RV_f = [];
epsf_LV_f = [];
epsf_SEP_f = [];
epsf_RV_f = [];
SLo_LV_f = [];
SLo_SEP_f = [];
SLo_RV_f = [];
sigmapas_LV_f = [];
sigmapas_SEP_f = [];
sigmapas_RV_f = [];
sigmaact_LV_f = [];
sigmaact_SEP_f = [];
sigmaact_RV_f = [];
sigmaf_LV_f = [];
sigmaf_SEP_f = [];
sigmaf_RV_f = [];
Tm_LV_f = [];
Tm_SEP_f = [];
Tm_RV_f = [];
Tx_LV_f = [];
Tx_SEP_f = [];
Tx_RV_f = [];
Ty_LV_f = [];
Ty_SEP_f = [];
Ty_RV_f = [];
P_LV_f = [];
P_RV_f = [];
P_SV_f = [];
P_PV_f = [];
P_SA_f = [];
P_PA_f = [];
P_Ao_f = [];

while k < p-1
    t_1 = t_store{k};
    t_2 = t_store{k+1};
    t_new = vertcat(t_1(2:end),t_2(2:end));
    t_f = vertcat(t_f,t_new);
    
    xm_LV1 = xm_LV_store{k};
    xm_LV2 = xm_LV_store{k+1};
    xm_LV_new = vertcat(xm_LV1(2:end),xm_LV2(2:end));
    xm_LV_f  =  vertcat(xm_LV_f,xm_LV_new);
    
    xm_SEP1 = xm_SEP_store{k};
    xm_SEP2 = xm_SEP_store{k+1};
    xm_SEP_new = vertcat(xm_SEP1(2:end),xm_SEP2(2:end));
    xm_SEP_f  =  vertcat(xm_SEP_f,xm_SEP_new);
    
    xm_RV1 = xm_RV_store{k};
    xm_RV2 = xm_RV_store{k+1};
    xm_RV_new = vertcat(xm_RV1(2:end),xm_RV2(2:end));
    xm_RV_f  =  vertcat(xm_RV_f,xm_RV_new);
    
    ym1 = ym_store{k};
    ym2 = ym_store{k+1};
    ym_new = vertcat(ym1(2:end),ym2(2:end));
    ym_f  =  vertcat(ym_f,ym_new);
    
    SL_LV1 = SL_LV_store{k};
    SL_LV2 = SL_LV_store{k+1};
    SL_LV_new = vertcat(SL_LV1(2:end),SL_LV2(2:end));
    SL_LV_f  =  vertcat(SL_LV_f,SL_LV_new);
    
    SL_SEP1 = SL_SEP_store{k};
    SL_SEP2 = SL_SEP_store{k+1};
    SL_SEP_new = vertcat(SL_SEP1(2:end),SL_SEP2(2:end));
    SL_SEP_f  =  vertcat(SL_SEP_f,SL_SEP_new);
    
    SL_RV1 = SL_RV_store{k};
    SL_RV2 = SL_RV_store{k+1};
    SL_RV_new = vertcat(SL_RV1(2:end),SL_RV2(2:end));
    SL_RV_f  =  vertcat(SL_RV_f,SL_RV_new);
    
    V_LV1 = V_LV_store{k};
    V_LV2 = V_LV_store{k+1};
    V_LV_new = vertcat(V_LV1(2:end),V_LV2(2:end));
    V_LV_f  =  vertcat(V_LV_f,V_LV_new);
    
    V_RV1 = V_RV_store{k};
    V_RV2 = V_RV_store{k+1};
    V_RV_new = vertcat(V_RV1(2:end),V_RV2(2:end));
    V_RV_f  =  vertcat(V_RV_f,V_RV_new);
    
    A1_LV1 = A1_LV_store{k};
    A1_LV2 = A1_LV_store{k+1};
    A1_LV_new = vertcat(A1_LV1(2:end),A1_LV2(2:end));
    A1_LV_f  =  vertcat(A1_LV_f,A1_LV_new);
    
    A2_LV1 = A2_LV_store{k};
    A2_LV2 = A2_LV_store{k+1};
    A2_LV_new = vertcat(A2_LV1(2:end),A2_LV2(2:end));
    A2_LV_f  =  vertcat(A2_LV_f,A2_LV_new);
    
    A3_LV1 = A3_LV_store{k};
    A3_LV2 = A3_LV_store{k+1};
    A3_LV_new = vertcat(A3_LV1(2:end),A3_LV2(2:end));
    A3_LV_f  =  vertcat(A3_LV_f,A3_LV_new);
    
    N_LV1 = N_LV_store{k};
    N_LV2 = N_LV_store{k+1};
    N_LV_new = vertcat(N_LV1(2:end),N_LV2(2:end));
    N_LV_f  =  vertcat(N_LV_f,N_LV_new);
    
    A1_SEP1 = A1_SEP_store{k};
    A1_SEP2 = A1_SEP_store{k+1};
    A1_SEP_new = vertcat(A1_SEP1(2:end),A1_SEP2(2:end));
    A1_SEP_f  =  vertcat(A1_SEP_f,A1_SEP_new);
    
    A2_SEP1 = A2_SEP_store{k};
    A2_SEP2 = A2_SEP_store{k+1};
    A2_SEP_new = vertcat(A2_SEP1(2:end),A2_SEP2(2:end));
    A2_SEP_f  =  vertcat(A2_SEP_f,A2_SEP_new);
    
    A3_SEP1 = A3_SEP_store{k};
    A3_SEP2 = A3_SEP_store{k+1};
    A3_SEP_new = vertcat(A3_SEP1(2:end),A3_SEP2(2:end));
    A3_SEP_f  =  vertcat(A3_SEP_f,A3_SEP_new);
    
    N_SEP1 = N_SEP_store{k};
    N_SEP2 = N_SEP_store{k+1};
    N_SEP_new = vertcat(N_SEP1(2:end),N_SEP2(2:end));
    N_SEP_f  =  vertcat(N_SEP_f,N_SEP_new);
    
    A1_RV1 = A1_RV_store{k};
    A1_RV2 = A1_RV_store{k+1};
    A1_RV_new = vertcat(A1_RV1(2:end),A1_RV2(2:end));
    A1_RV_f  =  vertcat(A1_RV_f,A1_RV_new);
    
    A2_RV1 = A2_RV_store{k};
    A2_RV2 = A2_RV_store{k+1};
    A2_RV_new = vertcat(A2_RV1(2:end),A2_RV2(2:end));
    A2_RV_f  =  vertcat(A2_RV_f,A2_RV_new);
    
    A3_RV1 = A3_RV_store{k};
    A3_RV2 = A3_RV_store{k+1};
    A3_RV_new = vertcat(A3_RV1(2:end),A3_RV2(2:end));
    A3_RV_f  =  vertcat(A3_RV_f,A3_RV_new);
    
    N_RV1 = N_RV_store{k};
    N_RV2 = N_RV_store{k+1};
    N_RV_new = vertcat(N_RV1(2:end),N_RV2(2:end));
    N_RV_f  =  vertcat(N_RV_f,N_RV_new);
    
    V_SV1 = V_SV_store{k};
    V_SV2 = V_SV_store{k+1};
    V_SV_new = vertcat(V_SV1(2:end),V_SV2(2:end));
    V_SV_f  =  vertcat(V_SV_f,V_SV_new);
    
    V_PV1 = V_PV_store{k};
    V_PV2 = V_PV_store{k+1};
    V_PV_new = vertcat(V_PV1(2:end),V_PV2(2:end));
    V_PV_f  =  vertcat(V_PV_f,V_PV_new);
    
    V_SA1 = V_SA_store{k};
    V_SA2 = V_SA_store{k+1};
    V_SA_new = vertcat(V_SA1(2:end),V_SA2(2:end));
    V_SA_f  =  vertcat(V_SA_f,V_LV_new);
    
    V_PA1 = V_PA_store{k};
    V_PA2 = V_PA_store{k+1};
    V_PA_new = vertcat(V_PA1(2:end),V_PA2(2:end));
    V_PA_f  =  vertcat(V_PA_f,V_PA_new);
    
    V_Ao1 = V_Ao_store{k};
    V_Ao2 = V_Ao_store{k+1};
    V_Ao_new = vertcat(V_Ao1(2:end),V_Ao2(2:end));
    V_Ao_f  =  vertcat(V_Ao_f,V_Ao_new);
    
    Vm_LV1 = Vm_LV_store{k};
    Vm_LV2 = Vm_LV_store{k+1};
    Vm_LV_new = vertcat(Vm_LV1(2:end),Vm_LV2(2:end));
    Vm_LV_f  =  vertcat(Vm_LV_f,Vm_LV_new);
    
    Vm_SEP1 = Vm_SEP_store{k};
    Vm_SEP2 = Vm_SEP_store{k+1};
    Vm_SEP_new = vertcat(Vm_SEP1(2:end),Vm_SEP2(2:end));
    Vm_SEP_f  =  vertcat(Vm_SEP_f,Vm_SEP_new);
    
    Vm_RV1 = Vm_RV_store{k};
    Vm_RV2 = Vm_RV_store{k+1};
    Vm_RV_new = vertcat(Vm_RV1(2:end),Vm_RV2(2:end));
    Vm_RV_f  =  vertcat(Vm_RV_f,Vm_RV_new);
    
    Am_LV1 = Am_LV_store{k};
    Am_LV2 = Am_LV_store{k+1};
    Am_LV_new = vertcat(Am_LV1(2:end),Am_LV2(2:end));
    Am_LV_f  =  vertcat(Am_LV_f,Am_LV_new);
    
    Am_SEP1 = Am_SEP_store{k};
    Am_SEP2 = Am_SEP_store{k+1};
    Am_SEP_new = vertcat(Am_SEP1(2:end),Am_SEP2(2:end));
    Am_SEP_f  =  vertcat(Am_SEP_f,Am_SEP_new);
    
    Am_RV1 = Am_RV_store{k};
    Am_RV2 = Am_RV_store{k+1};
    Am_RV_new = vertcat(Am_RV1(2:end),Am_RV2(2:end));
    Am_RV_f  =  vertcat(Am_RV_f,Am_RV_new);
    
    Cm_LV1 = Cm_LV_store{k};
    Cm_LV2 = Cm_LV_store{k+1};
    Cm_LV_new = vertcat(Cm_LV1(2:end),Cm_LV2(2:end));
    Cm_LV_f  =  vertcat(Cm_LV_f,Cm_LV_new);
    
    Cm_SEP1 = Cm_SEP_store{k};
    Cm_SEP2 = Cm_SEP_store{k+1};
    Cm_SEP_new = vertcat(Cm_SEP1(2:end),Cm_SEP2(2:end));
    Cm_SEP_f  =  vertcat(Cm_SEP_f,Cm_SEP_new);
    
    Cm_RV1 = Cm_RV_store{k};
    Cm_RV2 = Cm_RV_store{k+1};
    Cm_RV_new = vertcat(Cm_RV1(2:end),Cm_RV2(2:end));
    Cm_RV_f  =  vertcat(Cm_RV_f,Cm_RV_new);
    
    z_LV1 = z_LV_store{k};
    z_LV2 = z_LV_store{k+1};
    z_LV_new = vertcat(z_LV1(2:end),z_LV2(2:end));
    z_LV_f  =  vertcat(z_LV_f,z_LV_new);
    
    z_SEP1 = z_SEP_store{k};
    z_SEP2 = z_SEP_store{k+1};
    z_SEP_new = vertcat(z_SEP1(2:end),z_SEP2(2:end));
    z_SEP_f  =  vertcat(z_SEP_f,z_SEP_new);
    
    z_RV1 = z_RV_store{k};
    z_RV2 = z_RV_store{k+1};
    z_RV_new = vertcat(z_RV1(2:end),z_RV2(2:end));
    z_RV_f  =  vertcat(z_RV_f,z_RV_new);
    
    epsf_LV1 = epsf_LV_store{k};
    epsf_LV2 = epsf_LV_store{k+1};
    epsf_LV_new = vertcat(epsf_LV1(2:end),epsf_LV2(2:end));
    epsf_LV_f  =  vertcat(epsf_LV_f,epsf_LV_new);
    
    epsf_SEP1 = epsf_SEP_store{k};
    epsf_SEP2 = epsf_SEP_store{k+1};
    epsf_SEP_new = vertcat(epsf_SEP1(2:end),epsf_SEP2(2:end));
    epsf_SEP_f  =  vertcat(epsf_SEP_f,epsf_SEP_new);
    
    epsf_RV1 = epsf_RV_store{k};
    epsf_RV2 = epsf_RV_store{k+1};
    epsf_RV_new = vertcat(epsf_RV1(2:end),epsf_RV2(2:end));
    epsf_RV_f  =  vertcat(epsf_RV_f,epsf_RV_new);
    
    SLo_LV1 = SLo_LV_store{k};
    SLo_LV2 = SLo_LV_store{k+1};
    SLo_LV_new = vertcat(SLo_LV1(2:end),SLo_LV2(2:end));
    SLo_LV_f  =  vertcat(SLo_LV_f,SLo_LV_new);
    
    SLo_SEP1 = SLo_SEP_store{k};
    SLo_SEP2 = SLo_SEP_store{k+1};
    SLo_SEP_new = vertcat(SLo_SEP1(2:end),SLo_SEP2(2:end));
    SLo_SEP_f  =  vertcat(SLo_SEP_f,SLo_SEP_new);
    
    SLo_RV1 = SLo_RV_store{k};
    SLo_RV2 = SLo_RV_store{k+1};
    SLo_RV_new = vertcat(SLo_RV1(2:end),SLo_RV2(2:end));
    SLo_RV_f  =  vertcat(SLo_RV_f,SLo_RV_new);
    
    sigmapas_LV1 = sigmapas_LV_store{k};
    sigmapas_LV2 = sigmapas_LV_store{k+1};
    sigmapas_LV_new = vertcat(sigmapas_LV1(2:end),sigmapas_LV2(2:end));
    sigmapas_LV_f  =  vertcat(sigmapas_LV_f,sigmapas_LV_new);
    
    sigmapas_SEP1 = sigmapas_SEP_store{k};
    sigmapas_SEP2 = sigmapas_SEP_store{k+1};
    sigmapas_SEP_new = vertcat(sigmapas_SEP1(2:end),sigmapas_SEP2(2:end));
    sigmapas_SEP_f  =  vertcat(sigmapas_SEP_f,sigmapas_SEP_new);
    
    sigmapas_RV1 = sigmapas_RV_store{k};
    sigmapas_RV2 = sigmapas_RV_store{k+1};
    sigmapas_RV_new = vertcat(sigmapas_RV1(2:end),sigmapas_RV2(2:end));
    sigmapas_RV_f  =  vertcat(sigmapas_RV_f,sigmapas_RV_new);
    
    sigmaact_LV1 = sigmaact_LV_store{k};
    sigmaact_LV2 = sigmaact_LV_store{k+1};
    sigmaact_LV_new = vertcat(sigmaact_LV1(2:end),sigmaact_LV2(2:end));
    sigmaact_LV_f  =  vertcat(sigmaact_LV_f,sigmaact_LV_new);
    
    sigmaact_SEP1 = sigmaact_SEP_store{k};
    sigmaact_SEP2 = sigmaact_SEP_store{k+1};
    sigmaact_SEP_new = vertcat(sigmaact_SEP1(2:end),sigmaact_SEP2(2:end));
    sigmaact_SEP_f  =  vertcat(sigmaact_SEP_f,sigmaact_SEP_new);
    
    sigmaact_RV1 = sigmaact_RV_store{k};
    sigmaact_RV2 = sigmaact_RV_store{k+1};
    sigmaact_RV_new = vertcat(sigmaact_RV1(2:end),sigmaact_RV2(2:end));
    sigmaact_RV_f  =  vertcat(sigmaact_RV_f,sigmaact_RV_new);
    
    sigmaf_LV1 = sigmaf_LV_store{k};
    sigmaf_LV2 = sigmaf_LV_store{k+1};
    sigmaf_LV_new = vertcat(sigmaf_LV1(2:end),sigmaf_LV2(2:end));
    sigmaf_LV_f  =  vertcat(sigmaf_LV_f,sigmaf_LV_new);
    
    sigmaf_SEP1 = sigmaf_SEP_store{k};
    sigmaf_SEP2 = sigmaf_SEP_store{k+1};
    sigmaf_SEP_new = vertcat(sigmaf_SEP1(2:end),sigmaf_SEP2(2:end));
    sigmaf_SEP_f  =  vertcat(sigmaf_SEP_f,sigmaf_SEP_new);
    
    sigmaf_RV1 = sigmaf_RV_store{k};
    sigmaf_RV2 = sigmaf_RV_store{k+1};
    sigmaf_RV_new = vertcat(sigmaf_RV1(2:end),sigmaf_RV2(2:end));
    sigmaf_RV_f  =  vertcat(sigmaf_RV_f,sigmaf_RV_new);
    
    Tm_LV1 = Tm_LV_store{k};
    Tm_LV2 = Tm_LV_store{k+1};
    Tm_LV_new = vertcat(Tm_LV1(2:end),Tm_LV2(2:end));
    Tm_LV_f  =  vertcat(Tm_LV_f,Tm_LV_new);
    
    Tm_SEP1 = Tm_SEP_store{k};
    Tm_SEP2 = Tm_SEP_store{k+1};
    Tm_SEP_new = vertcat(Tm_SEP1(2:end),Tm_SEP2(2:end));
    Tm_SEP_f  =  vertcat(Tm_SEP_f,Tm_SEP_new);
    
    Tm_RV1 = Tm_RV_store{k};
    Tm_RV2 = Tm_RV_store{k+1};
    Tm_RV_new = vertcat(Tm_RV1(2:end),Tm_RV2(2:end));
    Tm_RV_f  =  vertcat(Tm_RV_f,Tm_RV_new);
    
    Tx_LV1 = Tx_LV_store{k};
    Tx_LV2 = Tx_LV_store{k+1};
    Tx_LV_new = vertcat(Tx_LV1(2:end),Tx_LV2(2:end));
    Tx_LV_f  =  vertcat(Tx_LV_f,Tx_LV_new);
    
    Tx_SEP1 = Tx_SEP_store{k};
    Tx_SEP2 = Tx_SEP_store{k+1};
    Tx_SEP_new = vertcat(Tx_SEP1(2:end),Tx_SEP2(2:end));
    Tx_SEP_f  =  vertcat(Tx_SEP_f,Tx_SEP_new);
    
    Tx_RV1 = Tx_RV_store{k};
    Tx_RV2 = Tx_RV_store{k+1};
    Tx_RV_new = vertcat(Tx_RV1(2:end),Tx_RV2(2:end));
    Tx_RV_f  =  vertcat(Tx_RV_f,Tx_RV_new);
    
    Ty_LV1 = Ty_LV_store{k};
    Ty_LV2 = Ty_LV_store{k+1};
    Ty_LV_new = vertcat(Ty_LV1(2:end),Ty_LV2(2:end));
    Ty_LV_f  =  vertcat(Ty_LV_f,Ty_LV_new);
    
    Ty_SEP1 = Ty_SEP_store{k};
    Ty_SEP2 = Ty_SEP_store{k+1};
    Ty_SEP_new = vertcat(Ty_SEP1(2:end),Ty_SEP2(2:end));
    Ty_SEP_f  =  vertcat(Ty_SEP_f,Ty_SEP_new);
    
    Ty_RV1 = Ty_RV_store{k};
    Ty_RV2 = Ty_RV_store{k+1};
    Ty_RV_new = vertcat(Ty_RV1(2:end),Ty_RV2(2:end));
    Ty_RV_f  =  vertcat(Ty_RV_f,Ty_RV_new);
    
    P_LV1 = P_LV_store{k};
    P_LV2 = P_LV_store{k+1};
    P_LV_new = vertcat(P_LV1(2:end),P_LV2(2:end));
    P_LV_f  =  vertcat(P_LV_f,P_LV_new);
    
    P_RV1 = P_RV_store{k};
    P_RV2 = P_RV_store{k+1};
    P_RV_new = vertcat(P_RV1(2:end),P_RV2(2:end));
    P_RV_f  =  vertcat(P_RV_f,P_RV_new);
    
    P_PV1 = P_PV_store{k};
    P_PV2 = P_PV_store{k+1};
    P_PV_new = vertcat(P_PV1(2:end),P_PV2(2:end));
    P_PV_f  =  vertcat(P_PV_f,P_PV_new);
    
    P_SV1 = P_SV_store{k};
    P_SV2 = P_SV_store{k+1};
    P_SV_new = vertcat(P_SV1(2:end),P_SV2(2:end));
    P_SV_f  =  vertcat(P_SV_f,P_SV_new);
    
    P_PA1 = P_PA_store{k};
    P_PA2 = P_PA_store{k+1};
    P_PA_new = vertcat(P_PA1(2:end),P_PA2(2:end));
    P_PA_f  =  vertcat(P_PA_f,P_PA_new);
    
    P_SA1 = P_SA_store{k};
    P_SA2 = P_SA_store{k+1};
    P_SA_new = vertcat(P_SA1(2:end),P_SA2(2:end));
    P_SA_f  =  vertcat(P_SA_f,P_SA_new);
    
    P_Ao1 = P_Ao_store{k};
    P_Ao2 = P_Ao_store{k+1};
    P_Ao_new = vertcat(P_Ao1(2:end),P_Ao2(2:end));
    P_Ao_f  =  vertcat(P_Ao_f,P_Ao_new);
    
    k = k+2;
end

beats  = 1:420;
ATPase_orig = zeros(1,120);
for i = 1:length(ATPase_orig)
    ATPase_orig(i) = ATPase_store(1);
end
ATPase_full = horzcat(ATPase_orig,ATPase_store);

MVO2_tissue_orig = zeros(1,120);
for i = 1:length(MVO2_tissue_orig)
    MVO2_tissue_orig(i) = MVO2_tissue_store(1);
end
MVO2_tissue_full = horzcat(MVO2_tissue_orig,MVO2_tissue_store);

PCrATP_orig = zeros(1,120);
for i = 1:length(PCrATP_orig)
    PCrATP_orig(i) = PCrATP_store(1);
end
PCrATP_full = horzcat(PCrATP_orig,PCrATP_store);

ATP_orig = zeros(1,120);
for i = 1:length(ATP_orig)
    ATP_orig(i) = ATP_store(1);
end
ATP_full = horzcat(ATP_orig,ATP_store);

ADP_orig = zeros(1,120);
for i = 1:length(ADP_orig)
    ADP_orig(i) = ADP_store(1);
end
ADP_full = horzcat(ADP_orig,ADP_store);

Pi_orig = zeros(1,120);
for i = 1:length(Pi_orig)
    Pi_orig(i) = Pi_store(1);
end
Pi_full = horzcat(Pi_orig,Pi_store);


% Get just last beat
last_beat = 419; 

last_beat_t = last_beat*stim_period;
for i = 1:length(t_f)
    if t_f(i) > last_beat_t
        idx_last = i;
        break
    end
end


%% Calculate
t_f_last = t_f(idx_last:end);
xm_LV_last = xm_LV_f(idx_last:end);
xm_SEP_last = xm_SEP_f(idx_last:end);
xm_RV_last = xm_RV_f(idx_last:end);
ym_last = ym_f(idx_last:end);
SL_LV_last = SL_LV_f(idx_last:end);
SL_SEP_last = SL_SEP_f(idx_last:end);
SL_RV_last = SL_RV_f(idx_last:end);
V_LV_last = V_LV_f(idx_last:end);
V_RV_last = V_RV_f(idx_last:end);
A1_LV_last = A1_LV_f(idx_last:end);
A2_LV_last = A2_LV_f(idx_last:end);
A3_LV_last = A3_LV_f(idx_last:end);
N_LV_last = N_LV_f(idx_last:end);
A1_SEP_last = A1_SEP_f(idx_last:end);
A2_SEP_last = A2_SEP_f(idx_last:end);
A3_SEP_last = A3_SEP_f(idx_last:end);
N_SEP_last = N_SEP_f(idx_last:end);
A1_RV_last = A1_RV_f(idx_last:end);
A2_RV_last = A2_RV_f(idx_last:end);
A3_RV_last = A3_RV_f(idx_last:end);
N_RV_last = N_RV_f(idx_last:end);
V_SV_last = V_SV_f(idx_last:end);
V_PV_last = V_PV_f(idx_last:end);
V_SA_last = V_SA_f(idx_last:end);
V_PA_last = V_PA_f(idx_last:end);
V_Ao_last = V_Ao_f(idx_last:end);
Vm_LV_last = Vm_LV_f(idx_last:end);
Vm_SEP_last = Vm_SEP_f(idx_last:end);
Vm_RV_last = Vm_RV_f(idx_last:end);
Am_LV_last = Am_LV_f(idx_last:end);
Am_SEP_last = Am_SEP_f(idx_last:end);
Am_RV_last = Am_RV_f(idx_last:end);
Cm_LV_last = Cm_LV_f(idx_last:end);
Cm_SEP_last = Cm_SEP_f(idx_last:end);
Cm_RV_last = Cm_RV_f(idx_last:end);
z_LV_last = z_LV_f(idx_last:end);
z_SEP_last = z_SEP_f(idx_last:end);
z_RV_last = z_RV_f(idx_last:end);
epsf_LV_last = epsf_LV_f(idx_last:end);
epsf_SEP_last = epsf_SEP_f(idx_last:end);
epsf_RV_last = epsf_RV_f(idx_last:end);
SLo_LV_last = SLo_LV_f(idx_last:end);
SLo_SEP_last = SLo_SEP_f(idx_last:end);
SLo_RV_last = SLo_RV_f(idx_last:end);
sigmapas_LV_last = sigmapas_LV_f(idx_last:end);
sigmapas_SEP_last = sigmapas_SEP_f(idx_last:end);
sigmapas_RV_last = sigmapas_RV_f(idx_last:end);
sigmaact_LV_last = sigmaact_LV_f(idx_last:end);
sigmaact_SEP_last = sigmaact_SEP_f(idx_last:end);
sigmaact_RV_last = sigmaact_RV_f(idx_last:end);
sigmaf_LV_last = sigmaf_LV_f(idx_last:end);
sigmaf_SEP_last = sigmaf_SEP_f(idx_last:end);
sigmaf_RV_last = sigmaf_RV_f(idx_last:end);
Tm_LV_last = Tm_LV_f(idx_last:end);
Tm_SEP_last = Tm_SEP_f(idx_last:end);
Tm_RV_last = Tm_RV_f(idx_last:end);
Tx_LV_last = Tx_LV_f(idx_last:end);
Tx_SEP_last = Tx_SEP_f(idx_last:end);
Tx_RV_last = Tx_RV_f(idx_last:end);
Ty_LV_last = Tx_LV_f(idx_last:end);
Ty_SEP_last = Ty_SEP_f(idx_last:end);
Ty_RV_last = Ty_RV_f(idx_last:end);
P_LV_last = P_LV_f(idx_last:end);
P_RV_last = P_RV_f(idx_last:end);
P_PV_last = P_PV_f(idx_last:end);
P_SV_last = P_SV_f(idx_last:end);
P_PA_last = P_PA_f(idx_last:end);
P_SA_last = P_SA_f(idx_last:end);
P_Ao_last = P_Ao_f(idx_last:end);

max_force = max(sigmaf_LV_last);
FS = ((min(SL_LV_last)-max(SL_LV_last))/max(SL_LV_last))*100;
CO = (max(V_LV_last)-min(V_LV_last))*HR;
EF = (max(V_LV_last)-min(V_LV_last))/max(V_LV_last);

max_P = max(P_LV_last);
min_P = min(P_LV_last);
LVDP = max_P-min_P;

for i = 1:length(P_LV_last)
    if P_LV_last(i) >= max_P
        idx_max = i;
        break
    end
end

work_rate_F = SV_LV_sim*MAP*HR/60;
ATP_F = ATP_store(300);
ADP_F = ADP_store(300);
Pi_F = Pi_store(300);
MVO2_F = MVO2_tissue_store(300);
PCrATP_F = PCrATP_store(300);
XB_turnover_F = rate_of_XB_turnover_ave;
ATPase_F =  ATPase_store(300);
work_per_beat = ((SV_LV_sim*MAP*HR/60)/HR)*60;
efficiency_F = (work_per_beat/ATPase_store(300))*1000;
V_LV_store_F = V_LV_last;
P_LV_store_F = P_LV_last;

%end

end