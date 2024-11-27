%%%% Main code to run myocyte model
% Adapted by Abby Teitgen based on code from: 

% Tewari, S. G., Bugenhagen, S. M., Palmer, B. M, Beard,
% D. A. (2016). Dynamics of corss-bridge cycling, ATP hydrolysis, force
% generation, and deformation in cardiac muscle. J. Mol. Cell Cardiol. 96:
% 11-25.

% Tewari,S. G., Bugenhagen, S. M., Vinnakota, K. C., Rice, J. J., Janssen,
% M. L., Beard, D. A. (2016). Influence of metabolic dysfuntion on cardiac
% mechanics in decompensated hypertrophy and heart failure. J. Mol.
% Cell Cardiol. 94: 162-175.

% Lopez, R. Marzban, B., Gao, X. Lauinger, E., Van den Bergh, F.,
% Whitesall, S. E., Converso-Baran, K., Burant, C. F., Michele, D. E.,
% Beard, D. A. (2020). Impaired myocardial energetics causes mechanical dysfunction
% in decompensated failing hearts. Function, 1(2): zqaa018

% Marzban, B., Lopez, R., Beard, D. A. (2020). Computational modeling of coupled
% energetics and mechanics in the rat ventricular myocardium. Physiome.

% Note: code requires parallel computing toolbox for use of parallel for
% loops (alternatively, this can be removed and normal for loops can be
% used)


function [T_final_XB, force_final, idx_XB, Shortening_final, SS_Ftotal_fpca] = myocyte_model(Ca_value, XB, dATP_percent, kf_i, k_f_i, kw_i, k_w_i, kp_i, k_p_i, kg_i, krecruit_i, k_on_i, k_off_i, k_coop_i, a_i, b_i, c_i)
%% Flags
% Calcium
Ca_flag = Ca_value; % 0 = ATP, 1 = dATP

% Protocol
XB_protocol = XB; % 0 = force pCa, 1 = twitch

% Percent dATP
dATP = dATP_percent; % dATP fraction
ATP = 100 - dATP; % ATP fraction

% Calcium transient parameters
a = a_i;
b = b_i;
c = c_i;

%% Parameters
% Crossbridge
kf = 250*(ATP/100) + kf_i*(dATP/100); % Myosin actin associaiton rate (P to A1) (s^-1)
k_f = 304.6708*(ATP/100) + k_f_i*(dATP/100); % Myosin actin dissociation rate (A1 to P) (s^-1) 
kw = 112.3727*(ATP/100) + kw_i*(dATP/100); % A1 to A2 transition rate (s^-1)
k_w = 21.296*(ATP/100) + k_w_i*(dATP/100); % A2 to A1 transition rate (s^-1) 
kp = 811.72*(ATP/100) + kp_i*(dATP/100); % A2 to A3 transition rate (s^-1)
k_p = 43.25*(ATP/100) + k_p_i*(dATP/100); % A3 to A2 transition rate (s^-1) 
kg = 144.5586*(ATP/100) + kg_i*(dATP/100); % A3 to P transition rate (s^-1)
krecruit = 0.2069*(ATP/100) + krecruit_i*(dATP/100); % Force dependence of transition to OFF state (N^-1 m^-1)
km = 15.4691375; % OFF to ON transion rate (N^-1 m^-1)
k_m = 50.032; % ON to OFF transition rate (s^-1)
k_on = k_on_i; % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = k_off_i; % Rate constant of Ca2+ unbinding from troponin C (s^-1)
k_coop = k_coop_i; % Strength of thin filament cooperativity

visc = 0.001; % Viscosity (mmHg*s/um)
kstiff1 = 1.8219e+04; % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = 4.7822e+05; % Stiffness constant due to working stroke of XBs (kPa/um)
k_passive = 0.1; % Passive stiffness constant (kPa/um)
Lthin = 1.200; % Length of thin filament (nm)
Lthick = 1.670; % Length of thick filament (nm)
Lbare = 0.100; % Bare length of thin filament (nm)
dr = 0.01; % Power-stroke size (um)
L_rest_pas = 1.51; % Length at which passive force = 0 (um)
stim_period = 1; % Hz

% Metabolite concentrations based on mean sham rat (Lopez et al. 2020)
MgATP = 7.8873; % mM
MgADP = 0.0501; % mM
Pi = 1.2308; % mM

%% Run model
%% Force pCa
if XB_protocol == 0 % Force pCa
flag = 0;
Kse = 50000; % Series element stiffness (mmHg/um) 
SL0_fpca = 2.25; % (um) based on experimental protocol (Regnier et al. 2004)

% For input into ODE solver
para = [MgATP, MgADP, Pi, kstiff1, kstiff2, k_passive, Kse, k_coop, k_on, k_off, km, krecruit, k_m, kf, k_f, kw, k_w, kp, k_p, kg, visc, stim_period, a, b, c];

% Defining time vector
tspan_fpca = 0:0.001:0.3;

% Defining Ca range
Ca_fraction_fpca = [0.1:0.1:100]; % pCa 7 to 4

init_fpca  = [zeros(1, 10), SL0_fpca,0]; % Initial conditions for the model
init_fpca(10) = 1; % Setting the initial value for nonpermissible state equal to 1
SS_Ftotal_fpca = zeros(1, length(Ca_fraction_fpca)); % Initialize storage

parfor k = 1:length(Ca_fraction_fpca) % Parallel for loop
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 5e-3);
    [t_fpca, Y_fpca] = ode15s(@dXdT_myocyte_mechanics, tspan_fpca, init_fpca, options, para, Ca_fraction_fpca(k), flag, Ca_flag);
    
    // file_name = append('pCa_states_marcus_',num2str(Ca_fraction_fpca(k)),'.xls');
    // % writematrix(Y_fpca, file_name);
    % Get states
    p1_0_fpca = Y_fpca(:,1);
    A1_fpca{k} = p1_0_fpca;
    p2_0_fpca = Y_fpca(:,4);
    A2_fpca{k} = p2_0_fpca;
    p2_1_fpca = Y_fpca(:,5);
    p3_0_fpca = Y_fpca(:,7);
    A3_fpca{k} = p3_0_fpca;
    p3_1_fpca = Y_fpca(:,8);
    SL_fpca = Y_fpca(:,11);
    N_fpca = Y_fpca(:,10);
    N_fpca_store{k} = N_fpca;
    P_fpca{k} = 1 - N_fpca - p1_0_fpca - p2_0_fpca - p3_0_fpca;
    
    % Overlap function
    OV_Zaxis_fpca = min(Lthick/2, SL_fpca/2); % Overlap region closest to Z-axis (nm)
    OV_Mline_fpca = max(SL_fpca/2 - (SL_fpca - Lthin), Lbare/2); % Overlap region closest to M-line (nm)
    LOV_fpca = OV_Zaxis_fpca - OV_Mline_fpca; % Length of overlap (nm)
    N_overlap_thick_fpca = LOV_fpca*2/(Lthick - Lbare); % Fraction of thick filament overlap

    % Active force
    B_process_fpca = kstiff2*dr*p3_0_fpca;   % Force due to XB ratcheting
    C_process_fpca = kstiff1*(p2_1_fpca + p3_1_fpca); % Force due to stretching of XBs
    F_XB_fpca = N_overlap_thick_fpca.*(B_process_fpca + C_process_fpca); % Active force

    % Passive force
    gamma = 8;
    F_passive_fpca = k_passive.*(SL_fpca - L_rest_pas).^gamma; 
    
    % Total force
    Ftotal_fpca = F_XB_fpca + F_passive_fpca; 
    SS_Ftotal_fpca(k) = Ftotal_fpca(end)/10; % Store force output
end
// writematrix(SS_Ftotal_fpca', 'Marcus_force_pCa.xls');
end        

if XB_protocol == 1 % Twitch
flag = 1;
Kse = 35; %35.5; %34; % Series element stiffness (mmHg/um) 
SL0_twitch = 1.84; % (um) from Korte et al. 2011 experimental protocol

% For input into ODE solver
para = [MgATP, MgADP, Pi, kstiff1, kstiff2, k_passive, Kse, k_coop, k_on, k_off, km, krecruit, k_m, kf, k_f, kw, k_w, kp, k_p, kg, visc, stim_period, a, b, c];

beats = 3; % Number of beats
tspan_twitch = 0:0.00001:beats; % Time vector (each beat is 1000 ms)
last_beat = beats*1000-1000; % Point at which last beat starts (for plotting)
y_XB = [zeros(1, 10), SL0_twitch, 0]; % Initialize
y_XB(10) = 1; % Setting the initial value for nonpermissible state equal to 1
init_twitch = y_XB; % Initial conditions
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
Ca_i = 0;
[t_twitch, Y_twitch] = ode15s(@dXdT_myocyte_mechanics, tspan_twitch, init_twitch, options, para, Ca_i, flag, Ca_flag);

% Get states
SL_twitch = Y_twitch(:,11);
p1_0_twitch = Y_twitch(:,1);
A1_twitch = p1_0_twitch;
p2_0_twitch = Y_twitch(:,4);
A2_twitch = p2_0_twitch;
p2_1_twitch = Y_twitch(:,5);
p3_0_twitch = Y_twitch(:,7);
A3_twitch = p3_0_twitch;
p3_1_twitch = Y_twitch(:,8);
N_twitch = Y_twitch(:,10);
P_twitch = 1-N_twitch-A1_twitch-A2_twitch-A3_twitch;

% Overlap function
OV_Zaxis_twitch = min(Lthick/2, SL_twitch/2); % Overlap region closest to Z-axis (nm)
OV_Mline_twitch = max(SL_twitch/2 - (SL_twitch - Lthin), Lbare/2); % Overal region closest to M-line (nm)
LOV_twitch = OV_Zaxis_twitch - OV_Mline_twitch; % Length of overlap (nm)
N_overlap_thick_twitch = LOV_twitch*2/(Lthick - Lbare); % Fraction of thick filament overlap

% Active force
B_process_twitch = kstiff2 * dr * p3_0_twitch;   % Force due to XB ratcheting
C_process_twitch = kstiff1 * (p2_1_twitch + p3_1_twitch ); % Force due to stretching of XBs
F_XB_twitch = N_overlap_thick_twitch.*(B_process_twitch + C_process_twitch ); % Active force

% Passive force
gamma = 8;
F_passive_twitch = k_passive.*(SL_twitch - L_rest_pas).^gamma; 

% Total force
Ftotal_twitch = F_XB_twitch + F_passive_twitch; 

T_final_XB = t_twitch*1000; % Convert to ms
force_final = Ftotal_twitch/10; % Store force output
Shortening_final = SL_twitch; % Store shortening output

%% Get just last beat
for i = 1:length(T_final_XB)
    if T_final_XB(i) > last_beat
        idx_XB = i;
        break
    end
end
end


%% Calculate
if XB_protocol == 0 % Force pCa
[hill, ec50_n] = pCa_calculate((Ca_fraction_fpca')*10^(-6),(SS_Ftotal_fpca./max(SS_Ftotal_fpca))')
Ca50 = -log10(ec50_n)
peak_SS_force = max(SS_Ftotal_fpca)
        
else if XB_protocol == 1 % Twitch
force_final_last = force_final(idx_XB:end);
T_XB_final_last = T_final_XB(idx_XB:end);
Shortening_final_last = Shortening_final(idx_XB:end);

% Twitch
[min_force,idx_min] = min(force_final_last)
[max_force,idx_max] = max(force_final_last)
TTP_twitch = T_XB_final_last(idx_max) - last_beat
RT50 = (max(force_final_last)-min(force_final_last))*0.5 + min(force_final_last);
time_RT50_twitch = 0;
for i = idx_max:length(force_final_last)
    if force_final_last(i) < RT50
        time_RT50_twitch = ((T_XB_final_last(i) + T_XB_final_last(i-1))/2);
        time_RT50_twitch = time_RT50_twitch - last_beat - TTP_twitch
        break
    end
end

RT90 = (max(force_final_last)-min(force_final_last))*0.1 + min(force_final_last);
time_RT90_twitch = 0;
for i = idx_max:length(force_final_last)
    if force_final_last(i) < RT90
        time_RT90_twitch = ((T_XB_final_last(i) + T_XB_final_last(i-1))/2);
        time_RT90_twitch = time_RT90_twitch - last_beat - TTP_twitch
        break
    end
end

% Shortening
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i) + T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last) - min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i) + T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening
        break
    end
end

end
end


% For output
if XB_protocol == 0 % force pCa
Ftotal_ktr = 0;
t_ktr = 0;
T_final_XB = 0;
idx_XB = 0;
Shortening_final = 0;
force_final = 0;  

else if XB_protocol == 1 % Twitch    
SS_Ftotal_fpca = 0;
Ftotal_ktr = 0;
t_ktr = 0;
end
end


end