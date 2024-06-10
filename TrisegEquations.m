%%%% Calculates equations for triseg model

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

function f = TrisegEquations(x,Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV, V_LV, V_RV,Amref_LV,Amref_SEP,Amref_RV)

xm_LV  = x(1); % LV heart geometry variable, cm
xm_SEP = x(2); % Septum heart geometry variable, cm
xm_RV  = x(3); % RV heart geometry variable, cm
ym     = x(4); % Heart geometry variable, cm

Lsref = 1.9; % Resting SL, micron
Kse = 10000; 

%% Ventricular mechanics
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
SLo_LV = Lsref*exp(epsf_LV); SLo_SEP = Lsref*exp(epsf_SEP); SLo_RV = Lsref*exp(epsf_RV);

% Total forces
sigmaf_LV = Kse*(SLo_LV - SL_LV);
sigmaf_SEP = Kse*(SLo_SEP - SL_SEP);
sigmaf_RV = Kse*(SLo_RV - SL_RV);

% Tension
Tm_LV = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5); % LV midwall tension
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5); % SEP midwall tension
Tm_RV = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5); % RV midwall tension

Tx_LV = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2); % LV axial tension
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2); % SEP axial tension
Tx_RV = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2); % RV axial tension

Ty_LV = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2); % LV radial tension
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2); % SEP radial tension
Ty_RV = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2); % RV radial tension

f(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; % xm_LV
f(2) = (Tx_LV + Tx_SEP + Tx_RV); % xm_SEP
f(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; % xm_RV
f(4) = (Ty_LV + Ty_SEP + Ty_RV);  % ym



