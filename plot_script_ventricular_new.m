%%%% Plots ventricular level figures from:
% Teitgen, A. E., Hock, M. T., McCabe, K. J., Childers, M. C., Huber, G. A.,
% Marzban, B., Beard, D. A., McCammon, J. A., Regnier, M., McCulloch, A. D. (2024).
% Multiscale modeling shows how 2'-deoxy-ATP rescues ventricular function
% in heart failure. PNAS.


% Model inputs:
%Ca, HF, dATP_percent, kf_i, k_f_i, kw_i, k_w_i, kp_i, k_p_i, kg_i, k_recruit_i, k_on_i, k_off_i, k_coop_i, a, b, c

% Calcium transient parameters
a_ATP = 0.2968;
b_ATP = 1.5900;
c_ATP = 1.7171;

a_dATP = 0.2218;
b_dATP = 1.4980;
c_dATP = 1.9989;

%% Figure 5
% ATP
[V_LV_store_F_ATP, P_LV_store_F_ATP, max_force_ATP, FS_ATP, CO_ATP, EF_ATP, LVDP_ATP, work_rate_F_ATP, ATP_F_ATP, ADP_F_ATP, Pi_F_ATP, MVO2_F_ATP, PCrATP_F_ATP, XB_turnover_F_ATP, ATPase_F_ATP, efficiency_F_ATP] = CardiovascularMechanics(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+
[V_LV_store_F_dATP_kf, P_LV_store_F_dATP_kf, max_force_dATP_kf, FS_dATP_kf, CO_dATP_kf, EF_dATP_kf, LVDP_dATP_kf, work_rate_F_dATP_kf, ATP_F_dATP_kf, ADP_F_dATP_kf, Pi_F_dATP_kf, MVO2_F_dATP_kf, PCrATP_F_dATP_kf, XB_turnover_F_dATP_kf, ATPase_F_dATP_kf, efficiency_F_dATP_kf] = CardiovascularMechanics(0, 0, 1, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+
[V_LV_store_F_dATP_kf_k_f_kw, P_LV_store_F_dATP_kf_k_f_kw, max_force_dATP_kf_k_f_kw, FS_dATP_kf_k_f_kw, CO_dATP_kf_k_f_kw, EF_dATP_kf_k_f_kw, LVDP_dATP_kf_k_f_kw, work_rate_F_dATP_kf_k_f_kw, ATP_F_dATP_kf_k_f_kw, ADP_F_dATP_kf_k_f_kw, Pi_F_dATP_kf_k_f_kw, MVO2_F_dATP_kf_k_f_kw, PCrATP_F_dATP_kf_k_f_kw, XB_turnover_F_dATP_kf_k_f_kw, ATPase_F_dATP_kf_k_f_kw, efficiency_F_dATP_kf_k_f_kw] = CardiovascularMechanics(0, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP krecruit
[V_LV_store_F_dATP_krecruit, P_LV_store_F_dATP_krecruit, max_force_dATP_krecruit, FS_dATP_krecruit, CO_dATP_krecruit, EF_dATP_krecruit, LVDP_dATP_krecruit, work_rate_F_dATP_krecruit, ATP_F_dATP_krecruit, ADP_F_dATP_krecruit, Pi_F_dATP_krecruit, MVO2_F_dATP_krecruit, PCrATP_F_dATP_krecruit, XB_turnover_F_dATP_krecruit, ATPase_F_dATP_krecruit, efficiency_F_dATP_krecruit] = CardiovascularMechanics(0, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+ + krecruit
[V_LV_store_F_dATP_kf_k_f_kw_krecruit, P_LV_store_F_dATP_kf_k_f_kw_krecruit, max_force_dATP_kf_k_f_kw_krecruit, FS_dATP_kf_k_f_kw_krecruit, CO_dATP_kf_k_f_kw_krecruit, EF_dATP_kf_k_f_kw_krecruit, LVDP_dATP_kf_k_f_kw_krecruit, work_rate_F_dATP_kf_k_f_kw_krecruit, ATP_F_dATP_kf_k_f_kw_krecruit, ADP_F_dATP_kf_k_f_kw_krecruit, Pi_F_dATP_kf_k_f_kw_krecruit, MVO2_F_dATP_kf_k_f_kw_krecruit, PCrATP_F_dATP_kf_k_f_kw_krecruit, XB_turnover_F_dATP_kf_k_f_kw_krecruit, ATPase_F_dATP_kf_k_f_kw_krecruit, efficiency_F_dATP_kf_k_f_kw_krecruit] = CardiovascularMechanics(0, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP Ca
[V_LV_store_F_dATP_Ca, P_LV_store_F_dATP_Ca, max_force_dATP_Ca, FS_dATP_Ca, CO_dATP_Ca, EF_dATP_Ca, LVDP_dATP_Ca, work_rate_F_dATP_Ca, ATP_F_dATP_Ca, ADP_F_dATP_Ca, Pi_F_dATP_Ca, MVO2_F_dATP_Ca, PCrATP_F_dATP_Ca, XB_turnover_F_dATP_Ca, ATPase_F_dATP_Ca, efficiency_F_dATP_Ca] = CardiovascularMechanics(1, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP krecruit + Ca
[V_LV_store_F_dATP_krecruit_Ca, P_LV_store_F_dATP_krecruit_Ca, max_force_dATP_krecruit_Ca, FS_dATP_krecruit_Ca, CO_dATP_krecruit_Ca, EF_dATP_krecruit_Ca, LVDP_dATP_krecruit_Ca, work_rate_F_dATP_krecruit_Ca, ATP_F_dATP_krecruit_Ca, ADP_F_dATP_krecruit_Ca, Pi_F_dATP_krecruit_Ca, MVO2_F_dATP_krecruit_Ca, PCrATP_F_dATP_krecruit_Ca, XB_turnover_F_dATP_krecruit_Ca, ATPase_F_dATP_krecruit_Ca, efficiency_F_dATP_krecruit_Ca] = CardiovascularMechanics(1, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP kf+ + kf- + kw+ + krecruit + Ca
[V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca, P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca, max_force_dATP_kf_k_f_kw_krecruit_Ca, FS_dATP_kf_k_f_kw_krecruit_Ca, CO_dATP_kf_k_f_kw_krecruit_Ca, EF_dATP_kf_k_f_kw_krecruit_Ca, LVDP_dATP_kf_k_f_kw_krecruit_Ca, work_rate_F_dATP_kf_k_f_kw_krecruit_Ca, ATP_F_dATP_kf_k_f_kw_krecruit_Ca, ADP_F_dATP_kf_k_f_kw_krecruit_Ca, Pi_F_dATP_kf_k_f_kw_krecruit_Ca, MVO2_F_dATP_kf_k_f_kw_krecruit_Ca, PCrATP_F_dATP_kf_k_f_kw_krecruit_Ca, XB_turnover_F_dATP_kf_k_f_kw_krecruit_Ca, ATPase_F_dATP_kf_k_f_kw_krecruit_Ca, efficiency_F_dATP_kf_k_f_kw_krecruit_Ca] = CardiovascularMechanics(1, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);


figure
hold on
plot(V_LV_store_F_ATP,P_LV_store_F_ATP,'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca,P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca,'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','dATP')
set(gca, 'fontsize', 14)

figure
hold on
bar(0, ((EF_dATP_kf-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(1, ((EF_dATP_kf_k_f_kw-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(2, ((EF_dATP_krecruit-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(3, ((EF_dATP_kf_k_f_kw_krecruit-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(4, ((EF_dATP_Ca-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(5, ((EF_dATP_krecruit_Ca-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
bar(6, ((EF_dATP_kf_k_f_kw_krecruit_Ca-EF_ATP)./EF_ATP)*100,'facecolor',[0.5 0.5 0.5])
yline(((EF_ATP-EF_ATP)./EF_ATP)*100,':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(14,':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('EF')
set(gca,'FontSize',18,'XTick',[],'XTickLabel',{'','','','','','','',''})


%% Figure 6
[V_LV_store_F_ATP, P_LV_store_F_ATP, max_force_ATP, FS_ATP, CO_ATP, EF_ATP, LVDP_ATP, work_rate_F_ATP, ATP_F_ATP, ADP_F_ATP, Pi_F_ATP, MVO2_F_ATP, PCrATP_F_ATP, XB_turnover_F_ATP, ATPase_F_ATP, efficiency_F_ATP] = CardiovascularMechanics(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

dATP_sim = 1; %[0 1 2 20];

for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F] = CardiovascularMechanics(1, 1, dATP, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);
    V_LV_store_F_failure{i} = V_LV_store_F;
    P_LV_store_F_failure{i} = P_LV_store_F;
    max_force_failure{i} = max_force;
    FS_failure{i} = FS;
    CO_failure{i} = CO;
    EF_failure{i} = EF;
    LVDP_failure{i} = LVDP;
    work_rate_F_failure{i} = work_rate_F;
    ATP_F_failure{i} = ATP_F; 
    ADP_F_failure{i} = ADP_F;
    Pi_F_failure{i} = Pi_F;
    MVO2_F_failure{i} = MVO2_F;
    PCrATP_F_failure{i} = PCrATP_F;
    XB_turnover_F_failure{i} = XB_turnover_F;
    ATPase_F_failure{i} = ATPase_F;
    efficiency_F_failure{i} = efficiency_F;
end

figure
hold on
plot(V_LV_store_F_ATP,P_LV_store_F_ATP,':','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(V_LV_store_F_failure{1},P_LV_store_F_failure{1},':','linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(V_LV_store_F_failure{2},P_LV_store_F_failure{2},'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(V_LV_store_F_failure{3},P_LV_store_F_failure{3},'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(V_LV_store_F_failure{4},P_LV_store_F_failure{4},'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','HF','HF + 1% dATP','HF + 2% dATP','HF + 20% dATP')
set(gca, 'fontsize', 14)


figure
subplot(4,3,1)
hold on
bar(0,CO_ATP/344.75,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,CO_failure{1}/344.75,'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,CO_failure{2}/344.75,'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,CO_failure{3}/344.75,'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,CO_failure{4}/344.75,'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('SV (mL/beat)')
set(gca, 'fontsize', 16)

subplot(4,3,2)
hold on
bar(0,EF_ATP*100,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,EF_failure{1}*100,'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,EF_failure{2}*100,'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,EF_failure{3}*100,'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,EF_failure{4}*100,'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('EF (%)')
set(gca, 'fontsize', 16)

subplot(4,3,3)
hold on
bar(0,CO_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,CO_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,CO_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,CO_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,CO_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('CO (mL/min)')
set(gca, 'fontsize', 16)

subplot(4,3,4)
hold on
bar(0,LVDP_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,LVDP_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,LVDP_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,LVDP_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,LVDP_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('LVDevP (mmHg)')
set(gca, 'fontsize', 16)

subplot(4,3,5)
hold on
bar(0,work_rate_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,work_rate_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,work_rate_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,work_rate_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,work_rate_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('Work Rate (mmHg*mL/s)')
set(gca, 'fontsize', 16)

subplot(4,3,6)
hold on
bar(0,ATP_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,ATP_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,ATP_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,ATP_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,ATP_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ATP (mM)')
set(gca, 'fontsize', 16)

subplot(4,3,7)
hold on
bar(0,ADP_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,ADP_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,ADP_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,ADP_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,ADP_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ADP (\muM)')
set(gca, 'fontsize', 16)

subplot(4,3,8)
hold on
bar(0,Pi_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,Pi_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,Pi_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,Pi_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,Pi_F_failure{4},'facecolor',[0.6667+0.0299*7 0.2667+0.0888*6 0.6000-0.0222*7])
xlabel('Percent dATP')
ylabel('Pi (mM)')
set(gca, 'fontsize', 16)

subplot(4,3,9)
hold on
bar(0,MVO2_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,MVO2_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,MVO2_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,MVO2_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,MVO2_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('MVO_{2} (uM O_{2}/min/g tissue)')
set(gca, 'fontsize', 16)

subplot(4,3,10)
hold on
bar(0,PCrATP_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,PCrATP_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,PCrATP_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,PCrATP_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,PCrATP_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('CrP/ATP ratio')
set(gca, 'fontsize', 16)

subplot(4,3,11)
hold on
bar(0,ATPase_F_ATP,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,ATPase_F_failure{1},'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,ATPase_F_failure{2},'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,ATPase_F_failure{3},'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,ATPase_F_failure{4},'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ATPase (M/s/L)')
set(gca, 'fontsize', 16)

subplot(4,3,12)
hold on
bar(0,efficiency_F_ATP./1000,'facecolor',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
bar(1,efficiency_F_failure{1}./1000,'facecolor',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
bar(2,efficiency_F_failure{2}./1000,'facecolor',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
bar(3,efficiency_F_failure{3}./1000,'facecolor',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
bar(4,efficiency_F_failure{4}./1000,'facecolor',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('Efficiency (mL^{2}*mmHg*ms/M)')
set(gca, 'fontsize', 16)


%% Figure 7
% Healthy
% ATP
[V_LV_store_F_ATP, P_LV_store_F_ATP, max_force_ATP, FS_ATP, CO_ATP, EF_ATP, LVDP_ATP, work_rate_F_ATP, ATP_F_ATP, ADP_F_ATP, Pi_F_ATP, MVO2_F_ATP, PCrATP_F_ATP, XB_turnover_F_ATP, ATPase_F_ATP, efficiency_F_ATP] = CardiovascularMechanics(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+
[V_LV_store_F_dATP_kf, P_LV_store_F_dATP_kf, max_force_dATP_kf, FS_dATP_kf, CO_dATP_kf, EF_dATP_kf, LVDP_dATP_kf, work_rate_F_dATP_kf, ATP_F_dATP_kf, ADP_F_dATP_kf, Pi_F_dATP_kf, MVO2_F_dATP_kf, PCrATP_F_dATP_kf, XB_turnover_F_dATP_kf, ATPase_F_dATP_kf, efficiency_F_dATP_kf] = CardiovascularMechanics(0, 0, 1, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+
[V_LV_store_F_dATP_kf_k_f_kw, P_LV_store_F_dATP_kf_k_f_kw, max_force_dATP_kf_k_f_kw, FS_dATP_kf_k_f_kw, CO_dATP_kf_k_f_kw, EF_dATP_kf_k_f_kw, LVDP_dATP_kf_k_f_kw, work_rate_F_dATP_kf_k_f_kw, ATP_F_dATP_kf_k_f_kw, ADP_F_dATP_kf_k_f_kw, Pi_F_dATP_kf_k_f_kw, MVO2_F_dATP_kf_k_f_kw, PCrATP_F_dATP_kf_k_f_kw, XB_turnover_F_dATP_kf_k_f_kw, ATPase_F_dATP_kf_k_f_kw, efficiency_F_dATP_kf_k_f_kw] = CardiovascularMechanics(0, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP krecruit
[V_LV_store_F_dATP_krecruit, P_LV_store_F_dATP_krecruit, max_force_dATP_krecruit, FS_dATP_krecruit, CO_dATP_krecruit, EF_dATP_krecruit, LVDP_dATP_krecruit, work_rate_F_dATP_krecruit, ATP_F_dATP_krecruit, ADP_F_dATP_krecruit, Pi_F_dATP_krecruit, MVO2_F_dATP_krecruit, PCrATP_F_dATP_krecruit, XB_turnover_F_dATP_krecruit, ATPase_F_dATP_krecruit, efficiency_F_dATP_krecruit] = CardiovascularMechanics(0, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+ + krecruit
[V_LV_store_F_dATP_kf_k_f_kw_krecruit, P_LV_store_F_dATP_kf_k_f_kw_krecruit, max_force_dATP_kf_k_f_kw_krecruit, FS_dATP_kf_k_f_kw_krecruit, CO_dATP_kf_k_f_kw_krecruit, EF_dATP_kf_k_f_kw_krecruit, LVDP_dATP_kf_k_f_kw_krecruit, work_rate_F_dATP_kf_k_f_kw_krecruit, ATP_F_dATP_kf_k_f_kw_krecruit, ADP_F_dATP_kf_k_f_kw_krecruit, Pi_F_dATP_kf_k_f_kw_krecruit, MVO2_F_dATP_kf_k_f_kw_krecruit, PCrATP_F_dATP_kf_k_f_kw_krecruit, XB_turnover_F_dATP_kf_k_f_kw_krecruit, ATPase_F_dATP_kf_k_f_kw_krecruit, efficiency_F_dATP_kf_k_f_kw_krecruit] = CardiovascularMechanics(0, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP Ca
[V_LV_store_F_dATP_Ca, P_LV_store_F_dATP_Ca, max_force_dATP_Ca, FS_dATP_Ca, CO_dATP_Ca, EF_dATP_Ca, LVDP_dATP_Ca, work_rate_F_dATP_Ca, ATP_F_dATP_Ca, ADP_F_dATP_Ca, Pi_F_dATP_Ca, MVO2_F_dATP_Ca, PCrATP_F_dATP_Ca, XB_turnover_F_dATP_Ca, ATPase_F_dATP_Ca, efficiency_F_dATP_Ca] = CardiovascularMechanics(1, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP krecruit + Ca
[V_LV_store_F_dATP_krecruit_Ca, P_LV_store_F_dATP_krecruit_Ca, max_force_dATP_krecruit_Ca, FS_dATP_krecruit_Ca, CO_dATP_krecruit_Ca, EF_dATP_krecruit_Ca, LVDP_dATP_krecruit_Ca, work_rate_F_dATP_krecruit_Ca, ATP_F_dATP_krecruit_Ca, ADP_F_dATP_krecruit_Ca, Pi_F_dATP_krecruit_Ca, MVO2_F_dATP_krecruit_Ca, PCrATP_F_dATP_krecruit_Ca, XB_turnover_F_dATP_krecruit_Ca, ATPase_F_dATP_krecruit_Ca, efficiency_F_dATP_krecruit_Ca] = CardiovascularMechanics(1, 0, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP kf+ + kf- + kw+ + krecruit + Ca
[V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca, P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca, max_force_dATP_kf_k_f_kw_krecruit_Ca, FS_dATP_kf_k_f_kw_krecruit_Ca, CO_dATP_kf_k_f_kw_krecruit_Ca, EF_dATP_kf_k_f_kw_krecruit_Ca, LVDP_dATP_kf_k_f_kw_krecruit_Ca, work_rate_F_dATP_kf_k_f_kw_krecruit_Ca, ATP_F_dATP_kf_k_f_kw_krecruit_Ca, ADP_F_dATP_kf_k_f_kw_krecruit_Ca, Pi_F_dATP_kf_k_f_kw_krecruit_Ca, MVO2_F_dATP_kf_k_f_kw_krecruit_Ca, PCrATP_F_dATP_kf_k_f_kw_krecruit_Ca, XB_turnover_F_dATP_kf_k_f_kw_krecruit_Ca, ATPase_F_dATP_kf_k_f_kw_krecruit_Ca, efficiency_F_dATP_kf_k_f_kw_krecruit_Ca] = CardiovascularMechanics(1, 0, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);


% Failure
% ATP
[V_LV_store_F_ATP_HF, P_LV_store_F_ATP_HF, max_force_ATP_HF, FS_ATP_HF, CO_ATP_HF, EF_ATP_HF, LVDP_ATP_HF, work_rate_F_ATP_HF, ATP_F_ATP_HF, ADP_F_ATP_HF, Pi_F_ATP_HF, MVO2_F_ATP_HF, PCrATP_F_ATP_HF, XB_turnover_F_ATP_HF, ATPase_F_ATP_HF, efficiency_F_ATP_HF] = CardiovascularMechanics(0, 1, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+
[V_LV_store_F_dATP_kf_HF, P_LV_store_F_dATP_kf_HF, max_force_dATP_kf_HF, FS_dATP_kf_HF, CO_dATP_kf_HF, EF_dATP_kf_HF, LVDP_dATP_kf_HF, work_rate_F_dATP_kf_HF, ATP_F_dATP_kf_HF, ADP_F_dATP_kf_HF, Pi_F_dATP_kf_HF, MVO2_F_dATP_kf_HF, PCrATP_F_dATP_kf_HF, XB_turnover_F_dATP_kf_HF, ATPase_F_dATP_kf_HF, efficiency_F_dATP_kf_HF] = CardiovascularMechanics(0, 1, 1, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+
[V_LV_store_F_dATP_kf_k_f_kw_HF, P_LV_store_F_dATP_kf_k_f_kw_HF, max_force_dATP_kf_k_f_kw_HF, FS_dATP_kf_k_f_kw_HF, CO_dATP_kf_k_f_kw_HF, EF_dATP_kf_k_f_kw_HF, LVDP_dATP_kf_k_f_kw_HF, work_rate_F_dATP_kf_k_f_kw_HF, ATP_F_dATP_kf_k_f_kw_HF, ADP_F_dATP_kf_k_f_kw_HF, Pi_F_dATP_kf_k_f_kw_HF, MVO2_F_dATP_kf_k_f_kw_HF, PCrATP_F_dATP_kf_k_f_kw_HF, XB_turnover_F_dATP_kf_k_f_kw_HF, ATPase_F_dATP_kf_k_f_kw_HF, efficiency_F_dATP_kf_k_f_kw_HF] = CardiovascularMechanics(0, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP krecruit
[V_LV_store_F_dATP_krecruit_HF, P_LV_store_F_dATP_krecruit_HF, max_force_dATP_krecruit_HF, FS_dATP_krecruit_HF, CO_dATP_krecruit_HF, EF_dATP_krecruit_HF, LVDP_dATP_krecruit_HF, work_rate_F_dATP_krecruit_HF, ATP_F_dATP_krecruit_HF, ADP_F_dATP_krecruit_HF, Pi_F_dATP_krecruit_HF, MVO2_F_dATP_krecruit_HF, PCrATP_F_dATP_krecruit_HF, XB_turnover_F_dATP_krecruit_HF, ATPase_F_dATP_krecruit_HF, efficiency_F_dATP_krecruit_HF] = CardiovascularMechanics(0, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP kf+ + kf- + kw+ + krecruit
[V_LV_store_F_dATP_kf_k_f_kw_krecruit_HF, P_LV_store_F_dATP_kf_k_f_kw_krecruit_HF, max_force_dATP_kf_k_f_kw_krecruit_HF, FS_dATP_kf_k_f_kw_krecruit_HF, CO_dATP_kf_k_f_kw_krecruit_HF, EF_dATP_kf_k_f_kw_krecruit_HF, LVDP_dATP_kf_k_f_kw_krecruit_HF, work_rate_F_dATP_kf_k_f_kw_krecruit_HF, ATP_F_dATP_kf_k_f_kw_krecruit_HF, ADP_F_dATP_kf_k_f_kw_krecruit_HF, Pi_F_dATP_kf_k_f_kw_krecruit_HF, MVO2_F_dATP_kf_k_f_kw_krecruit_HF, PCrATP_F_dATP_kf_k_f_kw_krecruit_HF, XB_turnover_F_dATP_kf_k_f_kw_krecruit_HF, ATPase_F_dATP_kf_k_f_kw_krecruit_HF, efficiency_F_dATP_kf_k_f_kw_krecruit_HF] = CardiovascularMechanics(0, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% dATP Ca
[V_LV_store_F_dATP_Ca_HF, P_LV_store_F_Ca_HF, max_force_dATP_Ca_HF, FS_dATP_Ca_HF, CO_dATP_Ca_HF, EF_dATP_Ca_HF, LVDP_dATP_Ca_HF, work_rate_F_dATP_Ca_HF, ATP_F_dATP_Ca_HF, ADP_F_dATP_Ca_HF, Pi_F_dATP_Ca_HF, MVO2_F_dATP_Ca_HF, PCrATP_F_dATP_Ca_HF, XB_turnover_F_dATP_Ca_HF, ATPase_F_dATP_Ca_HF, efficiency_F_dATP_Ca_HF] = CardiovascularMechanics(1, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP krecruit + Ca
[V_LV_store_F_dATP_krecruit_Ca_HF, P_LV_store_F_dATP_krecruit_Ca_HF, max_force_dATP_krecruit_Ca_HF, FS_dATP_krecruit_Ca_HF, CO_dATP_krecruit_Ca_HF, EF_dATP_krecruit_Ca_HF, LVDP_dATP_krecruit_Ca_HF, work_rate_F_dATP_krecruit_Ca_HF, ATP_F_dATP_krecruit_Ca_HF, ADP_F_dATP_krecruit_Ca_HF, Pi_F_dATP_krecruit_Ca_HF, MVO2_F_dATP_krecruit_Ca_HF, PCrATP_F_dATP_krecruit_Ca_HF, XB_turnover_F_dATP_krecruit_Ca_HF, ATPase_F_dATP_krecruit_Ca_HF, efficiency_F_dATP_krecruit_Ca_HF] = CardiovascularMechanics(1, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);

%dATP kf+ + kf- + kw+ + krecruit + Ca
[V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca_HF, P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca_HF, max_force_dATP_kf_k_f_kw_krecruit_Ca_HF, FS_dATP_kf_k_f_kw_krecruit_Ca_HF, CO_dATP_kf_k_f_kw_krecruit_Ca_HF, EF_dATP_kf_k_f_kw_krecruit_Ca_HF, LVDP_dATP_kf_k_f_kw_krecruit_Ca_HF, work_rate_F_dATP_kf_k_f_kw_krecruit_Ca_HF, ATP_F_dATP_kf_k_f_kw_krecruit_Ca_HF, ADP_F_dATP_kf_k_f_kw_krecruit_Ca_HF, Pi_F_dATP_kf_k_f_kw_krecruit_Ca_HF, MVO2_F_dATP_kf_k_f_kw_krecruit_Ca_HF, PCrATP_F_dATP_kf_k_f_kw_krecruit_Ca_HF, XB_turnover_F_dATP_kf_k_f_kw_krecruit_Ca_HF, ATPase_F_dATP_kf_k_f_kw_krecruit_Ca_HF, efficiency_F_dATP_kf_k_f_kw_krecruit_Ca_HF] = CardiovascularMechanics(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);


Y1 = [((EF_dATP_kf-EF_ATP)./EF_ATP)*100, ((EF_dATP_kf_HF-EF_ATP_HF)./EF_ATP_HF)*100
    ((EF_dATP_kf_k_f_kw-EF_ATP)./EF_ATP)*100, ((EF_dATP_kf_k_f_kw_HF-EF_ATP_HF)./EF_ATP_HF)*100 
    ((EF_dATP_krecruit-EF_ATP)./EF_ATP)*100, ((EF_dATP_krecruit_HF-EF_ATP_HF)./EF_ATP_HF)*100
    ((EF_dATP_kf_k_f_kw_krecruit-EF_ATP)./EF_ATP)*100, ((EF_dATP_kf_k_f_kw_krecruit_HF-EF_ATP_HF)./EF_ATP_HF)*100
    ((EF_dATP_Ca-EF_ATP)./EF_ATP)*100, ((EF_dATP_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100
    ((EF_dATP_krecruit_Ca-EF_ATP)./EF_ATP)*100, ((EF_dATP_krecruit_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100
    ((EF_dATP_kf_k_f_kw_krecruit_Ca-EF_ATP)./EF_ATP)*100, ((EF_dATP_kf_k_f_kw_krecruit_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100];  
figure
subplot(4,2,1)
hold on
bar(Y1)
ylabel('EF')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y2 = [((CO_dATP_kf-CO_ATP)./CO_ATP)*100, ((CO_dATP_kf_HF-CO_ATP_HF)./CO_ATP_HF)*100 
    ((CO_dATP_kf_k_f_kw-CO_ATP)./CO_ATP)*100, ((CO_dATP_kf_k_f_kw_HF-CO_ATP_HF)./CO_ATP_HF)*100
    ((CO_dATP_kf_k_f_kw_krecruit-CO_ATP)./CO_ATP)*100, ((CO_dATP_kf_k_f_kw_krecruit_HF-CO_ATP_HF)./CO_ATP_HF)*100
    ((CO_dATP_krecruit-CO_ATP)./CO_ATP)*100, ((CO_dATP_krecruit_HF-CO_ATP_HF)./CO_ATP_HF)*100 
    ((CO_dATP_Ca-CO_ATP)./CO_ATP)*100, ((CO_dATP_Ca_HF-CO_ATP_HF)./CO_ATP_HF)*100
    ((CO_dATP_krecruit_Ca-CO_ATP)./CO_ATP)*100, ((CO_dATP_krecruit_Ca_HF-CO_ATP_HF)./CO_ATP_HF)*100
    ((CO_dATP_kf_k_f_kw_krecruit_Ca-CO_ATP)./CO_ATP)*100, ((CO_dATP_kf_k_f_kw_krecruit_Ca_HF-CO_ATP_HF)./CO_ATP_HF)*100];
subplot(4,2,2)
hold on
bar(Y2)
ylabel('CO')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y3 = [((LVDP_dATP_kf-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_kf_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*10 
    ((LVDP_dATP_kf_k_f_kw-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_kf_k_f_kw_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100
    ((LVDP_dATP_krecruit-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_krecruit_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100
    ((LVDP_dATP_kf_k_f_kw_krecruit-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_kf_k_f_kw_krecruit_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100
    ((LVDP_dATP_Ca-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_Ca_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100
    ((LVDP_dATP_krecruit_Ca-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_krecruit_Ca_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100
    ((LVDP_dATP_kf_k_f_kw_krecruit_Ca-LVDP_ATP)./LVDP_ATP)*100, ((LVDP_dATP_kf_k_f_kw_krecruit_Ca_HF-LVDP_ATP_HF)./LVDP_ATP_HF)*100];
subplot(4,2,3)
hold on
bar(Y3)
ylabel('LVDevP')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y4 = [((work_rate_F_dATP_kf-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_kf_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100
    ((work_rate_F_dATP_kf_k_f_kw-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_kf_k_f_kw_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100 
    ((work_rate_F_dATP_kf_k_f_kw_krecruit-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_kf_k_f_kw_krecruit_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100 
    ((work_rate_F_dATP_krecruit-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_krecruit_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100
    ((work_rate_F_dATP_Ca_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100, ((work_rate_F_dATP_Ca-work_rate_F_ATP)./work_rate_F_ATP)*100 
    ((work_rate_F_dATP_krecruit_Ca-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_krecruit_Ca_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100
    ((work_rate_F_dATP_kf_k_f_kw_krecruit_Ca-work_rate_F_ATP)./work_rate_F_ATP)*100, ((work_rate_F_dATP_kf_k_f_kw_krecruit_Ca_HF-work_rate_F_ATP_HF)./work_rate_F_ATP_HF)*100];
subplot(4,2,4)
hold on
bar(Y4)
ylabel('Work Rate')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y5 = [((MVO2_F_dATP_kf-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_kf_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100 
    ((MVO2_F_dATP_kf_k_f_kw-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_kf_k_f_kw_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100
    ((MVO2_F_dATP_kf_k_f_kw_krecruit-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_kf_k_f_kw_krecruit_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100
    ((MVO2_F_dATP_krecruit-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_krecruit_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100
    ((MVO2_F_dATP_Ca-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_Ca_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100
    ((MVO2_F_dATP_krecruit_Ca-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_krecruit_Ca_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100
    ((MVO2_F_dATP_kf_k_f_kw_krecruit_Ca-MVO2_F_ATP)./MVO2_F_ATP)*100, ((MVO2_F_dATP_kf_k_f_kw_krecruit_Ca_HF-MVO2_F_ATP_HF)./MVO2_F_ATP_HF)*100];
subplot(4,2,5)
hold on
bar(Y5)
ylabel('MVO_{2}')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y6 = [((PCrATP_F_dATP_kf-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_kf_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100 
    ((PCrATP_F_dATP_kf_k_f_kw-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_kf_k_f_kw_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100
    ((PCrATP_F_dATP_kf_k_f_kw_krecruit-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_kf_k_f_kw_krecruit_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100
    ((PCrATP_F_dATP_krecruit-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_krecruit_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100
    ((PCrATP_F_dATP_Ca-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_Ca_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100
    ((PCrATP_F_dATP_krecruit_Ca-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_krecruit_Ca_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100
    ((PCrATP_F_dATP_kf_k_f_kw_krecruit_Ca-PCrATP_F_ATP)./PCrATP_F_ATP)*100, ((PCrATP_F_dATP_kf_k_f_kw_krecruit_Ca_HF-PCrATP_F_ATP_HF)./PCrATP_F_ATP_HF)*100];
subplot(4,2,6)
hold on
bar(Y6)
ylabel('CrP/ATP ratio')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y7 = [((ATPase_F_dATP_kf-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_kf_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100 
    ((ATPase_F_dATP_kf_k_f_kw-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_kf_k_f_kw_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100
    ((ATPase_F_dATP_kf_k_f_kw_krecruit-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_kf_k_f_kw_krecruit_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100
    ((ATPase_F_dATP_krecruit-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_krecruit_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100
    ((ATPase_F_dATP_Ca-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_Ca_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100
    ((ATPase_F_dATP_krecruit_Ca-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_krecruit_Ca_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100
    ((ATPase_F_dATP_kf_k_f_kw_krecruit_Ca-ATPase_F_ATP)./ATPase_F_ATP)*100, ((ATPase_F_dATP_kf_k_f_kw_krecruit_Ca_HF-ATPase_F_ATP_HF)./ATPase_F_ATP_HF)*100];
subplot(4,2,7)
hold on
bar(Y7)
ylabel('ATPase')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

Y8 = [((efficiency_F_dATP_kf-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_kf_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_kf_k_f_kw-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_kf_k_f_kw_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_kf_k_f_kw_krecruit-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_kf_k_f_kw_krecruit_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_krecruit-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_krecruit_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_Ca-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_Ca_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_krecruit_Ca-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_krecruit_Ca_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100
    ((efficiency_F_dATP_kf_k_f_kw_krecruit_Ca-efficiency_F_ATP)./efficiency_F_ATP)*100, ((efficiency_F_dATP_kf_k_f_kw_krecruit_Ca_HF-efficiency_F_ATP_HF)./efficiency_F_ATP_HF)*100];
subplot(4,2,8)
hold on
bar(Y8)
ylabel('Efficiency')
set(gca,'FontSize',16,'XTick',[],'XTickLabel',{'','','','','','','',''})

%% Figure S19 (Run simulations from Figure 5 to generate)
figure
hold on
plot(V_LV_store_F_ATP,P_LV_store_F_ATP,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(V_LV_store_F_dATP_kf,P_LV_store_F_dATP_kf,'linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_dATP_kf_k_f_kw,P_LV_store_F_dATP_kf_k_f_kw,':','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_dATP_krecruit,P_LV_store_F_dATP_krecruit,'-.','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_dATP_Ca,P_LV_store_F_dATP_Ca,'linewidth',3,'color',[0.8 0.8 0.8])
plot(V_LV_store_F_dATP_krecruit_Ca,P_LV_store_F_dATP_krecruit_Ca,'--','linewidth',3,'color',[0.5 0.5 0.5])
plot(V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca,P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca, 'fontsize', 14)
legend('ATP','dATP k_f', 'dATP k_f^+ + k_f^- + k_w^+','dATP k_{recruit}','dATP Ca','dATP k_{recruit} + Ca','dATP k_f^+ + k_f^- + k_w^+ + k_{recruit} + Ca')

%% Figure S20 (run simulations from Figure 7 to generate)
figure
hold on
plot(V_LV_store_F_ATP_HF,P_LV_store_F_ATP_HF,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(V_LV_store_F_dATP_kf_HF,P_LV_store_F_dATP_kf_HF,'linewidth',3,'color',[0.5 0.5 0.5])
plot(V_LV_store_F_dATP_kf_k_f_kw_HF,P_LV_store_F_dATP_kf_k_f_kw_HF,':','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_dATP_krecruit_HF,P_LV_store_F_dATP_krecruit_HF,'-.','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_dATP_Ca_HF,P_LV_store_F_Ca_HF,'linewidth',3,'color',[0.8 0.8 0.8])
plot(V_LV_store_F_dATP_krecruit_Ca_HF,P_LV_store_F_dATP_krecruit_Ca_HF,'--','linewidth',3,'color',[0.5 0.5 0.5])
plot(V_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca_HF,P_LV_store_F_dATP_kf_k_f_kw_krecruit_Ca_HF,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca, 'fontsize', 14)
legend('HF','dATP k_f', 'dATP k_f^+ + k_f^- + k_w^+','dATP k_{recruit}','dATP Ca','dATP k_{recruit} + Ca','dATP k_f^+ + k_f^- + k_w^+ + k_{recruit} + Ca')

figure
hold on
bar(0, ((EF_dATP_kf_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(1, ((EF_dATP_kf_k_f_kw_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(2, ((EF_dATP_krecruit_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(3, ((EF_dATP_kf_k_f_kw_krecruit_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(4, ((EF_dATP_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(5, ((EF_dATP_krecruit_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
bar(6, ((EF_dATP_kf_k_f_kw_krecruit_Ca_HF-EF_ATP_HF)./EF_ATP_HF)*100,'facecolor',[0.5 0.5 0.5])
yline(((EF_ATP_HF-EF_ATP_HF)./EF_ATP_HF)*100,':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
ylabel('EF')
set(gca,'FontSize',18,'XTick',[],'XTickLabel',{'','','','','','','',''})

%% Figure S21
% Healthy
dATP_sim = [0 1 2 4 6 8 10 20];
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F] = CardiovascularMechanics(1, 0, dATP, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);
    V_LV_store_F_healthy{i} = V_LV_store_F;
    P_LV_store_F_healthy{i} = P_LV_store_F;
    max_force_healthy{i} = max_force;
    FS_healthy{i} = FS;
    CO_healthy{i} = CO;
    EF_healthy{i} = EF;
    LVDP_healthy{i} = LVDP;
    work_rate_F_healthy{i} = work_rate_F;
    ATP_F_healthy{i} = ATP_F; 
    ADP_F_healthy{i} = ADP_F;
    Pi_F_healthy{i} = Pi_F;
    MVO2_F_healthy{i} = MVO2_F;
    PCrATP_F_healthy{i} = PCrATP_F;
    XB_turnover_F_healthy{i} = XB_turnover_F;
    ATPase_F_healthy{i} = ATPase_F;
    efficiency_F_healthy{i} = efficiency_F;
end

% Failue
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F] = CardiovascularMechanics(1, 1, dATP, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 200, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP);
    V_LV_store_F_failure{i} = V_LV_store_F;
    P_LV_store_F_failure{i} = P_LV_store_F;
    max_force_failure{i} = max_force;
    FS_failure{i} = FS;
    CO_failure{i} = CO;
    EF_failure{i} = EF;
    LVDP_failure{i} = LVDP;
    work_rate_F_failure{i} = work_rate_F;
    ATP_F_failure{i} = ATP_F; 
    ADP_F_failure{i} = ADP_F;
    Pi_F_failure{i} = Pi_F;
    MVO2_F_failure{i} = MVO2_F;
    PCrATP_F_failure{i} = PCrATP_F;
    XB_turnover_F_failure{i} = XB_turnover_F;
    ATPase_F_failure{i} = ATPase_F;
    efficiency_F_failure{i} = efficiency_F;
end

figure
hold on
plot(V_LV_store_F_healthy{1},P_LV_store_F_healthy{1},':','linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(V_LV_store_F_healthy{2},P_LV_store_F_healthy{2},'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(V_LV_store_F_healthy{3},P_LV_store_F_healthy{3},'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(V_LV_store_F_healthy{4},P_LV_store_F_healthy{4},'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(V_LV_store_F_healthy{5},P_LV_store_F_healthy{5},'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(V_LV_store_F_healthy{6},P_LV_store_F_healthy{6},'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(V_LV_store_F_healthy{7},P_LV_store_F_healthy{7},'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(V_LV_store_F_healthy{8},P_LV_store_F_healthy{8},'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(V_LV_store_F_failure{1},P_LV_store_F_failure{1},':','linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(V_LV_store_F_failure{2},P_LV_store_F_failure{2},'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(V_LV_store_F_failure{3},P_LV_store_F_failure{3},'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(V_LV_store_F_failure{4},P_LV_store_F_failure{4},'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(V_LV_store_F_failure{5},P_LV_store_F_failure{5},'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(V_LV_store_F_failure{6},P_LV_store_F_failure{6},'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(V_LV_store_F_failure{7},P_LV_store_F_failure{7},'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(V_LV_store_F_failure{8},P_LV_store_F_failure{8},'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])

xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','1% dATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','HF','HF + 1% dATP','HF + 2% dATP','HF + 4% dATP','HF + 6% dATP','HF + 8% dATP','HF + 10% dATP','HF + 20% dATP')
set(gca, 'fontsize', 14)


percent_dATP = [1 2 4 6 8 10 20];
figure
subplot(4,3,1)
hold on
plot(0,CO_healthy{1}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),CO_healthy{2}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),CO_healthy{3}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),CO_healthy{4}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),CO_healthy{5}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),CO_healthy{6}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),CO_healthy{7}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),CO_healthy{8}/344.75,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,CO_failure{1}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),CO_failure{2}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),CO_failure{3}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),CO_failure{4}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),CO_failure{5}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),CO_failure{6}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),CO_failure{7}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),CO_failure{8}/344.75,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('SV (mL/beat)')
set(gca, 'fontsize', 16)

subplot(4,3,2)
hold on
plot(0,EF_healthy{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),EF_healthy{2}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),EF_healthy{3}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),EF_healthy{4}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),EF_healthy{5}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),EF_healthy{6}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),EF_healthy{7}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),EF_healthy{8}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,EF_failure{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),EF_failure{2}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),EF_failure{3}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),EF_failure{4}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),EF_failure{5}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),EF_failure{6}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),EF_failure{7}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),EF_failure{8}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('EF (%)')
set(gca, 'fontsize', 16)

subplot(4,3,3)
hold on
plot(0,CO_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),CO_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),CO_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),CO_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),CO_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),CO_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),CO_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),CO_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,CO_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),CO_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),CO_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),CO_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),CO_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),CO_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),CO_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),CO_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('CO (mL/min)')
set(gca, 'fontsize', 16)

subplot(4,3,4)
hold on
plot(0,LVDP_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),LVDP_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),LVDP_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),LVDP_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),LVDP_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),LVDP_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),LVDP_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),LVDP_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])

plot(0,LVDP_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),LVDP_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),LVDP_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),LVDP_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),LVDP_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),LVDP_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),LVDP_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),LVDP_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('LVDevP (mmHg)')
set(gca, 'fontsize', 16)

subplot(4,3,5)
hold on
plot(0,work_rate_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),work_rate_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),work_rate_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),work_rate_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),work_rate_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),work_rate_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),work_rate_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),work_rate_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,work_rate_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),work_rate_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),work_rate_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),work_rate_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),work_rate_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),work_rate_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),work_rate_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),work_rate_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('Work Rate (mmHg*mL/s)')
set(gca, 'fontsize', 16)

subplot(4,3,6)
hold on
plot(0,ATP_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),ATP_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),ATP_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),ATP_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),ATP_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),ATP_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),ATP_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),ATP_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,ATP_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),ATP_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),ATP_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),ATP_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),ATP_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),ATP_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),ATP_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),ATP_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ATP (mM)')
set(gca, 'fontsize', 16)

subplot(4,3,7)
hold on
plot(0,ADP_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),ADP_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),ADP_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),ADP_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),ADP_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),ADP_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),ADP_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),ADP_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,ADP_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),ADP_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),ADP_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),ADP_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),ADP_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),ADP_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),ADP_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),ADP_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ADP (\muM)')
set(gca, 'fontsize', 16)

subplot(4,3,8)
hold on
plot(0,Pi_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),Pi_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),Pi_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),Pi_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),Pi_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),Pi_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),Pi_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),Pi_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,Pi_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),Pi_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),Pi_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),Pi_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),Pi_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),Pi_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),Pi_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),Pi_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*7 0.2667+0.0888*6 0.6000-0.0222*7])
xlabel('Percent dATP')
ylabel('Pi (mM)')
set(gca, 'fontsize', 16)

subplot(4,3,9)
hold on
plot(0,MVO2_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),MVO2_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),MVO2_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),MVO2_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),MVO2_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),MVO2_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),MVO2_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),MVO2_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*7 0.1333+0.0889*6 0.5333+0.0111*7])

plot(0,MVO2_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),MVO2_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),MVO2_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),MVO2_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),MVO2_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),MVO2_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),MVO2_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),MVO2_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('MVO_{2} (uM O_{2}/min/g tissue)')
set(gca, 'fontsize', 16)

subplot(4,3,10)
hold on
plot(0,PCrATP_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),PCrATP_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),PCrATP_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),PCrATP_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),PCrATP_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),PCrATP_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),PCrATP_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),PCrATP_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,PCrATP_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),PCrATP_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),PCrATP_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),PCrATP_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),PCrATP_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),PCrATP_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),PCrATP_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),PCrATP_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('CrP/ATP ratio')
set(gca, 'fontsize', 16)

subplot(4,3,11)
hold on
plot(0,ATPase_F_healthy{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),ATPase_F_healthy{2},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),ATPase_F_healthy{3},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),ATPase_F_healthy{4},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),ATPase_F_healthy{5},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),ATPase_F_healthy{6},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),ATPase_F_healthy{7},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),ATPase_F_healthy{8},'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,ATPase_F_failure{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),ATPase_F_failure{2},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),ATPase_F_failure{3},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),ATPase_F_failure{4},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),ATPase_F_failure{5},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),ATPase_F_failure{6},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),ATPase_F_failure{7},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),ATPase_F_failure{8},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('ATPase (M/s/L)')
set(gca, 'fontsize', 16)

subplot(4,3,12)
hold on
plot(0,efficiency_F_healthy{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*0 0.1333+0.0767*0 0.5333+0.0095*0])
plot(percent_dATP(1),efficiency_F_healthy{2}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*1 0.1333+0.0767*1 0.5333+0.0095*1])
plot(percent_dATP(2),efficiency_F_healthy{3}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*2 0.1333+0.0767*2 0.5333+0.0095*2])
plot(percent_dATP(3),efficiency_F_healthy{4}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*3 0.1333+0.0767*3 0.5333+0.0095*3])
plot(percent_dATP(4),efficiency_F_healthy{5}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*4 0.1333+0.0767*4 0.5333+0.0095*4])
plot(percent_dATP(5),efficiency_F_healthy{6}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*5 0.1333+0.0767*5 0.5333+0.0095*5])
plot(percent_dATP(6),efficiency_F_healthy{7}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*6 0.1333+0.0767*6 0.5333+0.0095*6])
plot(percent_dATP(7),efficiency_F_healthy{8}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.01*7 0.1333+0.0767*7 0.5333+0.0095*7])

plot(0,efficiency_F_failure{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*0 0.2667+0.0762*0 0.6000-0.0190*0])
plot(percent_dATP(1),efficiency_F_failure{2}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*1 0.2667+0.0762*1 0.6000-0.0190*1])
plot(percent_dATP(2),efficiency_F_failure{3}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*2 0.2667+0.0762*2 0.6000-0.0190*2])
plot(percent_dATP(3),efficiency_F_failure{4}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*3 0.2667+0.0762*3 0.6000-0.0190*3])
plot(percent_dATP(4),efficiency_F_failure{5}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*4 0.2667+0.0762*4 0.6000-0.0190*4])
plot(percent_dATP(5),efficiency_F_failure{6}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*5 0.2667+0.0762*5 0.6000-0.0190*5])
plot(percent_dATP(6),efficiency_F_failure{7}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*6 0.2667+0.0762*6 0.6000-0.0190*6])
plot(percent_dATP(7),efficiency_F_failure{8}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6000-0.0190*7])
xlabel('Percent dATP')
ylabel('Efficiency (mL^{2}*mmHg*ms/M)')
set(gca, 'fontsize', 16)

