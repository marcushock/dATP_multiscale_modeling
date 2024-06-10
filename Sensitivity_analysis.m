%% Myocyte
% k_recruit_range = [0.2 1 10 20 30 37 40 50 100 200 300 400 500 600 700 800];
% a_range = [0.1208    0.1160    0.1114    0.1068    0.1022    0.0972    0.0917    0.0868    0.0864     0.0821    0.0773    0.0736    0.0689    0.0636];
% b_range = [0.6651    0.6649    0.6669    0.6696    0.6707    0.6680    0.6621    0.6517    0.5971     0.6635    0.6526    0.6870    0.6846    0.6763];
% c_range = [1.7374   1.7786     1.8201    1.8634    1.9084    1.9571    2.0116    2.0601    2.0301     2.1234    2.1772    2.2444    2.3104    2.3858];
% 
% 
% last_beat = 2000;
% 
% for i = 1:length(k_recruit_range)
%     k_recruit = k_recruit_range(i);
%     for j = 1:length(a_range)
%     a = a_range(j);
%     b = b_range(j);
%     c = c_range(j);
%     
%     [T_final_XB, force_final, idx_XB, Shortening_final, ~] = myocyte_model(0, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, k_recruit, 50, 723.8520, 9.6846, a, b, c);
% 
%     force_final_last = force_final(idx_XB:end);
%     T_XB_final_last = T_final_XB(idx_XB:end);
%     Shortening_final_last = Shortening_final(idx_XB:end);
%     peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
%     FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
%     [min_SL,idx_min] = min(Shortening_final_last);
%     TTP_shortening = T_XB_final_last(idx_min) - last_beat;
% 
%     RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
%     time_RT50_shortening = 0;
%     for p = idx_min:length(Shortening_final_last)
%         if Shortening_final_last(p) > RT50
%             time_RT50_shortening = ((T_XB_final_last(p)+T_XB_final_last(p-1))/2);
%             time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
%             break
%         end
%     end
% 
%     RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
%     time_RT90_shortening = 0;
%     for p = idx_min:length(Shortening_final_last)
%         if Shortening_final_last(p) > RT90
%             time_RT90_shortening = ((T_XB_final_last(p)+T_XB_final_last(p-1))/2);
%             time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
%             break
%         end
%     end
%     
%     FS_SA{i,j} = FS;
%     RT50_SA{i,j} = time_RT50_shortening;
%     RT90_SA{i,j} = time_RT90_shortening;
%     
%     i
%     j
%     end
% end
% 
% 
% % FS
% figure
% FS_new = cell2mat(FS_SA);
% columnlabels = {'0.2 (ATP)','1','10','20','30','37 (dATP)','40','50','100','200','300','400','500','600','700','800'};
% rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','25,31 (dATP)','28%','32%','36%','40%','44%','48%'};h = heatmap(rowlabels, columnlabels, FS_new);
% h.Colormap = parula;
% h.XLabel = 'Decrease in RT50 and RT90';
% h.YLabel = 'k_{recruit}';
% 
% % RT50
% figure
% RT50_new = cell2mat(RT50_SA);
% columnlabels = {'0.2 (ATP)','1','10','20','30','37 (dATP)','40','50','100','200','300','400','500','600','700','800'};
% rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','25,31 (dATP)','28%','32%','36%','40%','44%','48%'};h = heatmap(rowlabels, columnlabels, RT50_new);
% h.Colormap = parula;
% h.XLabel = 'Decrease in RT50 and RT90';
% h.YLabel = 'k_{recruit}';
% 
% % RT90
% figure
% RT90_new = cell2mat(RT90_SA);
% columnlabels = {'0.2 (ATP)','1','10','20','30','37 (dATP)','40','50','100','200','300','400','500','600','700','800'};
% rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','25,31 (dATP)','28%','32%','36%','40%','44%','48%'};
% h = heatmap(rowlabels, columnlabels, RT90_new);
% h.Colormap = parula;
% h.XLabel = 'Decrease in RT50 and RT90';
% h.YLabel = 'k_{recruit}';

%% Ventricular

%k_recruit_range = 37; %[2, 4, 6, 8, 10]; %[0.2, 0.4, 0.6, 0.8, 1];
%2, 4, 6, 8, 10, 12];
%12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34
     %36, 38, 40, 42, 44, 46, 48, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800];

k_recruit_range = [0.2 1 10 20 30 37 40 50 100 200 300 400 500 600 700 800]; %[0.2 1 10 20 30 40 50 100 200 300 400 500 600 700 800];
%[0.2 1 10 20 30 37 40 50 100 200 300 400 500 600 700 800];

% Ca range
a_range = 0.2218; %[0.2968 0.2920 0.2821 0.2759 0.2664 0.2502 0.2409 0.2204 0.2113 0.1972];


    %[0.2968    0.2977    0.2920    0.2876    0.2821    0.2804    0.2759    0.2720    0.2664    0.2600    0.2502    0.2447    0.2409    0.2261...
    %0.2204    0.2108    0.2113    0.2049    0.1972    0.1926    0.1801    0.1721    0.1733    0.1774    0.1614    0.1600];

b_range = 1.4980; %[1.5900 1.6557 1.6472 1.7059 1.8181 1.6987 1.8327 1.6532 1.7356 1.7537];


    %[1.5900    1.6178    1.6557    1.6504    1.6472    1.6997    1.7059    1.7688    1.8181    1.7989    1.6987    1.7681    1.8327    1.5972...
    %1.6532    1.5730    1.7356    1.7741    1.7537    1.8401    1.7001    1.6373    1.9017    2.3100    1.9127    1.9867];

c_range = 1.9989; %[1.7171 1.7400 1.7733 1.8013 1.8448 1.8985 1.9464 2.0207 2.0714 2.1416];


    %[1.7171    1.7167    1.7400    1.7546    1.7733    1.7841    1.8013    1.8203    1.8448    1.8682    1.8985    1.9259    1.9464    1.9903...
    %2.0207    2.0576    2.0714    2.1037    2.1416    2.1709    2.2266    2.2670    2.2821    2.2818    2.3533    2.3677];

for i = 1:length(k_recruit_range)
    k_recruit = k_recruit_range(i);
    for j = 1:length(a_range)
    a = a_range(j);
    b = b_range(j);
    c = c_range(j);

    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F] = CardiovascularMechanics(1, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, k_recruit, 200, 723.8520, 9.6846, a, b, c);  
    V_LV_store_F_SA_31{i,j} = V_LV_store_F;
    P_LV_store_F_SA_31{i,j} = P_LV_store_F;
    max_force_SA_31{i,j} = max_force;
    FS_SA_31{i,j} = FS;
    CO_SA_31{i,j} = CO;
    EF_SA_31{i,j} = EF;
    LVDP_SA_31{i,j} = LVDP;
    work_rate_F_SA_31{i,j} = work_rate_F;
    ATP_F_SA_31{i,j} = ATP_F; 
    ADP_F_SA_31{i,j} = ADP_F;
    Pi_F_SA_31{i,j} = Pi_F;
    MVO2_F_SA_31{i,j} = MVO2_F;
    PCrATP_F_SA_31{i,j} = PCrATP_F;
    XB_turnover_F_SA_31{i,j} = XB_turnover_F;
    ATPase_F_SA_31{i,j} = ATPase_F;
    efficiency_F_SA_31{i,j} = efficiency_F;
    i
    j
    end
end




%EF
EF_new = cell2mat(EF_SA);
EF_new_37 = cell2mat(EF_SA_37);
EF_new_31 = cell2mat(EF_SA_31);
EF_full = [EF_new(1:5,:); EF_new_37; EF_new(6:end,:)];
EF_full = [EF_full(:,1:7), EF_new_31, EF_full(:,8:end)];

figure
%EF_new = cell2mat(EF_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, EF_full*100);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%CO
CO_new = cell2mat(CO_SA);
CO_new_37 = cell2mat(CO_SA_37);
CO_new_31 = cell2mat(CO_SA_31);
CO_full = [CO_new(1:5,:); CO_new_37; CO_new(6:end,:)];
CO_full = [CO_full(:,1:7), CO_new_31, CO_full(:,8:end)];

figure
%CO_new = cell2mat(CO_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, CO_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%LVDevP
LVDP_new = cell2mat(LVDP_SA);
LVDP_new_37 = cell2mat(LVDP_SA_37);
LVDP_new_31 = cell2mat(LVDP_SA_31);
LVDP_full = [LVDP_new(1:5,:); LVDP_new_37; LVDP_new(6:end,:)];
LVDP_full = [LVDP_full(:,1:7), LVDP_new_31, LVDP_full(:,8:end)];

figure
%LVDP_new = cell2mat(LVDP_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, LVDP_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%Work rate
work_rate_F_new = cell2mat(work_rate_F_SA);
work_rate_F_new_37 = cell2mat(work_rate_F_SA_37);
work_rate_F_new_31 = cell2mat(work_rate_F_SA_31);
work_rate_F_full = [work_rate_F_new(1:5,:); work_rate_F_new_37; work_rate_F_new(6:end,:)];
work_rate_F_full = [work_rate_F_full(:,1:7), work_rate_F_new_31, work_rate_F_full(:,8:end)];

figure
%work_rate_F_new = cell2mat(work_rate_F_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, work_rate_F_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%MVO2
MVO2_F_new = cell2mat(MVO2_F_SA);
MVO2_F_new_37 = cell2mat(MVO2_F_SA_37);
MVO2_F_new_31 = cell2mat(MVO2_F_SA_31);
MVO2_F_full = [MVO2_F_new(1:5,:); MVO2_F_new_37; MVO2_F_new(6:end,:)];
MVO2_F_full = [MVO2_F_full(:,1:7), MVO2_F_new_31, MVO2_F_full(:,8:end)];

figure
%MVO2_F_new = cell2mat(MVO2_F_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, MVO2_F_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%CrPATP
PCrATP_F_new = cell2mat(PCrATP_F_SA);
PCrATP_F_new_37 = cell2mat(PCrATP_F_SA_37);
PCrATP_F_new_31 = cell2mat(PCrATP_F_SA_31);
PCrATP_F_full = [PCrATP_F_new(1:5,:); PCrATP_F_new_37; PCrATP_F_new(6:end,:)];
PCrATP_F_full = [PCrATP_F_full(:,1:7), PCrATP_F_new_31, PCrATP_F_full(:,8:end)];

figure
%PCrATP_F_new = cell2mat(PCrATP_F_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, PCrATP_F_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%ATPase
ATPase_F_new = cell2mat(ATPase_F_SA);
ATPase_F_new_37 = cell2mat(ATPase_F_SA_37);
ATPase_F_new_31 = cell2mat(ATPase_F_SA_31);
ATPase_F_full = [ATPase_F_new(1:5,:); ATPase_F_new_37; ATPase_F_new(6:end,:)];
ATPase_F_full = [ATPase_F_full(:,1:7), ATPase_F_new_31, ATPase_F_full(:,8:end)];

figure
%ATPase_F_new = cell2mat(ATPase_F_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, ATPase_F_full);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';

%Efficiency
efficiency_F_new = cell2mat(efficiency_F_SA);
efficiency_F_new_37 = cell2mat(efficiency_F_SA_37);
efficiency_F_new_31 = cell2mat(efficiency_F_SA_31);
efficiency_F_full = [efficiency_F_new(1:5,:); efficiency_F_new_37; efficiency_F_new(6:end,:)];
efficiency_F_full = [efficiency_F_full(:,1:7), efficiency_F_new_31, efficiency_F_full(:,8:end)];

figure
%efficiency_F_new = cell2mat(efficiency_F_SA);
columnlabels = {'0.2 (ATP)','1','10','20','30','37','40','50','100','200','300','400','500','600','700','800'};
rowlabels = {'0% (ATP)','4%','8%','12%','16%','20%','24%','dATP','28%','32%','36%'};
h = heatmap(rowlabels, columnlabels, efficiency_F_full/1000);
h.Colormap = parula;
h.XLabel = 'Decrease in RT50 and RT90';
h.YLabel = 'k_{recruit}';






