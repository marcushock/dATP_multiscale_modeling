%%%% Plots myocyte level figures from:
% Teitgen, A. E., Hock, M. T., McCabe, K. J., Childers, M. C., Huber, G. A.,
% Marzban, B., Beard, D. A., McCammon, J. A., Regnier, M., McCulloch, A. D. (2024).
% Multiscale modeling shows how 2'-deoxy-ATP rescues ventricular function
% in heart failure. PNAS.

% Calcium transient parameters
a_ATP = 0.1208;
b_ATP = 0.6651;
c_ATP = 1.7374;

a_dATP = 0.0864;
b_dATP = 0.5971;
c_dATP = 2.0301;

%% Digitized experimental data (interpolating at specified points)
% ***Run this section before any of other figure sections
tspan = 0:0.001:2;

% Force pCa (Regnier et al. 2004)
Ca_fraction_fpca = [0.1:0.1:100];
ATP_points = csvread('ATP_points.csv');
dATP_points = csvread('dATP_points.csv');
ATP_Hill = csvread('ATP_Hill.csv');
ATP_Hill1 = ATP_Hill(:,1);
ATP_Hill2 = ATP_Hill(:,2);
dATP_Hill = csvread('dATP_Hill.csv');
dATP_Hill1 = dATP_Hill(:,1);
dATP_Hill2 = dATP_Hill(:,2);
[x, index] = unique(ATP_Hill1); 
ATP_Hill_interp = interp1(x, ATP_Hill2(index), Ca_fraction_fpca);
ATP_Hill_interp(1:40) = 0.9818;
ATP_Hill_interp(58:end) = 0.0061;
[x, index] = unique(dATP_Hill1); 
dATP_Hill_interp = interp1(x, dATP_Hill2(index), Ca_fraction_fpca);

% Model inputs:
% Ca_value, XB, dATP_percent, kf_i, k_f_i, kw_i, k_w_i, kp_i, k_p_i, kg_i, krecruit_i, k_on_i, k_off_i, k_coop_i

last_beat = 2000;

%% Figure 4
% ATP
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~] = myocyte_model(0, 1, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
force_final_last = force_final_ATP(idx_XB_ATP:end);
T_XB_final_last = T_final_XB_ATP(idx_XB_ATP:end);
Shortening_final_last = Shortening_final_ATP(idx_XB_ATP:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(1) = FS;
TTP_store(1) = TTP_shortening;
RT50_store(1) = time_RT50_shortening;
RT90_store(1) = time_RT90_shortening;

% dATP kf+
[T_final_XB_dATP_kf, force_final_dATP_kf, idx_XB_dATP_kf, Shortening_final_dATP_kf, ~] = myocyte_model(0, 1, 1, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
force_final_last = force_final_dATP_kf(idx_XB_dATP_kf:end);
T_XB_final_last = T_final_XB_dATP_kf(idx_XB_dATP_kf:end);
Shortening_final_last = Shortening_final_dATP_kf(idx_XB_dATP_kf:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(2) = FS;
TTP_store(2) = TTP_shortening;
RT50_store(2) = time_RT50_shortening;
RT90_store(2) = time_RT90_shortening;

% dATP kf+ + kf- + kw+
[T_final_XB_dATP_kf_k_f_kw, force_final_dATP_kf_k_f_kw, idx_XB_dATP_kf_k_f_kw, Shortening_final_dATP_kf_k_f_kw, ~] = myocyte_model(0, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
force_final_last = force_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end);
T_XB_final_last = T_final_XB_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end);
Shortening_final_last = Shortening_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(3) = FS;
TTP_store(3) = TTP_shortening;
RT50_store(3) = time_RT50_shortening;
RT90_store(3) = time_RT90_shortening;

% dATP krecruit
[T_final_XB_dATP_krecruit, force_final_dATP_krecruit, idx_XB_dATP_krecruit, Shortening_final_dATP_krecruit, ~] = myocyte_model(0, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
force_final_last = force_final_dATP_krecruit(idx_XB_dATP_krecruit:end);
T_XB_final_last = T_final_XB_dATP_krecruit(idx_XB_dATP_krecruit:end);
Shortening_final_last = Shortening_final_dATP_krecruit(idx_XB_dATP_krecruit:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(4) = FS;
TTP_store(4) = TTP_shortening;
RT50_store(4) = time_RT50_shortening;
RT90_store(4) = time_RT90_shortening;

% dATP kf+ + kf- + kw+ + krecruit
[T_final_XB_dATP_kf_k_f_kw_krecruit, force_final_dATP_kf_k_f_kw_krecruit, idx_XB_dATP_kf_k_f_kw_krecruit, Shortening_final_dATP_kf_k_f_kw_krecruit, ~] = myocyte_model(0, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
force_final_last = force_final_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end);
T_XB_final_last = T_final_XB_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end);
Shortening_final_last = Shortening_final_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(5) = FS;
TTP_store(5) = TTP_shortening;
RT50_store(5) = time_RT50_shortening;
RT90_store(5) = time_RT90_shortening;

% dATP Ca
[T_final_XB_dATP_Ca, force_final_dATP_Ca, idx_XB_dATP_Ca, Shortening_final_dATP_Ca, ~] = myocyte_model(1, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP); 
force_final_last = force_final_dATP_Ca(idx_XB_dATP_Ca:end);
T_XB_final_last = T_final_XB_dATP_Ca(idx_XB_dATP_Ca:end);
Shortening_final_last = Shortening_final_dATP_Ca(idx_XB_dATP_Ca:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(6) = FS;
TTP_store(6) = TTP_shortening;
RT50_store(6) = time_RT50_shortening;
RT90_store(6) = time_RT90_shortening;

% dATP krecruit + Ca
[T_final_XB_dATP_krecruit_Ca, force_final_dATP_krecruit_Ca, idx_XB_dATP_krecruit_Ca, Shortening_final_dATP_krecruit_Ca, ~] = myocyte_model(1, 1, 1, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 37, 50, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP); 
force_final_last = force_final_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end);
T_XB_final_last = T_final_XB_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end);
Shortening_final_last = Shortening_final_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(7) = FS;
TTP_store(7) = TTP_shortening;
RT50_store(7) = time_RT50_shortening;
RT90_store(7) = time_RT90_shortening;

% dATP kf+ + kf- + krecruit + Ca
[T_final_XB_dATP_kf_k_f_kw_krecruit_Ca, force_final_dATP_kf_k_f_kw_krecruit_Ca, idx_XB_dATP_kf_k_f_kw_krecruit_Ca, Shortening_final_dATP_kf_k_f_kw_krecruit_Ca, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 50, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP); 
force_final_last = force_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end);
T_XB_final_last = T_final_XB_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end);
Shortening_final_last = Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(8) = FS;
TTP_store(8) = TTP_shortening;
RT50_store(8) = time_RT50_shortening;
RT90_store(8) = time_RT90_shortening;

% From experimental data (% change)
FS_store(9) = 0.34; %0.55; %0.26; 
TTP_store(9) = 0; %-0.03; %-0.07; 
RT50_store(9) = -0.30; %-0.29; %-0.10; 
RT90_store(9) = -0.28; %-0.43; %-0.16; 

figure
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat/1000,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)/1000-last_beat/1000,Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)./max(Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative shortening')
legend('ATP','dATP')
set(gca,'FontSize',14)


x_axis = [0 1 2 3 4 5 6];

figure
subplot(3,1,1)
hold on
bar(x_axis,(FS_store(2:8)-FS_store(1))./FS_store(1),'facecolor',[0.5, 0.5, 0.5])
yline((FS_store(1)-FS_store(1))./FS_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(FS_store(9),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('FS')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

subplot(3,1,2)
hold on
bar(x_axis,(RT50_store(2:8)-RT50_store(1))./RT50_store(1),'facecolor',[0.5, 0.5, 0.5])
yline((RT50_store(1)-RT50_store(1))./RT50_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(RT50_store(9),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('RT50')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

subplot(3,1,3)
hold on
bar(x_axis,(RT90_store(2:8)-RT90_store(1))./RT90_store(1),'facecolor',[0.5, 0.5, 0.5])
yline((RT90_store(1)-RT90_store(1))./RT90_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(RT90_store(9),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('RT90')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

%% Figure S12
[~, ~, ~, ~, SS_Ftotal_fpca_ATP] = myocyte_model(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_kf] = myocyte_model(0, 0, 100, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_kf_k_f_kw] = myocyte_model(0, 0, 100, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_kf./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.5])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative force');
legend('ATP model','dATP model','ATP data','dATP data')

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_kf_k_f_kw./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.5])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative force');
legend('ATP model','dATP model','ATP data','dATP data')

%% Figure S13
[~, ~, ~, ~, SS_Ftotal_fpca_ATP] = myocyte_model(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % ATP force pCa

% kf+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_10] = myocyte_model(0, 0, 100, 250*0.1, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf+ 10%
[~, ~, ~, ~, SS_Ftotal_fpca_kf_50] = myocyte_model(0, 0, 100, 250*0.5, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf+ 50%
[~, ~, ~, ~, SS_Ftotal_fpca_kf_150] = myocyte_model(0, 0, 100, 250*1.5, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf+ 150%
[~, ~, ~, ~, SS_Ftotal_fpca_kf_200] = myocyte_model(0, 0, 100, 250*2, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf+ 200%

% kf-
[~, ~, ~, ~, SS_Ftotal_fpca_k_f_10] = myocyte_model(0, 0, 100, 250, 304.6708*0.1, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf- 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k_f_50] = myocyte_model(0, 0, 100, 250, 304.6708*0.5, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf- 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k_f_150] = myocyte_model(0, 0, 100, 250, 304.6708*1.5, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf- 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k_f_200] = myocyte_model(0, 0, 100, 250, 304.6708*2, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kf- 200%

% kw+
[~, ~, ~, ~, SS_Ftotal_fpca_kw_10] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727*0.1, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw+ 10%
[~, ~, ~, ~, SS_Ftotal_fpca_kw_50] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727*0.5, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw+ 50%
[~, ~, ~, ~, SS_Ftotal_fpca_kw_150] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727*1.5, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw+ 150%
[~, ~, ~, ~, SS_Ftotal_fpca_kw_200] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727*2, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw+ 200%

% kw-
[~, ~, ~, ~, SS_Ftotal_fpca_k_w_10] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296*0.1, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw- 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k_w_50] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296*0.5, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw- 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k_w_150] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296*1.5, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw- 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k_w_200] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296*2, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kw- 200%

% kp+
[~, ~, ~, ~, SS_Ftotal_fpca_kp_10] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72*0.1, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp+ 10%
[~, ~, ~, ~, SS_Ftotal_fpca_kp_50] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72*0.5, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp+ 50%
[~, ~, ~, ~, SS_Ftotal_fpca_kp_150] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72*1.5, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp+ 150%
[~, ~, ~, ~, SS_Ftotal_fpca_kp_200] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72*2, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp+ 200%

% kp-
[~, ~, ~, ~, SS_Ftotal_fpca_k_p_10] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25*0.1, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp- 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k_p_50] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25*0.5, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp- 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k_p_150] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25*1.5, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp- 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k_p_200] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25*2, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kp- 200%

% kg+
[~, ~, ~, ~, SS_Ftotal_fpca_kg_10] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586*0.1, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kg+ 10%
[~, ~, ~, ~, SS_Ftotal_fpca_kg_50] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586*0.5, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kg+ 50%
[~, ~, ~, ~, SS_Ftotal_fpca_kg_150] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586*1.5, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kg+ 150%
[~, ~, ~, ~, SS_Ftotal_fpca_kg_200] = myocyte_model(0, 0, 100, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586*2, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % force pCa scaling kg+ 200%


figure
subplot(2,4,1)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_f^+')

subplot(2,4,2)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_f_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_f_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_f_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_f_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_f^-')

subplot(2,4,3)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kw_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kw_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kw_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kw_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_w^+')

subplot(2,4,4)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_w_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_w_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_w_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_w_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_w^-')

subplot(2,4,5)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kp_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kp_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kp_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kp_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_p^+')

subplot(2,4,6)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_p_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_p_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_p_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_p_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_p^-')

subplot(2,4,7)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kg_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kg_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kg_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kg_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_g^+')

%% Figure S14
% ATP
[~, ~, ~, ~, SS_Ftotal_fpca_ATP] = myocyte_model(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % ATP force pCa

% kf+
[~, ~, ~, ~, SS_Ftotal_fpca_kf] = myocyte_model(0, 0, 100, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

% kf+ + kf-
[~, ~, ~, ~, SS_Ftotal_fpca_kf_k_f] = myocyte_model(0, 0, 100, 478, 460, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

% kf+ + kw+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_kw] = myocyte_model(0, 0, 100, 478, 304.6708, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

% kf+ + kg+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_kg] = myocyte_model(0, 0, 100, 478, 304.6708, 112.3727, 21.296, 811.72, 43.25, 180, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

% kg+ + kf- + kw+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_k_f_kw] = myocyte_model(0, 0, 100, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% kf+ + kw+ + kg+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_kw_kg] = myocyte_model(0, 0, 100, 478, 304.6708, 170, 21.296, 811.72, 43.25, 180, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 

% kf+ + kf- + kg+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_k_f_kg] = myocyte_model(0, 0, 100, 478, 460, 112.3727, 21.296, 811.72, 43.25, 150, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

% kf+ + kf- + kw+ + kg+
[~, ~, ~, ~, SS_Ftotal_fpca_kf_k_f_kw_kg] = myocyte_model(0, 0, 100, 478, 460, 170, 21.296, 811.72, 43.25, 150, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP);

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf/max(SS_Ftotal_fpca_ATP),':','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_k_f/max(SS_Ftotal_fpca_ATP),'-.','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_kw/max(SS_Ftotal_fpca_ATP),'--','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_kg/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.5, 0.5, 0.5])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','dATP k_f^+','dATP k_f^+ + k_f^-','dATP k_f^+ + k_w^+','dATP k_f^+ + k_g^+','ATP data','dATP data')

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf/max(SS_Ftotal_fpca_ATP),':','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_k_f_kw/max(SS_Ftotal_fpca_ATP),'-.','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_kw_kg/max(SS_Ftotal_fpca_ATP),'--','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_k_f_kg/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kf_k_f_kw_kg/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.8, 0.8, 0.8])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','dATP k_f^+','dATP k_f^+ + k_f^- + k_w^+','dATP k_f^+ + k_w+ + k_g^+','dATP k_f^+ + k_f^- + k_g^+','dATP k_f^+ + k_f^- + k_w^+ + k_g^+','ATP data','dATP data')

%% Figure S15 (run simulations for Figure 4 to generate)
figure
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat/1000,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_kf(idx_XB_dATP_kf:end)/1000-last_beat/1000,Shortening_final_dATP_kf(idx_XB_dATP_kf:end)./max(Shortening_final_dATP_kf(idx_XB_dATP_kf:end)),'linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)/1000-last_beat/1000,Shortening_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)./max(Shortening_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)),':','linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_krecruit(idx_XB_dATP_krecruit:end)/1000-last_beat/1000,Shortening_final_dATP_krecruit(idx_XB_dATP_krecruit:end)./max(Shortening_final_dATP_krecruit(idx_XB_dATP_krecruit:end)),'-.','linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_Ca(idx_XB_dATP_Ca:end)/1000-last_beat/1000,Shortening_final_dATP_Ca(idx_XB_dATP_Ca:end)./max(Shortening_final_dATP_Ca(idx_XB_dATP_Ca:end)),'linewidth',3,'color',[0.8 0.8 0.8])
plot(T_final_XB_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end)/1000-last_beat/1000,Shortening_final_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end)./max(Shortening_final_dATP_krecruit_Ca(idx_XB_dATP_krecruit_Ca:end)),'--','linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)/1000-last_beat/1000,Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)./max(Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative shortening')
legend('ATP','dATP k_f', 'dATP k_f^+ + k_f^- + k_w^+','dATP k_{recruit}','dATP Ca','dATP k_{recruit} + Ca','dATP k_f^+ + k_f^- + k_w^+ + k_{recruit} + Ca')
set(gca,'FontSize',14)

%% Figure S16
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~] = myocyte_model(0, 1, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); 
Max_shortening_ATP = min(Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)));

k_recruit_range = [0.2069 0.3 0.4 0.5 0.6 0.7 0.8 0.8 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 60 80 100 200 500 1000];
for i = 1:length(k_recruit_range)
    [T_final_XB_dATP_kf_k_f_kw_krecruit_Ca, force_final_dATP_kf_k_f_kw_krecruit_Ca, idx_XB_dATP_kf_k_f_kw_krecruit_Ca, Shortening_final_dATP_kf_k_f_kw_krecruit_Ca, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, k_recruit_range(i), 50, 723.8520, 9.6846, a_dATP, b_dATP, c_dATP); 
    max_shortening_store(i) = min(Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)./max(Shortening_final_dATP_kf_k_f_kw_krecruit_Ca(idx_XB_dATP_kf_k_f_kw_krecruit_Ca:end)));
end

figure
hold on
plot(k_recruit_range, max_shortening_store,'o')
yline(0.905,'--','color',[0.8 0.8 0.8])
xlabel('k_{recruit}')
ylabel('Relative Max Shortening')


% Shortening
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~] = myocyte_model(0, 1, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % ATP shortening
[T_final_XB_dATP_kf_k_f_kw, force_final_dATP_kf_k_f_kw, idx_XB_dATP_kf_k_f_kw, Shortening_final_dATP_kf_k_f_kw, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP kf+ kf- kw+ shortening
[T_final_XB_dATP_kf_k_f_kw_kp, force_final_dATP_kf_k_f_kw_kp, idx_XB_dATP_kf_k_f_kw_kp, Shortening_final_dATP_kf_k_f_kw_kp, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 100000, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP kf+ kf- kw+ kp+ shortening
[T_final_XB_dATP_kf_k_f_kw_kg, force_final_dATP_kf_k_f_kw_kg, idx_XB_dATP_kf_k_f_kw_kg, Shortening_final_dATP_kf_k_f_kw_kg, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 0.1, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP kf+ kf- kw+ kg+ shortening
[T_final_XB_dATP_kf_k_f_kw_krecruit, force_final_dATP_kf_k_f_kw_krecruit, idx_XB_dATP_kf_k_f_kw_krecruit, Shortening_final_dATP_kf_k_f_kw_krecruit, ~] = myocyte_model(1, 1, 1, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 37, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP kf+ kf- kw+ krecruit shortening

last_beat = 2;

figure
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)/1000-last_beat,Shortening_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)./max(Shortening_final_dATP_kf_k_f_kw(idx_XB_dATP_kf_k_f_kw:end)),':','linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_kf_k_f_kw_kp(idx_XB_dATP_kf_k_f_kw_kp:end)/1000-last_beat,Shortening_final_dATP_kf_k_f_kw_kp(idx_XB_dATP_kf_k_f_kw_kp:end)./max(Shortening_final_dATP_kf_k_f_kw_kp(idx_XB_dATP_kf_k_f_kw_kp:end)),'-.','linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_kf_k_f_kw_kg(idx_XB_dATP_kf_k_f_kw_kg:end)/1000-last_beat,Shortening_final_dATP_kf_k_f_kw_kg(idx_XB_dATP_kf_k_f_kw_kg:end)./max(Shortening_final_dATP_kf_k_f_kw_kg(idx_XB_dATP_kf_k_f_kw_kg:end)),'linewidth',3,'color',[0.5 0.5 0.5])
plot(T_final_XB_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end)/1000-last_beat,Shortening_final_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end)./max(Shortening_final_dATP_kf_k_f_kw_krecruit(idx_XB_dATP_kf_k_f_kw_krecruit:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP','dATP, k_f^+ + k_f^- + k_w^+','dATP, k_f^+ + k_f^- + k_w^+, increasing k_p^+', 'dATP, k_f^+ + k_f^- + k_w^+, increasing k_g^+', 'dATP, k_f^+ + k_f^- + k_w^+, increasing k_{recruit}')
ylim([0.87 1])
set(gca,'FontSize',14)
ylim([0.9 1])

%% Figure S18
[~, ~, ~, ~, SS_Ftotal_fpca_ATP] = myocyte_model(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % ATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP] = myocyte_model(0, 0, 100, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 50, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_ATP_shifted] = myocyte_model(0, 0, 0, 250, 304.6708, 112.3727, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % ATP force pCa shifted
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_shifted] = myocyte_model(0, 0, 100, 478, 460, 170, 21.296, 811.72, 43.25, 144.5586, 0.2069, 200, 723.8520, 9.6846, a_ATP, b_ATP, c_ATP); % dATP force pCa shifted

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP_shifted/max(SS_Ftotal_fpca_ATP_shifted),'--','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_shifted./max(SS_Ftotal_fpca_ATP_shifted),'--','linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model in vitro','dATP model in vitro','ATP model in vivo','dATP model in vivo','ATP data','dATP data')
