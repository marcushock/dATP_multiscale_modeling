%%%% Script used to fit calcium transient parameters to experimental data

stim_period = 1000; %174; %1000 %ms
t = linspace(0,stim_period,stim_period+1);

% Korte myocyte
% a = 0.106;
% b = 0.5635;
% c = 1.8017;
% Ca0 = 0;
% 
% a = 0.0534;
% b = 0.5484;
% c = 2.4732;
% Ca0 = 0;
% 
% Korte ventricular
% a = 0.2968;
% b = 1.5900;
% c = 1.7171;
% Ca0 = 0.1622;
% 
% a = 0.1447;
% b = 2.1011;
% c = 2.4775;
% Ca0 = 0.1622;

% Average myocyte
a = 0.1208;
b = 0.6651;
c = 1.7374;
Ca0 = 0;

a = 0.0864;
b = 0.5971;
c = 2.0301;
Ca0 = 0;

% Average ventricular
% a = 0.2968;
% b = 1.5900;
% c = 1.7171;
% Ca0 = 0.1622;
% 
% a = 0.2218;
% b = 1.4980;
% c = 1.9989;
% Ca0 = 0.1622;

% Myocyte
Ca0 = 0;
% Percent change
% 4
a_store(1) = 0.1160;    b_store(1) = 0.6649;    c_store(1) = 1.7786;
% 8
a_store(2) = 0.1114;    b_store(2) = 0.6669;    c_store(2) = 1.8201;
% 12
a_store(3) = 0.1068;    b_store(3) = 0.6696;    c_store(3) = 1.8634;
% 16
a_store(4) = 0.1022;    b_store(4) = 0.6707;    c_store(4) = 1.9084;
% 20
a_store(5) = 0.0972;    b_store(5) = 0.6680;    c_store(5) = 1.9571;
% 24
a_store(6) = 0.0917;    b_store(6) = 0.6621;    c_store(6) = 2.0116;
% 28
a_store(7) = 0.0868;    b_store(7) = 0.6517;    c_store(7) = 2.0601;
% 32
a_store(8) = 0.0821;    b_store(8) = 0.6635;    c_store(8) = 2.1234;
% 36
a_store(9) = 0.0773;    b_store(9) = 0.6526;    c_store(9) = 2.1772;
% 40
a_store(10) = 0.0736;    b_store(10) = 0.6870;    c_store(10) = 2.2444;
% 44
a_store(11) = 0.0689;    b_store(11) = 0.6846;    c_store(11) = 2.3104;
% 48
a_store(12) = 0.0636;    b_store(12) = 0.6763;    c_store(12) = 2.3858;


figure
hold on
h = legend('show','location','best');
percent_change = 0;
name = num2str(percent_change);
plot(t, Ca_i_ATP,'DisplayName',name,'linewidth',2,'color',[0.2,0.1333,0.5333])
plot(t, Ca_i_dATP,':','DisplayName',name,'linewidth',2,'color',[0.5 0.5 0.5])
% 0.26667, 0.6667, 0.6


for j = 1:12
    x_final = [a_store(j), b_store(j), c_store(j)];

    for i=1:stim_period+1
    phi = mod(t(i)+0.001,stim_period)/stim_period;
    Ca_i(i) = (x_final(1)/phi)*exp(-x_final(2)*(log(phi)+x_final(3))^2) + Ca0;
    end

    percent_change = j*4;
    name = num2str(percent_change);
    plot(t, Ca_i,'DisplayName',name,'linewidth',2,'color',[0.2+0.00555*j,0.1333+0.0444*j,0.5333+0.00555*j])
end




% 
% % 
% % 
% % % Ventricular
% % Ca0 = 0.1622;
% % % Percent change
% % % 2
% % a_store(1) = 0.2977;    b_store(1) = 1.6178;    c_store(1) = 1.7167;
% % % 4
% % a_store(2) = 0.2920;    b_store(2) = 1.6557;    c_store(2) = 1.7400;
% % % 6
% % a_store(3) = 0.2876;    b_store(3) = 1.6504;    c_store(3) = 1.7546;
% % % 8
% % a_store(4) = 0.2821;    b_store(4) = 1.6472;    c_store(4) = 1.7733;
% % % 10
% % a_store(5) = 0.2804;    b_store(5) = 1.6997;    c_store(5) = 1.7841;
% % % 12
% % a_store(6) = 0.2759;    b_store(6) = 1.7059;    c_store(6) = 1.8013;
% % % 14
% % a_store(7) = 0.2720;    b_store(7) = 1.7688;    c_store(7) = 1.8203;
% % % 16
% % a_store(8) = 0.2664;    b_store(8) = 1.8181;    c_store(8) = 1.8448;
% % % 18
% % a_store(9) = 0.2600;    b_store(9) = 1.7989;    c_store(9) = 1.8682;
% % % 20
% % a_store(10) = 0.2502;    b_store(10) = 1.6987;    c_store(10) = 1.8985;
% % % 22
% % a_store(11) = 0.2447;    b_store(11) = 1.7681;    c_store(11) = 1.9259;
% % % 24
% % a_store(12) = 0.2409;    b_store(12) = 1.8327;    c_store(12) = 1.9464;
% % % 26
% % a_store(13) = 0.2261;    b_store(13) = 1.5972;    c_store(13) = 1.9903;
% % % 28
% % a_store(14) = 0.2204;    b_store(14) = 1.6532;    c_store(14) = 2.0207;
% % % 30
% % a_store(15) = 0.2108;    b_store(15) = 1.5730;    c_store(15) = 2.0576;
% % % 32
% % a_store(16) = 0.2113;    b_store(16) = 1.7356;    c_store(16) = 2.0714;
% % % 34
% % a_store(17) = 0.2049;    b_store(17) = 1.7741;    c_store(17) = 2.1037;
% % % 36
% % a_store(18) = 0.1972;    b_store(18) = 1.7537;    c_store(18) = 2.1416;
% % % 38
% % a_store(19) = 0.1926;    b_store(19) = 1.8401;    c_store(19) = 2.1709;
% % % 40
% % a_store(20) = 0.1801;    b_store(20) = 1.7001;    c_store(20) = 2.2266;
% % % 42
% % a_store(21) = 0.1721;    b_store(21) = 1.6373;    c_store(21) = 2.2670;
% % % 44
% % a_store(22) = 0.1733;    b_store(22) = 1.9017;    c_store(22) = 2.2821;
% % % 46
% % a_store(23) = 0.1774;    b_store(23) = 2.3100;    c_store(23) = 2.2818;
% % % 48
% % a_store(24) = 0.1614;    b_store(24) = 1.9127;    c_store(24) = 2.3533;
% % % 50
% % a_store(25) = 0.1600;    b_store(25) = 1.9867;    c_store(25) = 2.3677;
% 
% figure
% hold on
% h = legend('show','location','best');
% percent_change = 0;
% name = num2str(percent_change);
% plot(t, Ca_i_ATP,'DisplayName',name,'linewidth',2,'color',[0.2,0.1333,0.5333])
% plot(t, Ca_i_dATP,':','DisplayName',name,'linewidth',2,'color',[0.5 0.5 0.5])
% 
% % 0.26667, 0.6667, 0.6
% 
% % for j = 1:length(a_store)
% %     x_final = [a_store(j), b_store(j), c_store(j)];
% % 
% %     for i=1:stim_period+1
% %     phi = mod(t(i)+0.001,stim_period)/stim_period;
% %     Ca_i(i) = (x_final(1)/phi)*exp(-x_final(2)*(log(phi)+x_final(3))^2) + Ca0;
% %     end
% % 
% %     [min_Ca,idx_min] = min(Ca_i);
% %     [max_Ca,idx_max] = max(Ca_i);
% %     TTP = t(idx_max);
% %     DT50 = (max(Ca_i)-min(Ca_i))*0.5 + min(Ca_i);
% %     time_DT50 = 0;
% %     for l = idx_max:stim_period+1
% %         if Ca_i(l) < DT50
% %             time_DT50 = (t(l)+t(l-1))/2;
% %             time_DT50 = time_DT50;
% %             break
% %         end
% %     end
% % 
% %     DT90 = (max(Ca_i)-min(Ca_i))*0.1 + min(Ca_i);
% %     time_DT90 = 0;
% %     for l = idx_max:stim_period+1
% %         if Ca_i(l) < DT90
% %             time_DT90 = (t(l)+t(l-1))/2;
% %             time_DT90 = time_DT90;
% %             break
% %         end
% %     end
% %     DT50_store(j) = time_DT50;
% %     DT90_store(j) = time_DT90;
% % 
% % 
% %     percent_change = j*2;
% %     name = num2str(percent_change);
% %     plot(t, Ca_i,'DisplayName',name,'linewidth',2,'color',[0.2+0.00256*j,0.1333+0.0205*j,0.5333+0.00256*j])
% % end
% 
% p = 0;
% for j = 2:2:18
%     p = p+1
%     x_final = [a_store(j), b_store(j), c_store(j)];
% 
%     for i=1:stim_period+1
%     phi = mod(t(i)+0.001,stim_period)/stim_period;
%     Ca_i(i) = (x_final(1)/phi)*exp(-x_final(2)*(log(phi)+x_final(3))^2) + Ca0;
%     end
% % 
% %     [min_Ca,idx_min] = min(Ca_i);
% %     [max_Ca,idx_max] = max(Ca_i);
% %     TTP = t(idx_max);
% %     DT50 = (max(Ca_i)-min(Ca_i))*0.5 + min(Ca_i);
% %     time_DT50 = 0;
% %     for l = idx_max:stim_period+1
% %         if Ca_i(l) < DT50
% %             time_DT50 = (t(l)+t(l-1))/2;
% %             time_DT50 = time_DT50;
% %             break
% %         end
% %     end
% % 
% %     DT90 = (max(Ca_i)-min(Ca_i))*0.1 + min(Ca_i);
% %     time_DT90 = 0;
% %     for l = idx_max:stim_period+1
% %         if Ca_i(l) < DT90
% %             time_DT90 = (t(l)+t(l-1))/2;
% %             time_DT90 = time_DT90;
% %             break
% %         end
% %     end
% %     DT50_store(j) = time_DT50;
% %     DT90_store(j) = time_DT90;
% 
% 
%     percent_change = j*2;
%     name = num2str(percent_change);
%     plot(t, Ca_i,'DisplayName',name,'linewidth',2,'color',[0.2+0.00677*p,0.1333+0.05334*p,0.5333+0.006676*p])
% end
% 
% 


%fit a, b, and c in Ca2+ equation
x0 = [a, b, c];
x_final = [a, b, c];
%x_final = fminsearch(@objective_function,x0)


for i=1:stim_period+1
    phi = mod(t(i)+0.001,stim_period)/stim_period;
    Ca_i_dATP(i) = (x_final(1)/phi)*exp(-x_final(2)*(log(phi)+x_final(3))^2) + Ca0;
end

% figure
% hold on
% plot(t,Ca_i_ATP)
% plot(t,Ca_i_dATP)
% legend('ATP','dATP')


%calculate DT50 and DT90 for fitted Ca2+ transient
[min_Ca,idx_min] = min(Ca_i)
[max_Ca,idx_max] = max(Ca_i)
TTP = t(idx_max)
DT50 = (max(Ca_i)-min(Ca_i))*0.5 + min(Ca_i);
time_DT50 = 0;
for i = idx_max:stim_period+1
    if Ca_i(i) < DT50
        time_DT50 = (t(i)+t(i-1))/2;
        time_DT50 = time_DT50
        break
    end
end

DT90 = (max(Ca_i)-min(Ca_i))*0.1 + min(Ca_i);
time_DT90 = 0;
for i = idx_max:stim_period+1
    if Ca_i(i) < DT90
        time_DT90 = (t(i)+t(i-1))/2;
        time_DT90 = time_DT90
        break
    end
end


% calculate objective function 
function SSE = objective_function(x)
SSE = 0;
stim_period = 1000; 
Ca0 = 0; 
for t=1:stim_period+1
    phi = mod(t+0.001,stim_period)/stim_period;
    Ca_i(t) = (x(1)/phi)*exp(-x(2)*(log(phi)+x(3))^2) + Ca0;
end

t = linspace(0,stim_period,stim_period+1);
[min_Ca,idx_min] = min(Ca_i);
[max_Ca,idx_max] = max(Ca_i);
TTP = t(idx_max);
time_DT50 = 0;
DT50 = (max(Ca_i)-min(Ca_i))*0.5 + min(Ca_i);
for i = idx_max:stim_period+1
    if Ca_i(i) < DT50
        time_DT50 = (t(i)+t(i-1))/2;
        time_DT50 = time_DT50;
        break
    end
end

DT90 = (max(Ca_i)-min(Ca_i))*0.1 + min(Ca_i);
time_DT90 = 0;
for i = idx_max:stim_period+1
    if Ca_i(i) < DT90
        time_DT90 = (t(i)+t(i-1))/2;
        time_DT90 = time_DT90;
        break
    end
end

%SSE = (((time_DT50 - 147.3)^2) + ((time_DT90 - 477.3)^2) + ((max_Ca - 0.34)^2));
SSE = (((time_DT50 - 119)^2) + ((time_DT90 - 278)^2) + ((max_Ca - 1)^2));


end