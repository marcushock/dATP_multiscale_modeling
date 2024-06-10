%% Digitized experimental data (interpolating at specified points)
% ***Run this section before any of other figure sections
pCa = [7 6.9 6.8 6.7 6.6 6.5 6.4 6.3 6.2 6.1 6 5.9 5.8 5.7 5.6 5.5 5.4 5.3 5.2 5.1 5 4.9 4.8 4.7 4.6 4.5 4.4 4.3 4.2 4.1 4];

% Instructions for running simulations to generate model output for figures
% can be found in github readme
%% Figure 3
ATP = [0	0	0.0000000480769	0	0	0	0	0.000000620188	0	0.0000449956	0.0000370382	0.000029602	0.000191319	0.000443113	0.00107737	0.00206548	0.00505044	0.00706438	0.00777289	0.0086795	0.00916643	0.00895904	0.00944259	0.00941371	0.00930686	0.00892618	0.00963093	0.00947331	0.00954831	0.0101342	0.00994253];
dATP_kf = [0	0	0	0	0	0	0	0	0	0.0000255961	0.0000203373	0.0000471051	0.000158947	0.000336889	0.000737122	0.00253851	0.00415322	0.0055948	0.00696991	0.00825973	0.00923599	0.0085854	0.00985886	0.00931049	0.0100148	0.0100014	0.010324	0.0102593	0.0104779	0.0090965	0.0103356];
dATP_kf_kp_kg = [0	0	0	0	0	0	0	0	0	0	0.0000366064	0.0000651911	0.00014268	0.000232423	0.000898265	0.00256978	0.00471865	0.00659074	0.00751959	0.00844057	0.00871045	0.00918688	0.0092527	0.00947122	0.00986729	0.00936085	0.00901622	0.0102164	0.0102386	0.0105029	0.0109011];
dATP_krecruit = [0	0	0	0	0	0	0	0	0.00000642817	0.0000507169	0.0000359075	0.0000209634	0.000204544	0.000599876	0.00134884	0.00369611	0.00587672	0.0084652	0.00931444	0.00980871	0.0105115	0.0114194	0.0121317	0.0111479	0.0120782	0.0124053	0.0119912	0.0120563	0.0125642	0.0118705	0.0129944];
dATP_kf_kp_kg_krecruit = [0	0	0	0	0	0	0	0	0.000000699512	0.00000233863	0.00000115442	0.0000616602	0.000138807	0.000469877	0.00144128	0.00336739	0.00607202	0.00758845	0.00920166	0.010582	0.0112417	0.01176	0.0115009	0.0115396	0.0124816	0.0119165	0.0117041	0.0119955	0.012157	0.012921	0.0122413];
dATP_kf_kp_kg_krecruit_gamma_B = [0.00383517	0.00521872	0.00540495	0.00632056	0.007565	0.00850565	0.00932302	0.0100549	0.0100657	0.0108907	0.0119236	0.0117683	0.0122565	0.0120423	0.0122047	0.0129137	0.0127484	0.0142768	0.0128392	0.0136256	0.0128389	0.0130955	0.0138736	0.0128101	0.0129578	0.0131323	0.0131138	0.0142769	0.012893	0.0129391	0.0129192];
dATP_kf_kp_kg_krecruit_mu_M = [0	0	0	0	0	0	0	0	0	0.0000680974	0.0000166834	0.0000357209	0.000172647	0.00058834	0.00102954	0.00337169	0.00518987	0.00716667	0.0081496	0.0089742	0.00923048	0.00853074	0.00994725	0.0107558	0.00932373	0.00983532	0.00978726	0.00967046	0.00976285	0.0107841	0.00944938];

[hill_ATP ec50_n_ATP] = pCa_calculate(pCa,ATP);
A_ATP = min(ATP);
D_ATP = max(ATP);
x = pCa;
cf_ATP = D_ATP+(A_ATP-D_ATP)./(1+(x/ec50_n_ATP).^hill_ATP);

[hill_dATP_kf ec50_n_dATP_kf] = pCa_calculate(pCa,dATP_kf);
A_dATP_kf = min(dATP_kf);
D_dATP_kf = max(dATP_kf);
x = pCa;
cf_dATP_kf = D_dATP_kf+(A_dATP_kf-D_dATP_kf)./(1+(x/ec50_n_dATP_kf).^hill_dATP_kf);

[hill_dATP_kf_kp_kg ec50_n_dATP_kf_kp_kg] = pCa_calculate(pCa,dATP_kf_kp_kg);
A_dATP_kf_kp_kg = min(dATP_kf_kp_kg);
D_dATP_kf_kp_kg = max(dATP_kf_kp_kg);
x = pCa;
cf_dATP_kf_kp_kg = D_dATP_kf_kp_kg+(A_dATP_kf_kp_kg-D_dATP_kf_kp_kg)./(1+(x/ec50_n_dATP_kf_kp_kg).^hill_dATP_kf_kp_kg);

[hill_dATP_kf_kp_kg_krecruit ec50_n_dATP_kf_kp_kg_krecruit] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit);
A_dATP_kf_kp_kg_krecruit = min(dATP_kf_kp_kg_krecruit);
D_dATP_kf_kp_kg_krecruit = max(dATP_kf_kp_kg_krecruit);
x = pCa;
cf_dATP_kf_kp_kg_krecruit = D_dATP_kf_kp_kg_krecruit+(A_dATP_kf_kp_kg_krecruit-D_dATP_kf_kp_kg_krecruit)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit).^hill_dATP_kf_kp_kg_krecruit);

[hill_dATP_krecruit ec50_n_dATP_krecruit] = pCa_calculate(pCa,dATP_krecruit);
A_dATP_krecruit = min(dATP_krecruit);
D_dATP_krecruit = max(dATP_krecruit);
x = pCa;
cf_dATP_krecruit = D_dATP_krecruit+(A_dATP_krecruit-D_dATP_krecruit)./(1+(x/ec50_n_dATP_krecruit).^hill_dATP_krecruit);

[hill_dATP_kf_kp_kg_krecruit_gamma_B ec50_n_dATP_kf_kp_kg_krecruit_gamma_B] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit_gamma_B);
A_dATP_kf_kp_kg_krecruit_gamma_B = min(dATP_kf_kp_kg_krecruit_gamma_B);
D_dATP_kf_kp_kg_krecruit_gamma_B = max(dATP_kf_kp_kg_krecruit_gamma_B);
x = pCa;
cf_dATP_kf_kp_kg_krecruit_gamma_B = D_dATP_kf_kp_kg_krecruit_gamma_B+(A_dATP_kf_kp_kg_krecruit_gamma_B-D_dATP_kf_kp_kg_krecruit_gamma_B)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit_gamma_B).^hill_dATP_kf_kp_kg_krecruit_gamma_B);

[hill_dATP_kf_kp_kg_krecruit_mu_M ec50_n_dATP_kf_kp_kg_krecruit_mu_M] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit_mu_M);
A_dATP_kf_kp_kg_krecruit_mu_M = min(dATP_kf_kp_kg_krecruit_mu_M);
D_dATP_kf_kp_kg_krecruit_mu_M = max(dATP_kf_kp_kg_krecruit_mu_M);
x = pCa;
cf_dATP_kf_kp_kg_krecruit_mu_M = D_dATP_kf_kp_kg_krecruit_mu_M+(A_dATP_kf_kp_kg_krecruit_mu_M-D_dATP_kf_kp_kg_krecruit_mu_M)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit_mu_M).^hill_dATP_kf_kp_kg_krecruit_mu_M);


max_force_store = [max(cf_ATP), max(cf_dATP_kf), max(cf_dATP_kf_kp_kg), max(cf_dATP_krecruit), max(cf_dATP_kf_kp_kg_krecruit)];

x_axis = [0 1 2 3];

figure
hold on
plot(pCa,cf_ATP./max(ATP),'LineWidth',3,'color',[0.2, 0.1333, 0.5333])
plot(pCa,cf_dATP_kf_kp_kg_krecruit./max(ATP),'LineWidth',3,'color',[0.2667, 0.6667, 0.6])
plot(pCa,cf_dATP_kf_kp_kg_krecruit_gamma_B./max(ATP),'LineWidth',3,'color',[0.5 0.5 0.5])
plot(pCa,cf_dATP_kf_kp_kg_krecruit_mu_M./max(ATP),'LineWidth',3,'color',[0.5 0.5 0.5])
legend('ATP','dATP','dATP - \gamma_B','dATP - \mu_M')
set(gca, 'xdir', 'reverse')
xlim([4 7])
ylim([0 1.4])
set(gca,'FontSize',14)
xlabel('pCa')
ylabel('Relative Force')

figure
hold on
bar(x_axis,(max_force_store(2:5)-max_force_store(1))./max_force_store(1),'facecolor',[0.5, 0.5, 0.5])
yline((max_force_store(1)-max_force_store(1))./max_force_store(1),':','color',[0.2, 0.1333, 0.5333],'LineWidth',3);
yline(0.4,':','color',[0.2667, 0.6667, 0.6],'LineWidth',3);
ylabel('Max Steady State Force')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

%% Figure S8
ATP = [0	0	0.0000000480769	0	0	0	0	0.000000620188	0	0.0000449956	0.0000370382	0.000029602	0.000191319	0.000443113	0.00107737	0.00206548	0.00505044	0.00706438	0.00777289	0.0086795	0.00916643	0.00895904	0.00944259	0.00941371	0.00930686	0.00892618	0.00963093	0.00947331	0.00954831	0.0101342	0.00994253];
dATP_kf = [0	0	0	0	0	0	0	0	0	0.0000255961	0.0000203373	0.0000471051	0.000158947	0.000336889	0.000737122	0.00253851	0.00415322	0.0055948	0.00696991	0.00825973	0.00923599	0.0085854	0.00985886	0.00931049	0.0100148	0.0100014	0.010324	0.0102593	0.0104779	0.0090965	0.0103356];
dATP_kf_kp_kg = [0	0	0	0	0	0	0	0	0	0	0.0000366064	0.0000651911	0.00014268	0.000232423	0.000898265	0.00256978	0.00471865	0.00659074	0.00751959	0.00844057	0.00871045	0.00918688	0.0092527	0.00947122	0.00986729	0.00936085	0.00901622	0.0102164	0.0102386	0.0105029	0.0109011];
dATP_krecruit = [0	0	0	0	0	0	0	0	0.00000642817	0.0000507169	0.0000359075	0.0000209634	0.000204544	0.000599876	0.00134884	0.00369611	0.00587672	0.0084652	0.00931444	0.00980871	0.0105115	0.0114194	0.0121317	0.0111479	0.0120782	0.0124053	0.0119912	0.0120563	0.0125642	0.0118705	0.0129944];
dATP_kf_kp_kg_krecruit = [0	0	0	0	0	0	0	0	0.000000699512	0.00000233863	0.00000115442	0.0000616602	0.000138807	0.000469877	0.00144128	0.00336739	0.00607202	0.00758845	0.00920166	0.010582	0.0112417	0.01176	0.0115009	0.0115396	0.0124816	0.0119165	0.0117041	0.0119955	0.012157	0.012921	0.0122413];
dATP_kf_kp_kg_krecruit_gamma_B = [0.00383517	0.00521872	0.00540495	0.00632056	0.007565	0.00850565	0.00932302	0.0100549	0.0100657	0.0108907	0.0119236	0.0117683	0.0122565	0.0120423	0.0122047	0.0129137	0.0127484	0.0142768	0.0128392	0.0136256	0.0128389	0.0130955	0.0138736	0.0128101	0.0129578	0.0131323	0.0131138	0.0142769	0.012893	0.0129391	0.0129192];
dATP_kf_kp_kg_krecruit_gamma_M = [0.0000121656	0.0000311877	0.0000225621	0	0.0000134868	0.00000666545	0.0000347002	0.0000958632	0.0000720472	0.000119026	0.000235999	0.000397217	0.000641188	0.00148303	0.0028047	0.00508223	0.00745079	0.0104135	0.0112876	0.011497	0.0133508	0.0128069	0.0132829	0.0143916	0.0129132	0.0133899	0.0134939	0.0135051	0.013417	0.0134485	0.0139339];
dATP_kf_kp_kg_krecruit_mu_M = [0	0	0	0	0	0	0	0	0	0.0000680974	0.0000166834	0.0000357209	0.000172647	0.00058834	0.00102954	0.00337169	0.00518987	0.00716667	0.0081496	0.0089742	0.00923048	0.00853074	0.00994725	0.0107558	0.00932373	0.00983532	0.00978726	0.00967046	0.00976285	0.0107841	0.00944938];

[hill_ATP ec50_n_ATP] = pCa_calculate(pCa,ATP);
A_ATP = min(ATP);
D_ATP = max(ATP);
x = pCa;
cf_ATP = D_ATP+(A_ATP-D_ATP)./(1+(x/ec50_n_ATP).^hill_ATP);

[hill_dATP_kf ec50_n_dATP_kf] = pCa_calculate(pCa,dATP_kf);
A_dATP_kf = min(dATP_kf);
D_dATP_kf = max(dATP_kf);
x = pCa;
cf_dATP_kf = D_dATP_kf+(A_dATP_kf-D_dATP_kf)./(1+(x/ec50_n_dATP_kf).^hill_dATP_kf);

[hill_dATP_kf_kp_kg ec50_n_dATP_kf_kp_kg] = pCa_calculate(pCa,dATP_kf_kp_kg);
A_dATP_kf_kp_kg = min(dATP_kf_kp_kg);
D_dATP_kf_kp_kg = max(dATP_kf_kp_kg);
x = pCa;
cf_dATP_kf_kp_kg = D_dATP_kf_kp_kg+(A_dATP_kf_kp_kg-D_dATP_kf_kp_kg)./(1+(x/ec50_n_dATP_kf_kp_kg).^hill_dATP_kf_kp_kg);

[hill_dATP_kf_kp_kg_krecruit ec50_n_dATP_kf_kp_kg_krecruit] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit);
A_dATP_kf_kp_kg_krecruit = min(dATP_kf_kp_kg_krecruit);
D_dATP_kf_kp_kg_krecruit = max(dATP_kf_kp_kg_krecruit);
x = pCa;
cf_dATP_kf_kp_kg_krecruit = D_dATP_kf_kp_kg_krecruit+(A_dATP_kf_kp_kg_krecruit-D_dATP_kf_kp_kg_krecruit)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit).^hill_dATP_kf_kp_kg_krecruit);

[hill_dATP_krecruit ec50_n_dATP_krecruit] = pCa_calculate(pCa,dATP_krecruit);
A_dATP_krecruit = min(dATP_krecruit);
D_dATP_krecruit = max(dATP_krecruit);
x = pCa;
cf_dATP_krecruit = D_dATP_krecruit+(A_dATP_krecruit-D_dATP_krecruit)./(1+(x/ec50_n_dATP_krecruit).^hill_dATP_krecruit);

[hill_dATP_kf_kp_kg_krecruit_gamma_B ec50_n_dATP_kf_kp_kg_krecruit_gamma_B] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit_gamma_B);
A_dATP_kf_kp_kg_krecruit_gamma_B = min(dATP_kf_kp_kg_krecruit_gamma_B);
D_dATP_kf_kp_kg_krecruit_gamma_B = max(dATP_kf_kp_kg_krecruit_gamma_B);
x = pCa;
cf_dATP_kf_kp_kg_krecruit_gamma_B = D_dATP_kf_kp_kg_krecruit_gamma_B+(A_dATP_kf_kp_kg_krecruit_gamma_B-D_dATP_kf_kp_kg_krecruit_gamma_B)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit_gamma_B).^hill_dATP_kf_kp_kg_krecruit_gamma_B);

[hill_dATP_kf_kp_kg_krecruit_gamma_M ec50_n_dATP_kf_kp_kg_krecruit_gamma_M] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit_gamma_M);
A_dATP_kf_kp_kg_krecruit_gamma_M = min(dATP_kf_kp_kg_krecruit_gamma_M);
D_dATP_kf_kp_kg_krecruit_gamma_M = max(dATP_kf_kp_kg_krecruit_gamma_M);
x = pCa;
cf_dATP_kf_kp_kg_krecruit_gamma_M = D_dATP_kf_kp_kg_krecruit_gamma_M+(A_dATP_kf_kp_kg_krecruit_gamma_M-D_dATP_kf_kp_kg_krecruit_gamma_M)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit_gamma_M).^hill_dATP_kf_kp_kg_krecruit_gamma_M);

[hill_dATP_kf_kp_kg_krecruit_mu_M ec50_n_dATP_kf_kp_kg_krecruit_mu_M] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit_mu_M);
A_dATP_kf_kp_kg_krecruit_mu_M = min(dATP_kf_kp_kg_krecruit_mu_M);
D_dATP_kf_kp_kg_krecruit_mu_M = max(dATP_kf_kp_kg_krecruit_mu_M);
x = pCa;
cf_dATP_kf_kp_kg_krecruit_mu_M = D_dATP_kf_kp_kg_krecruit_mu_M+(A_dATP_kf_kp_kg_krecruit_mu_M-D_dATP_kf_kp_kg_krecruit_mu_M)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit_mu_M).^hill_dATP_kf_kp_kg_krecruit_mu_M);

figure
hold on
plot(pCa,cf_ATP./max(ATP),'LineWidth',3,'color',[0.2, 0.1333, 0.5333])
plot(pCa,cf_dATP_kf./max(ATP),':','LineWidth',3,'color',[0.5 0.5 0.5])
plot(pCa,cf_dATP_kf_kp_kg./max(ATP),'--','LineWidth',3,'color',[0.5 0.5 0.5])
plot(pCa,cf_dATP_krecruit./max(ATP),'-.','LineWidth',3,'color',[0.5 0.5 0.5])
plot(pCa,cf_dATP_kf_kp_kg_krecruit./max(ATP),'LineWidth',3,'color',[0.2667, 0.6667, 0.6])
legend('ATP','dATP, k_f^+', 'dATP, k_f^+ + k_p^+ + k_g^+', 'dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit}', 'dATP, k_{recruit}')
set(gca, 'xdir', 'reverse')
xlim([4 7])
ylim([0 1.4])
xlabel('pCa')
ylabel('Relative Force')
set(gca,'FontSize',14)

figure
hold on
plot(pCa,cf_ATP./max(ATP),'LineWidth',3,'color',[0.2, 0.1333, 0.5333])
plot(pCa,cf_dATP_kf_kp_kg_krecruit./max(ATP),'LineWidth',3,'color',[0.2667, 0.6667, 0.6])
plot(pCa,cf_dATP_kf_kp_kg_krecruit_gamma_B./max(ATP),'LineWidth',3,'color',[0.5 0.5 0.5])
plot(pCa,cf_dATP_kf_kp_kg_krecruit_gamma_M./max(ATP),'LineWidth',3,'color',[0.6667 0.2667 0.6])
plot(pCa,cf_dATP_kf_kp_kg_krecruit_mu_M./max(ATP),'LineWidth',3,'color',[0.6667+0.0286*7 0.2667+0.0762*7 0.6-0.0190*7])
set(gca, 'xdir', 'reverse')
legend('ATP', 'dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit}','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \gamma_B','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \gamma_M','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \mu_M')
xlabel('pCa')
ylabel('Relative Force')
set(gca,'FontSize',14)

%% Figure S9
k_force = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 710, 720, 730, 740, 750, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]; %0.9; %[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
k_plus_SR_ref = 16; 
k_minus_SR_ref = 15; 

l = 0;
for i = 1:length(k_force)
    for j = 1:length(k_plus_SR_ref)
        for k = 1:length(k_minus_SR_ref)
        k_force_str = compose("%1.6f",k_force(i));
        k_plus_SR_ref_str = compose("%1.6f",k_plus_SR_ref(j));
        k_minus_SR_ref_str = compose("%1.6f",k_minus_SR_ref(k));
        filepath = "SA_krecruit/Force_pCa_Optmz dATP 0.010000 k_force " + k_force_str + " k_plus_SR_ref " + k_plus_SR_ref_str + " k_minus_SR_ref " + k_minus_SR_ref_str + ".csv";
        name = "k_{force} = " + k_force_str + ", k_{plus, SR} = " + k_plus_SR_ref_str + ", k_{minus, SR} = " + k_minus_SR_ref_str;
        data = readmatrix(filepath);
        fold_change(i,j,k) = max(data(:,2))./max(ATP);
        end
    end
end

k_recruit_vals = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 710, 720, 730, 740, 750, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000];
figure
hold on
plot(k_recruit_vals,fold_change,'o')
xlabel('k_{recruit}')
ylabel('Relative Increase in Max Steady State Force')


ATP = [0	0	0.0000000480769	0	0	0	0	0.000000620188	0	0.0000449956	0.0000370382	0.000029602	0.000191319	0.000443113	0.00107737	0.00206548	0.00505044	0.00706438	0.00777289	0.0086795	0.00916643	0.00895904	0.00944259	0.00941371	0.00930686	0.00892618	0.00963093	0.00947331	0.00954831	0.0101342	0.00994253];
dATP_kf_kp_kg = [0	0	0	0	0	0	0	0	0	0	0.0000366064	0.0000651911	0.00014268	0.000232423	0.000898265	0.00256978	0.00471865	0.00659074	0.00751959	0.00844057	0.00871045	0.00918688	0.0092527	0.00947122	0.00986729	0.00936085	0.00901622	0.0102164	0.0102386	0.0105029	0.0109011];
dATP_kf_kp_kg_inckp = [0	0	0	0	0	0	0	0	0.00000314794	0.0000273239	0.0000282756	0.0000394855	0.000162758	0.00048906	0.00121303	0.00250622	0.00435949	0.00652864	0.0074768	0.00811639	0.00820373	0.00823955	0.00823949	0.00935598	0.00957472	0.00941191	0.00954726	0.00964436	0.00963114	0.00980181	0.00958236];
dATP_kf_kp_kg_inckg = [0	0	0	0	0	0	0	0.0000122915	0	0.00000773663	0.00000138548	0.0000209993	0.000181323	0.000401291	0.000897719	0.00226944	0.0044974	0.00640256	0.00742478	0.00722778	0.00813215	0.00916994	0.00850921	0.00881475	0.00886067	0.00959634	0.00919599	0.00967489	0.00961082	0.0088352	0.00974364];
dATP_kf_kp_kg_inckb = [0	0	0	0	0.00000244591	0.00000983205	0	0.0000550219	0.000180124	0.000644301	0.00201057	0.00450795	0.00647031	0.00851075	0.00907478	0.00988144	0.0102201	0.0101692	0.00968527	0.0105247	0.00981087	0.00999395	0.0112272	0.00985506	0.0102966	0.0106925	0.0105661	0.0107029	0.0104754	0.0116982	0.0110115];
dATP_kf_kp_kg_krecruit = [0	0	0	0	0	0	0	0	0.000000699512	0.00000233863	0.00000115442	0.0000616602	0.000138807	0.000469877	0.00144128	0.00336739	0.00607202	0.00758845	0.00920166	0.010582	0.0112417	0.01176	0.0115009	0.0115396	0.0124816	0.0119165	0.0117041	0.0119955	0.012157	0.012921	0.0122413];

[hill_ATP ec50_n_ATP] = pCa_calculate(pCa,ATP);
A_ATP = min(ATP);
D_ATP = max(ATP);
x = pCa;
cf_ATP = D_ATP+(A_ATP-D_ATP)./(1+(x/ec50_n_ATP).^hill_ATP);

[hill_dATP_kf ec50_n_dATP_kf] = pCa_calculate(pCa,dATP_kf);
A_dATP_kf = min(dATP_kf);
D_dATP_kf = max(dATP_kf);
x = pCa;
cf_dATP_kf = D_dATP_kf+(A_dATP_kf-D_dATP_kf)./(1+(x/ec50_n_dATP_kf).^hill_dATP_kf);

[hill_dATP_kf_kp_kg ec50_n_dATP_kf_kp_kg] = pCa_calculate(pCa,dATP_kf_kp_kg);
A_dATP_kf_kp_kg = min(dATP_kf_kp_kg);
D_dATP_kf_kp_kg = max(dATP_kf_kp_kg);
x = pCa;
cf_dATP_kf_kp_kg = D_dATP_kf_kp_kg+(A_dATP_kf_kp_kg-D_dATP_kf_kp_kg)./(1+(x/ec50_n_dATP_kf_kp_kg).^hill_dATP_kf_kp_kg);

[hill_dATP_kf_kp_kg_inckp ec50_n_dATP_kf_kp_kg_inckp] = pCa_calculate(pCa,dATP_kf_kp_kg_inckp);
A_dATP_kf_kp_kg_inckp = min(dATP_kf_kp_kg_inckp);
D_dATP_kf_kp_kg_inckp = max(dATP_kf_kp_kg_inckp);
x = pCa;
cf_dATP_kf_kp_kg_inckp = D_dATP_kf_kp_kg_inckp+(A_dATP_kf_kp_kg_inckp-D_dATP_kf_kp_kg_inckp)./(1+(x/ec50_n_dATP_kf_kp_kg_inckp).^hill_dATP_kf_kp_kg_inckp);

[hill_dATP_kf_kp_kg_inckg ec50_n_dATP_kf_kp_kg_inckg] = pCa_calculate(pCa,dATP_kf_kp_kg_inckg);
A_dATP_kf_kp_kg_inckg = min(dATP_kf_kp_kg_inckg);
D_dATP_kf_kp_kg_inckg = max(dATP_kf_kp_kg_inckg);
x = pCa;
cf_dATP_kf_kp_kg_inckg = D_dATP_kf_kp_kg_inckg+(A_dATP_kf_kp_kg_inckg-D_dATP_kf_kp_kg_inckg)./(1+(x/ec50_n_dATP_kf_kp_kg_inckg).^hill_dATP_kf_kp_kg_inckg);

[hill_dATP_kf_kp_kg_inckb ec50_n_dATP_kf_kp_kg_inckb] = pCa_calculate(pCa,dATP_kf_kp_kg_inckb);
A_dATP_kf_kp_kg_inckb = min(dATP_kf_kp_kg_inckb);
D_dATP_kf_kp_kg_inckb = max(dATP_kf_kp_kg_inckb);
x = pCa;
cf_dATP_kf_kp_kg_inckb = D_dATP_kf_kp_kg_inckb+(A_dATP_kf_kp_kg_inckb-D_dATP_kf_kp_kg_inckb)./(1+(x/ec50_n_dATP_kf_kp_kg_inckb).^hill_dATP_kf_kp_kg_inckb);

[hill_dATP_kf_kp_kg_krecruit ec50_n_dATP_kf_kp_kg_krecruit] = pCa_calculate(pCa,dATP_kf_kp_kg_krecruit);
A_dATP_kf_kp_kg_krecruit = min(dATP_kf_kp_kg_krecruit);
D_dATP_kf_kp_kg_krecruit = max(dATP_kf_kp_kg_krecruit);
x = pCa;
cf_dATP_kf_kp_kg_krecruit = D_dATP_kf_kp_kg_krecruit+(A_dATP_kf_kp_kg_krecruit-D_dATP_kf_kp_kg_krecruit)./(1+(x/ec50_n_dATP_kf_kp_kg_krecruit).^hill_dATP_kf_kp_kg_krecruit);

figure
hold on
plot(pCa,cf_ATP./max(ATP),'LineWidth',3,'color',[0.2, 0.1333, 0.5333])
plot(pCa,cf_dATP_kf_kp_kg./max(ATP),'LineWidth',3)
plot(pCa,cf_dATP_kf_kp_kg_inckp./max(ATP),'LineWidth',3)
plot(pCa,cf_dATP_kf_kp_kg_inckg./max(ATP),'LineWidth',3)
plot(pCa,cf_dATP_kf_kp_kg_inckb./max(ATP),'LineWidth',3)
plot(pCa,cf_dATP_kf_kp_kg_krecruit./max(ATP),'LineWidth',3,'color',[0.2667, 0.6667, 0.6])
legend('ATP','dATP, k_f^+ + k_p^+ + k_g^+','dATP, k_f^+ + k_p^+ + k_g^+, increasing k_p^+','dATP, k_f^+ + k_p^+ + k_g^+, increasing k_g^+','dATP, k_f^+ + k_p^+ + k_g^+, increasing k_b^+','dATP, k_f^+ + k_p^+ + k_g^+, increasing k_recruit^+')
set(gca, 'xdir', 'reverse')
xlim([4 7])
ylim([0 1.4])


%% Figure S4
States_ATP = readmatrix("States/States_ATP.csv");
States_dATP_kf_kp_kg = readmatrix("States/States_dATP_kf_kp_kg.csv");
States_dATP_kf_kp_kg_krecruit = readmatrix("States/States_dATP_kf_kp_kg_krecruit.csv");
States_dATP_kf_kp_kg_krecruit_gammaB = readmatrix("States/States_dATP_kf_kp_kg_krecruit_gammaB.csv");
States_dATP_kf_kp_kg_krecruit_gammaM = readmatrix("States/States_dATP_kf_kp_kg_krecruit_gammaM.csv");
States_dATP_kf_kp_kg_krecruit_muM = readmatrix("States/States_dATP_kf_kp_kg_krecruit_muM.csv");

Y = [min(States_ATP(100:end,152)) min(States_dATP_kf_kp_kg(100:end,152)) min(States_dATP_kf_kp_kg_krecruit(100:end,152)) min(States_dATP_kf_kp_kg_krecruit_gammaB(100:end,152)) min(States_dATP_kf_kp_kg_krecruit_gammaM(100:end,152)) min(States_dATP_kf_kp_kg_krecruit_muM(100:end,152)); ...
min(States_ATP(100:end,153)) min(States_dATP_kf_kp_kg(100:end,153)) min(States_dATP_kf_kp_kg_krecruit(100:end,153)) min(States_dATP_kf_kp_kg_krecruit_gammaB(100:end,153)) min(States_dATP_kf_kp_kg_krecruit_gammaM(100:end,153)) min(States_dATP_kf_kp_kg_krecruit_muM(100:end,153)); ...
min(States_ATP(100:end,154)) min(States_dATP_kf_kp_kg(100:end,154)) min(States_dATP_kf_kp_kg_krecruit(100:end,154)) min(States_dATP_kf_kp_kg_krecruit_gammaB(100:end,154)) min(States_dATP_kf_kp_kg_krecruit_gammaM(100:end,154)) min(States_dATP_kf_kp_kg_krecruit_muM(100:end,154)); ...
min(States_ATP(100:end,155)) min(States_dATP_kf_kp_kg(100:end,155)) min(States_dATP_kf_kp_kg_krecruit(100:end,155)) min(States_dATP_kf_kp_kg_krecruit_gammaB(100:end,155)) min(States_dATP_kf_kp_kg_krecruit_gammaM(100:end,155)) min(States_dATP_kf_kp_kg_krecruit_muM(100:end,155)); ...
min(States_ATP(100:end,156)) min(States_dATP_kf_kp_kg(100:end,156)) min(States_dATP_kf_kp_kg_krecruit(100:end,156)) min(States_dATP_kf_kp_kg_krecruit_gammaB(100:end,156)) min(States_dATP_kf_kp_kg_krecruit_gammaM(100:end,156)) min(States_dATP_kf_kp_kg_krecruit_muM(100:end,156)); ...
];

figure
hold on
bar(Y)
legend('ATP','dATP, k_f^+ + k_p^+ + k_g^+','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit}','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \gammaB','dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \gammaM', 'dATP, k_f^+ + k_p^+ + k_g^+ + k_{recruit} - \muM')
xlabel('State')
ylabel('State Occupancy')
set(gca,'FontSize',14)

%% Figure S12
States_ATP = readmatrix("States/States_ATP.csv");
figure
hold on
plot(States_ATP(:,1),States_ATP(:,152),'LineWidth',3)
plot(States_ATP(:,1),States_ATP(:,153),'LineWidth',3)
plot(States_ATP(:,1),States_ATP(:,154),'LineWidth',3)
plot(States_ATP(:,1),States_ATP(:,155),'LineWidth',3)
plot(States_ATP(:,1),States_ATP(:,156),'LineWidth',3)
legend('M1','M2','C','B','USR')

