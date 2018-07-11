% load the control initial conditions files
%cd('/Users/Air/Documents/PROFESSIONAL/Projects/CMRG/CaMKII/CaMKII model/Paper versions/Outputs')
load('TG_x0');
load('WT_x0');
load('NoEAD_APClamp');
% load('EAD_APClamp');
% Variables not used in the function calls but used to convert output currents to fluxes for plotting at end of script 

p.VJSR=0.12e-6;
p.VNSR=2.098e-6;
p.Vmyo=25.84e-6;
p.VSS = 1.458e-9;
p.Faraday=96.485;
p.Acap=1.534e-4;

% Beats for plotting

BEATS=[1,2,5,10,20,50];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%% WT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% param.beta = 0.30; % proportion of NCX current in the dyad (ratio, dimensionless)
% param.vNCXmax = 2.0; % Vmax for NCX (A/F)
% param.SERCAVmax = 0.2; % Vmax for SERCA (umol/ms)
% param.kmr = 3500; % luminal kd for SERCA (umol/L)
% param.kmf = 0.6; % cytosolic kd for SERCA (umol/L)
% param.GCa = 0.1729; % ICaL maximal conductance (A/F)
% param.SRsens = 600; % RyR luminal calcium sensitivity (kd)
% param.KmNai = 16650; % Na+-K+ pump kd for cytosolic sodium
% param.Kpcmax = 0.11662;  % Change in Bondarenko 2010
% param.koCa = 0.005;%
% param.kiCa = 0.22;
% param.kim = 1.2;
% param.kom = 0.06;
% param.MaxSR = 18;
% param.MinSR = 1; 
% %******************* Bondarenko et al 2004
% param.Gtof = 0.4067;
% param.GKss = 0.05;   
% param.GKur = 0.16;
% %******************* Bondarenko et al 2004
% param.GK1  = 0.2938;
% param.dPRyR = 0.02;
% % 
% % function call for WT model
% [T1,S1,loginfo1,x1,ioninfo1] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,50,param,WT_x0,'ioninfo');
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % % 
% params.beta = 0.30;
% %***************
% params.vNCXmax=2.0*2.5;
% %***************
% params.SERCAVmax = 0.2*0.7;
% params.kmr = 3500;
% params.kmf = 0.6;
% params.GCa = 0.1729*1.2; % Tg GCa increase
% params.SRsens = 300;
% params.KmNai = 16650; % Despa et al., Circ Res 2005
% params.Kpcmax = 0.11662*0.75;% Change in Bondarenko 2010 + Tg
% params.koCa = 0.005;
% params.kiCa = 0.22;
% params.kim = 1.2;
% params.kom = 0.06;
% params.MaxSR = 10;
% params.MinSR = 1;
% params.kleak = (0.1e-5)*5; % CaMKII leak is much higher (Maier 2003 for this value - spark frequency)
% %******************* Bondarenko et al 2004
% params.Gtof = 0.4067*0.65;   
% params.GKss = 0.05;   
% params.GKur = 0.16;
% %******************* Bondarenko et al 2004
% params.GK1  = 0.2938*0.5; % Just the Tg effect here
% params.dPRyR = 0.02;
% % Na channel parameters for Grandi 2007
% params.P2a5 = 7.8;
% params.P1b5 = 0.0084/1.8;
% params.P1b6 = 6.5449e-07;
% params.P1a7 = 3.377e-04;
% params.P1b7 = 1.868e-04;
% params.P1a8 = 6.5e-06/1000;
% params.P1b8 = 3.8e-03/1000;
% % 
% % % function call for Tg model
% [T2,S2,loginfo2,x2,ioninfo2] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,50,params,TG_x0,'ioninfo');
%  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG + Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cd('/Users/Andy/Dropbox/cmrg/CaMKII/CaMKII models/Final')
 
params2.beta = 0.30;
%******************
params2.vNCXmax=2.0*2.5;
%******************
params2.SERCAVmax = 0.2*0.7;
params2.kmr = 3500*0.55;
params2.kmf = 0.6*0.55;
%******************
params2.GCa = 0.1729*1.2*2.0; % Iso GCa increase
%******************
params2.SRsens = 300;
params2.KmNai = 16650*0.72; % Despa et al., Circ Res 2005
params2.Kpcmax = 0.11662*0.75;% Change in Bondarenko 2010 + match TGthalf (maier & kohlhaas)
params2.koCa = 0.005;
params2.kiCa = 0.22;
params2.kim = 1.2;
params2.kom = 0.06;
params2.MaxSR = 10;
params2.MinSR = 1;
params2.kleak = (0.1e-5)*5; % CaMKII leak is much higher (Maier 2003 for this value - spark frequency)
%******************* Bondarenko et al 2004
params2.Gtof = 0.4067*0.65;   
params2.GKss = 0.05;   
params2.GKur = 0.160*1.20;
%******************* Bondarenko et al 2004
params2.GK1  = 0.2938*0.5; % Just the Tg effect here
params2.dPRyR = 0.02;
params2.P2a5 = 7.8;
params2.P1b5 = 0.0084/1.8;
params2.P1b6 = 6.5449e-07;
params2.P1a7 = 3.377e-04;
params2.P1b7 = 1.868e-04;
params2.P1a8 = 6.5e-06/1000;
params2.P1b8 = 3.8e-03/1000;
params2.alphaV = 22;
params2.betaV = 22; 

params2.VClampAmp = NoEAD_APClamp(:,2)';
params2.VClampTimes = NoEAD_APClamp(:,1)';
% params2.VClampAmp = EAD_VClamp(:,2)';
% params2.VClampTimes = EAD_VClamp(:,1)';

% 
[T,S,loginfo,x,ioninfo] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,1,params2,TG_x0,'ioninfo');
 
% for i = 1:length(BEATS)
%     fluxJup{i}=cumtrapz(T1{BEATS(i)},loginfo1.ionfluxes{1,BEATS(i)}(:,5));
%     fluxSS_JNCX{i}=cumtrapz(T1{BEATS(i)},(loginfo1.ionfluxes{1,BEATS(i)}(:,6)*p.VSS/p.Vmyo));
%     fluxi_JNCX{i}=cumtrapz(T1{BEATS(i)},loginfo1.ionfluxes{1,BEATS(i)}(:,7));
%     fluxJpca{i}=cumtrapz(T1{BEATS(i)},loginfo1.ionfluxes{1,BEATS(i)}(:,4));
%     fluxJtotal{i}=cumtrapz(T1{BEATS(i)},loginfo1.ionfluxes{1,BEATS(i)}(:,3));
%     fluxJICa{i}=cumtrapz(T1{BEATS(i)},((-loginfo1.currents{1,BEATS(i)}(:,11))*p.Acap)*(p.VSS/p.Vmyo))/(2*p.Faraday*p.VSS);
%     TGfluxJup{i}=cumtrapz(T2{BEATS(i)},loginfo2.ionfluxes{1,BEATS(i)}(:,5));
%     TGfluxSS_JNCX{i}=cumtrapz(T2{BEATS(i)},(loginfo2.ionfluxes{1,BEATS(i)}(:,6)*p.VSS/p.Vmyo));
%     TGfluxi_JNCX{i}=cumtrapz(T2{BEATS(i)},loginfo2.ionfluxes{1,BEATS(i)}(:,7));
%     TGfluxJpca{i}=cumtrapz(T2{BEATS(i)},loginfo2.ionfluxes{1,BEATS(i)}(:,4));
%     TGfluxJtotal{i}=cumtrapz(T2{BEATS(i)},loginfo2.ionfluxes{1,BEATS(i)}(:,3));
%     TGfluxJICa{i}=cumtrapz(T2{BEATS(i)},((-loginfo2.currents{1,BEATS(i)}(:,11))*p.Acap)*(p.VSS/p.Vmyo))/(2*p.Faraday*p.VSS);
%     Iso_TGfluxJup{i}=cumtrapz(T3{BEATS(i)},loginfo3.ionfluxes{1,BEATS(i)}(:,5));
%     Iso_TGfluxSS_JNCX{i}=cumtrapz(T3{BEATS(i)},(loginfo3.ionfluxes{1,BEATS(i)}(:,6)*p.VSS/p.Vmyo));
%     Iso_TGfluxi_JNCX{i}=cumtrapz(T3{BEATS(i)},loginfo3.ionfluxes{1,BEATS(i)}(:,7));
%     Iso_TGfluxJpca{i}=cumtrapz(T3{BEATS(i)},loginfo3.ionfluxes{1,BEATS(i)}(:,4));
%     Iso_TGfluxJtotal{i}=cumtrapz(T3{BEATS(i)},loginfo3.ionfluxes{1,BEATS(i)}(:,3));
%     Iso_TGfluxJICa{i}=cumtrapz(T3{BEATS(i)},((-loginfo3.currents{1,BEATS(i)}(:,11))*p.Acap)*(p.VSS/p.Vmyo))/(2*p.Faraday*p.VSS);
% end
% 
% close all
% 
% varcols = [37,33,34,35,36];
% colors = {'b','r','g','k','m','c'};
% TGcolors = {'*b','*r','*g','*k','*m','*c'};
% IsoTGcolors = {'ob','or','og','ok','om','oc'};
% IsoWTcolors = {'+b','+r','+g','+k','+m','+c'};
% 
% for i = 1:length(BEATS)
%      CaSR{i}=((S1{BEATS(i)}(:,35)+loginfo1.buffers{BEATS(i)}(:,3))*p.VJSR/(p.VJSR+p.VNSR)+S1{BEATS(i)}(:,36)*p.VNSR/(p.VJSR+p.VNSR))*(p.VJSR+p.VNSR)/p.Vmyo;
%      TGCaSR{i} = ((S2{BEATS(i)}(:,35)+loginfo2.buffers{BEATS(i)}(:,3))*p.VJSR/(p.VJSR+p.VNSR)+S2{BEATS(i)}(:,36)*p.VNSR/(p.VJSR+p.VNSR))*(p.VJSR+p.VNSR)/p.Vmyo;
%      Iso_TGCaSR{i} = ((S3{BEATS(i)}(:,35)+loginfo3.buffers{BEATS(i)}(:,3))*p.VJSR/(p.VJSR+p.VNSR)+S3{BEATS(i)}(:,36)*p.VNSR/(p.VJSR+p.VNSR))*(p.VJSR+p.VNSR)/p.Vmyo;
% end
% for i = 1:length(BEATS)
%     for j = 1:length(S1{BEATS(i)}(:,20))
%         WT_PR{i}(j) = 1 - sum(S1{BEATS(i)}(j,20:22));
%     end
%     if length(WT_PR{i})>length(S1{BEATS(i)}(:,20))
%         WT_PR{i} = WT_PR{i}(1:length(S1{BEATS(i)}(:,20)));
%     end
%     for j = 1:length(S2{BEATS(i)}(:,20))
%         TG_PR{i}(j) = 1 - sum(S2{BEATS(i)}(j,20:22));
%     end
%     if length(TG_PR{i})>length(S2{BEATS(i)}(:,20))
%         TG_PR{i} = TG_PR{i}(1:length(S2{BEATS(i)}(:,20)));
%     end
%     for j = 1:length(S3{BEATS(i)}(:,20))
%         Iso_TG_PR{i}(j) = 1 - sum(S3{BEATS(i)}(j,20:22));
%     end
%     if length(Iso_TG_PR{i})>length(S3{BEATS(i)}(:,20))
%         Iso_TG_PR{i} = Iso_TG_PR{i}(1:length(S3{BEATS(i)}(:,20)));
%     end
% end
% 
% figure
% hold on
% plot(T1{BEATS(end)},S1{BEATS(end)}(:,20), 'b', T1{BEATS(end)},S1{BEATS(end)}(:,21), 'r', T1{BEATS(end)},S1{BEATS(end)}(:,22), 'g',T1{BEATS(end)},WT_PR{end}, 'k')
% plot(T2{BEATS(end)},S2{BEATS(end)}(:,20), '*b', T2{BEATS(end)},S2{BEATS(end)}(:,21), '*r', T2{BEATS(end)},S2{BEATS(end)}(:,22), '*g',T2{BEATS(end)},TG_PR{end}, '*k')
% plot(T3{BEATS(end)},S3{BEATS(end)}(:,20), 'ob', T3{BEATS(end)},S3{BEATS(end)}(:,21), 'or', T3{BEATS(end)},S3{BEATS(end)}(:,22), 'og',T3{BEATS(end)},Iso_TG_PR{end}, 'ok')
% pause
%     
% for i = 1: length(varcols)
%     figure
%     hold on
%     for j = 1: length(BEATS)
%         plot (T1{BEATS(j)},S1{BEATS(j)}(:,varcols(i)), char(colors(j)))
%         plot (T2{BEATS(j)},S2{BEATS(j)}(:,varcols(i)), char(TGcolors(j)))
%         plot (T3{BEATS(j)},S3{BEATS(j)}(:,varcols(i)), char(IsoTGcolors(j)))
%     end
%     pause
% end
% figure
% hold on
% for j = 1: length(BEATS)
%     plot (T1{BEATS(j)},CaSR{j}, char(colors(j)))
%     plot (T2{BEATS(j)},TGCaSR{j}, char(TGcolors(j)))
%     plot (T3{BEATS(j)},Iso_TGCaSR{j}, char(IsoTGcolors(j)))
% end
% pause
% figure
% hold on
% for j = 1: length(BEATS)
%     plot (T1{BEATS(j)},loginfo1.currents{BEATS(j)}(:,11), char(colors(j)))
%     plot (T2{BEATS(j)},loginfo2.currents{BEATS(j)}(:,11), char(TGcolors(j)))
%     plot (T3{BEATS(j)},loginfo3.currents{BEATS(j)}(:,11), char(IsoTGcolors(j)))
%     pause
% end
% figure
% hold on
% for j = 1: length(BEATS)
%     plot (T1{BEATS(j)},loginfo1.ionfluxes{BEATS(j)}(:,8), char(colors(j)))
%     plot (T2{BEATS(j)},loginfo2.ionfluxes{BEATS(j)}(:,8), char(TGcolors(j)))
%     plot (T3{BEATS(j)},loginfo3.ionfluxes{BEATS(j)}(:,8), char(IsoTGcolors(j)))
%     pause
% end
% figure    
% hold on
% for i = 2  
%     plot(T1{BEATS(i)},100*fluxJtotal{i}/max(fluxJtotal{i}), 'b', T1{BEATS(i)},100*fluxJup{i}/max(fluxJtotal{i}), 'r', T1{BEATS(i)},100*(fluxi_JNCX{i}+fluxSS_JNCX{i})/max(fluxJtotal{i}), 'g',T1{BEATS(i)},100*fluxSS_JNCX{i}/max(fluxJtotal{i}), 'k', T1{BEATS(i)},100*fluxJpca{i}/max(fluxJtotal{i}), 'y',T1{BEATS(i)},100*fluxJICa{i}/max(fluxJtotal{i}), 'm')
%     plot(T2{BEATS(i)},100*TGfluxJtotal{i}/max(TGfluxJtotal{i}), '*b', T2{BEATS(i)},100*TGfluxJup{i}/max(TGfluxJtotal{i}), '*r', T2{BEATS(i)},100*(TGfluxi_JNCX{i}+TGfluxSS_JNCX{i})/max(TGfluxJtotal{i}), '*g',T2{BEATS(i)},100*TGfluxSS_JNCX{i}/max(TGfluxJtotal{i}), '*k', T2{BEATS(i)},100*TGfluxJpca{i}/max(TGfluxJtotal{i}), '*y',T2{BEATS(i)},100*TGfluxJICa{i}/max(TGfluxJtotal{i}), '*m')
%     plot(T3{BEATS(i)},100*Iso_TGfluxJtotal{i}/max(Iso_TGfluxJtotal{i}), 'ob', T3{BEATS(i)},100*Iso_TGfluxJup{i}/max(Iso_TGfluxJtotal{i}), 'or', T3{BEATS(i)},100*(Iso_TGfluxi_JNCX{i}+Iso_TGfluxSS_JNCX{i})/max(Iso_TGfluxJtotal{i}), 'og',T3{BEATS(i)},100*Iso_TGfluxSS_JNCX{i}/max(Iso_TGfluxJtotal{i}), 'ok', T3{BEATS(i)},100*Iso_TGfluxJpca{i}/max(Iso_TGfluxJtotal{i}), 'oy',T3{BEATS(i)},100*Iso_TGfluxJICa{i}/max(Iso_TGfluxJtotal{i}), 'om')
%     pause
% end
% 
% for i = 1:length(BEATS)
%     [peak_Ca(i),peak_Ca_ind(i)] = max(S1{BEATS(i)}(:,33));
%     diastolic_Ca(i) = min(S1{BEATS(i)}(:,33));
%     delta_Ca(i) = peak_Ca(i)-diastolic_Ca(i);
%     SRmax(i) = max(CaSR{i});
%     FR(i) = (max(CaSR{i})-min(CaSR{i}))/max(CaSR{i});
%     x = T1{BEATS(i)}(peak_Ca_ind(i):end);
%     y = S1{BEATS(i)}(peak_Ca_ind(i):end,33);
%     caseval = 1;
%     fits = exp2fit(x,y,caseval);
%     tau(i) = fits(3)
%     [RMP(i),RMP_ind(i)] = min(S1{BEATS(i)}(:,37));
%     [peak_Vm(i),peak_Vm_ind(i)] = max(S1{BEATS(i)}(:,37));
%     delta_Vm(i) = peak_Vm(i) - RMP(i);
%     AP_50_ind(i) = find(S1{BEATS(i)}(peak_Vm_ind(i):end,37) < (peak_Vm(i)-0.5*delta_Vm(i)),1)+peak_Vm_ind(i);
%     AP_50(i) = T1{BEATS(i)}(AP_50_ind(i))-T1{BEATS(i)}(peak_Vm_ind(i));
%     AP_90_ind(i) = find(S1{BEATS(i)}(peak_Vm_ind(i):end,37) < (peak_Vm(i)-0.9*delta_Vm(i)),1)+peak_Vm_ind(i);
%     AP_90(i) = T1{BEATS(i)}(AP_90_ind(i))-T1{BEATS(i)}(peak_Vm_ind(i));
%     
%     [TG_peak_Ca(i),TG_peak_Ca_ind(i)] = max(S2{BEATS(i)}(:,33));
%     [TG_diastolic_Ca(i),TG_diastolic_Ca_ind(i)] = min(S2{BEATS(i)}(:,33));
%     TG_delta_Ca(i) = TG_peak_Ca(i)-TG_diastolic_Ca(i);
%     TGSRmax(i) = max(TGCaSR{i});
%     TG_FR(i) = (max(TGCaSR{i})-min(TGCaSR{i}))/max(TGCaSR{i});
%     TG_x = T2{BEATS(i)}(TG_peak_Ca_ind(i):end);
%     TG_y = S2{BEATS(i)}(TG_peak_Ca_ind(i):end,33);
%     caseval = 1;
%     TG_fits = exp2fit(TG_x,TG_y,caseval);
%     TG_tau(i) = TG_fits(3)
%     [TG_RMP(i),TG_RMP_ind(i)] = min(S2{BEATS(i)}(:,37));
%     [TG_peak_Vm(i),TG_peak_Vm_ind(i)] = max(S2{BEATS(i)}(:,37));
%     TG_delta_Vm(i) = TG_peak_Vm(i) - TG_RMP(i);
%     TG_AP_50_ind(i) = find(S2{BEATS(i)}(TG_peak_Vm_ind(i):end,37) < (TG_peak_Vm(i)-0.5*TG_delta_Vm(i)),1)+TG_peak_Vm_ind(i);
%     TG_AP_50(i) = T2{BEATS(i)}(TG_AP_50_ind(i))-T2{BEATS(i)}(TG_peak_Vm_ind(i));
%     TG_AP_90_ind(i) = find(S2{BEATS(i)}(TG_peak_Vm_ind(i):end,37) < (TG_peak_Vm(i)-0.9*TG_delta_Vm(i)),1)+TG_peak_Vm_ind(i);
%     TG_AP_90(i) = T2{BEATS(i)}(TG_AP_90_ind(i))-T2{BEATS(i)}(TG_peak_Vm_ind(i));
% %     
%     [Iso_TG_peak_Ca(i),Iso_TG_peak_Ca_ind(i)] = max(S3{BEATS(i)}(:,33));
%     [Iso_TG_diastolic_Ca(i),Iso_TG_diastolic_Ca_ind(i)] = min(S3{BEATS(i)}(:,33));
%     Iso_TG_delta_Ca(i) = Iso_TG_peak_Ca(i)-Iso_TG_diastolic_Ca(i);
%     Iso_TGSRmax(i) = max(Iso_TGCaSR{i});
%     Iso_TG_FR(i) = (max(Iso_TGCaSR{i})-min(Iso_TGCaSR{i}))/max(Iso_TGCaSR{i});
%     Iso_TG_x = T3{BEATS(i)}(Iso_TG_peak_Ca_ind(i):end);
%     Iso_TG_y = S3{BEATS(i)}(Iso_TG_peak_Ca_ind(i):end,33);
%     Iso_caseval = 1;
%     Iso_TG_fits = exp2fit(Iso_TG_x,Iso_TG_y,Iso_caseval);
%     Iso_TG_tau(i) = Iso_TG_fits(3)
%     [Iso_TG_RMP(i),Iso_TG_RMP_ind(i)] = min(S3{BEATS(i)}(:,37));
%     [Iso_TG_peak_Vm(i),Iso_TG_peak_Vm_ind(i)] = max(S3{BEATS(i)}(:,37));
%     Iso_TG_delta_Vm(i) = Iso_TG_peak_Vm(i) - Iso_TG_RMP(i);
%     Iso_TG_AP_50_ind(i) = find(S3{BEATS(i)}(Iso_TG_peak_Vm_ind(i):end,37) < (Iso_TG_peak_Vm(i)-0.5*Iso_TG_delta_Vm(i)),1)+Iso_TG_peak_Vm_ind(i);
%     Iso_TG_AP_50(i) = T3{BEATS(i)}(Iso_TG_AP_50_ind(i))-T3{BEATS(i)}(Iso_TG_peak_Vm_ind(i));
%     Iso_TG_AP_90_ind(i) = find(S3{BEATS(i)}(Iso_TG_peak_Vm_ind(i):end,37) < (Iso_TG_peak_Vm(i)-0.9*Iso_TG_delta_Vm(i)),1)+Iso_TG_peak_Vm_ind(i);
%     Iso_TG_AP_90(i) = T3{BEATS(i)}(Iso_TG_AP_90_ind(i))-T3{BEATS(i)}(Iso_TG_peak_Vm_ind(i));
% end
% 
% Titles = {'RMP','delta_Em','AP_50','AP_90','diastolic_Ca','delta_Ca','Ca_tau','FR','SR load'};
% 
% Outputs =[RMP(6),TG_RMP(6),Iso_TG_RMP(4);...
% delta_Vm(6),TG_delta_Vm(6),Iso_TG_delta_Vm(4);...
% AP_50(6),TG_AP_50(6),Iso_TG_AP_50(4);...
% AP_90(6),TG_AP_90(6),Iso_TG_AP_90(4);...
% diastolic_Ca(6),TG_diastolic_Ca(6),Iso_TG_diastolic_Ca(4);...
% delta_Ca(6),TG_delta_Ca(6),Iso_TG_delta_Ca(4);...
% tau(6),TG_tau(6),Iso_TG_tau(4);... 
% FR(6),TG_FR(6),Iso_TG_FR(4);...
% SRmax(6),TGSRmax(6),Iso_TGSRmax(4)];...
% 
% 
% cd('/Users/Air/Desktop/Model outs')
% 
% time = clock;
% year = num2str(time(1));
% month = num2str(time(2));
% day = num2str(time(3));
% hour = num2str(time(4));
% min = num2str(time(5));
% fid = fopen(strcat('BondMasterOut',month,day,year,'_',hour,':',min),'a');
% fprintf(fid,'%s\t','');
% fprintf(fid,'%s\t','WT');
% fprintf(fid,'%s\t','TG');
% fprintf(fid,'%s\n','TG+Iso');
% for i = 1:length(Titles)
%     fprintf(fid,'%s\t',Titles{i});
%     fprintf(fid,'%12.4f\t',Outputs(i,1));
%     fprintf(fid,'%12.4f\t',Outputs(i,2));
%     fprintf(fid,'%12.4f\n',Outputs(i,3)); 
% end
%     
% fprintf(fid,'%12.4f\n','');
% fprintf(fid,'%12.4f\n','WT');
% 
% paramfields = fields(param);
% 
% for i = 1:length(paramfields)
%     fprintf(fid,'%s\t',paramfields{i});
%     fprintf(fid,'%12.4f\n',param.(paramfields{i}));
% end
% % 
% fprintf(fid,'%12.4f\n','');
% fprintf(fid,'%12.4f\n','TG');
% 
% paramfields = fields(params);
% 
% for i = 1:length(paramfields)
%     fprintf(fid,'%s\t',paramfields{i});
%     fprintf(fid,'%12.4f\n',params.(paramfields{i}));
% end
% 
% fprintf(fid,'%12.4f\n','');
% fprintf(fid,'%12.4f\n','TG+Iso');
% 
% paramfields = fields(params2);
% 
% for i = 1:length(paramfields)
%     fprintf(fid,'%s\t',paramfields{i});
%     fprintf(fid,'%12.4f\n',params2.(paramfields{i}));
% end


 

