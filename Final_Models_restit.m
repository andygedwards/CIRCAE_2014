% load the control initial conditions files

load('TG_x0');

% Beats for plotting

BEATS=[1,2,5,10,20,50];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG + Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

params2.beta = 0.30;
%******************
params2.vNCXmax=2.0*2.6;
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
% 
[T,S,loginfo,x,ioninfo] = pacemodel(@Final_BondGrandiINa_SBSERCA,100,50,params2,TG_x0,'ioninfo');

[T1,S1,loginfo1,x1,ioninfo1] = pacemodel(@Final_BondGrandiINa_SBSERCA,200,50,params2,TG_x0,'ioninfo');

[T2,S2,loginfo2,x2,ioninfo2] = pacemodel(@Final_BondGrandiINa_SBSERCA,400,50,params2,TG_x0,'ioninfo');

[T3,S3,loginfo3,x3,ioninfo3] = pacemodel(@Final_BondGrandiINa_SBSERCA,600,50,params2,TG_x0,'ioninfo');

[T4,S4,loginfo4,x4,ioninfo4] = pacemodel(@Final_BondGrandiINa_SBSERCA,800,50,params2,TG_x0,'ioninfo');

[T5,S5,loginfo5,x5,ioninfo5] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,50,params2,TG_x0,'ioninfo');

plot(T{50},S{50}(:,37),'r')
hold on 
plot(T1{50},S1{50}(:,37),'y')
plot(T2{50},S2{50}(:,37),'g')
plot(T3{50},S3{50}(:,37),'c')
plot(T4{50},S4{50}(:,37),'b')
plot(T5{9},S5{9}(:,37),'m')
plot(T5{35},S5{35}(:,37),'k')
% 
[peak_100,peak_100ind] = max(S{50}(:,37));
peak_100_90repol_ind = find(S{50}(peak_100ind:end,37)<(S{49}(end,37)+0.1*(peak_100-(S{49}(end,37)))),1,'first')+peak_100ind-1;
peak_100_APD_90 = T{50}(peak_100_90repol_ind)-T{50}(peak_100ind);

[peak_200,peak_200ind] = max(S1{50}(:,37));
peak_200_90repol_ind = find(S1{50}(peak_200ind:end,37)<(S1{49}(end,37)+0.1*(peak_200-(S1{49}(end,37)))),1,'first')+peak_200ind-1;
peak_200_APD_90 = T1{50}(peak_200_90repol_ind)-T1{50}(peak_200ind);

[peak_400,peak_400ind] = max(S2{50}(:,37));
peak_400_90repol_ind = find(S2{50}(peak_400ind:end,37)<(S2{49}(end,37)+0.1*(peak_400-(S2{49}(end,37)))),1,'first')+peak_400ind-1;
peak_400_APD_90 = T2{50}(peak_400_90repol_ind)-T2{50}(peak_400ind);

[peak_600,peak_600ind] = max(S3{50}(:,37));
peak_600_90repol_ind = find(S3{50}(peak_600ind:end,37)<(S3{49}(end,37)+0.1*(peak_600-(S3{49}(end,37)))),1,'first')+peak_600ind-1;
peak_600_APD_90 = T3{50}(peak_600_90repol_ind)-T3{50}(peak_600ind);

[peak_800,peak_800ind] = max(S4{50}(:,37));
peak_800_90repol_ind = find(S4{50}(peak_800ind:end,37)<(S4{49}(end,37)+0.1*(peak_800-(S4{49}(end,37)))),1,'first')+peak_800ind-1;
peak_800_APD_90 = T4{50}(peak_800_90repol_ind)-T4{50}(peak_800ind);

[peak_1000,peak_1000ind] = max(S5{50}(:,37));
peak_1000_90repol_ind = find(S5{50}(peak_1000ind:end,37)<(S5{49}(end,37)+0.1*(peak_1000-(S5{49}(end,37)))),1,'first')+peak_1000ind-1;
peak_1000_APD_90 = T5{50}(peak_1000_90repol_ind)-T5{50}(peak_1000ind);

[peak_1000NoEad,peak_1000NoEadind] = max(S5{9}(:,37));
peak_1000NoEad_90repol_ind = find(S5{9}(peak_1000NoEadind:end,37)<(S5{8}(end,37)+0.1*(peak_1000NoEad-(S5{8}(end,37)))),1,'first')+peak_1000NoEadind-1;
peak_1000NoEad_APD_90 = T5{9}(peak_1000NoEad_90repol_ind)-T5{9}(peak_1000NoEadind);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG+Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
