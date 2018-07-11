%%%%%%%%%%%%%%%%%%%%%%%%%%%% TG + Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load('TG_ablateINa_x0');

 
params2.beta = 0.30;
%******************
params2.vNCXmax=2.0*2.5;
%******************
params2.SERCAVmax = 0.2*0.7;
params2.kmr = 3500*0.55;
params2.kmf = 0.6*0.55;

%******************
params2.GCa = 0.1729*1.2*2.0; % Iso GCa increase WITHOUT ADDITIONAL TG SPECIFIC INCREASE
%******************

params2.SRsens = 300;
params2.KmNai = 16650*0.72; % Despa et al., Circ Res 2005

%******************
params2.Kpcmax = 0.11662*0.75;% Change in Bondarenko 2010 % Iso GCa increase WITHOUT ADDITIONAL TG SPECIFIC INCREASE
%******************

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
params2.P2a5 = 7.2; % WT value
params2.P1b5 = 0.0084/1.9; % WT value
params2.P1b6 = 8.9554e-07; % WT value
params2.P1a7 = 4.8696e-05; % WT value
params2.P1b7 = 2.868e-04; % WT value
params2.P1a8 = 1e-08/1000; % WT value
params2.P1b8 = 9.8e-03/1000; % WT value
params2.alphaV = 22;
params2.betaV = 22;

% 
[T4,S4,loginfo4,x4,ioninfo4] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,1,params2,TG_ablateINa_x0,'ioninfo');