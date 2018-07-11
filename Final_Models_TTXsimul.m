load('TG_x0');
final_Iso_ind = 36;

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
%*********
params2.GNa = 13; % normal conductance 
%********* 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG+Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[T3,S3,loginfo3,x3,ioninfo3] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,final_Iso_ind,params2,TG_x0,'ioninfo');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TG+Iso %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
TG_ablate_INa_x0 = create_inits(S3,final_Iso_ind);

params3.beta = 0.30;
%******************
params3.vNCXmax=2.0*2.6;
%******************
params3.SERCAVmax = 0.2*0.7;
params3.kmr = 3500*0.55;
params3.kmf = 0.6*0.55;
%******************
params3.GCa = 0.1729*1.2*2.0; % Iso GCa increase
%******************
params3.SRsens = 300;
params3.KmNai = 16650*0.72; % Despa et al., Circ Res 2005
params3.Kpcmax = 0.11662*0.75;% Change in Bondarenko 2010 + match TGthalf (maier & kohlhaas)
params3.koCa = 0.005;
params3.kiCa = 0.22;
params3.kim = 1.2;
params3.kom = 0.06;
params3.MaxSR = 10;
params3.MinSR = 1;
params3.kleak = (0.1e-5)*5; % CaMKII leak is much higher (Maier 2003 for this value - spark frequency)
%******************* Bondarenko et al 2004
params3.Gtof = 0.4067*0.65;   
params3.GKss = 0.05;   
params3.GKur = 0.160*1.20;
%******************* Bondarenko et al 2004
params3.GK1  = 0.2938*0.5; % Just the Tg effect here
params3.dPRyR = 0.02;
params3.P2a5 = 7.8;
params3.P1b5 = 0.0084/1.8;
params3.P1b6 = 6.5449e-07;
params3.P1a7 = 3.377e-04;
params3.P1b7 = 1.868e-04;
params3.P1a8 = 6.5e-06/1000;
params3.P1b8 = 3.8e-03/1000;
params3.alphaV = 22;
params3.betaV = 22;  
%******************* Simulated 1 uM TTX
params3.GNa=13*0.5;
%******************* Simulated 1 uM TTX

[T4,S4,loginfo4,x4,ioninfo4] = pacemodel(@Final_BondGrandiINa_SBSERCA,1000,10,params3,TG_ablate_INa_x0,'ioninfo');

T(1:length(T3)) = T3;
T(end+1:end+length(T4)) = T4;
S(1:length(S3)) = S3;
S(end+1:end+length(S4)) = S4;

preintroT = [0,200];
preintroV = [S{length(S3)}(1,37),S{length(S3)}(1,37)];
postintroT = [0,200];
postintroV = [S{length(S3)+1}(1,37),S{length(S3)+1}(1,37)];

plot(preintroT,preintroV,'r');
hold on 
plot(postintroT,postintroV,'b');

hold on 
plot(T{length(T3)}+200,S{length(S3)}(:,37),'r')
plot(T{length(T3)+1}+200,S{length(S3)+1}(:,37),'b')

cat_T(1) = 1;
cat_T(2) = 200;
uno = struct2cell(TG_x0);
for i = 1:length(uno)
    cat_S(1,i) = uno{i};
end
cat_S(2,:) = cat_S(1,:);

for i = 1:length(T)
    cat_T(end+1:end+length(T{i})) = T{i}+cat_T(end);
    cat_S(length(cat_S(:,1))+1:length(S{i}(:,1))+length(cat_S(:,1)),:) = S{i};
end
    


