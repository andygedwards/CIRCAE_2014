function [dx,varargout]=Final_BondGrandiINa_SBSERCA(t,x,p,loginfo)
% Usage:
%BondInit
%[T,S]=ode15slog(@Bond,[0,60],x0,[],p);

%-------- Mouse model from Bondarenko et al. 2010, with noted modifications
%-------- Implemented by Andy Edwards and Kevin Vincent. January 2013      

%----------Initialize vector for Diff eqs. --------------------------

% -----INa-------
CNa2    = x(1) ;
CNa1    = x(2) ;
ONa     = x(3) ;
IFNa    = x(4) ;
I1Na    = x(5) ;
CNa3    = x(6);
ICNa2   = x(7) ;
ICNa3   = x(8) ;
LONa    = x(9);
LCNa1   = x(10);
LCNa2   = x(11);
LCNa3   = x(12);
% -------ICa------
O       = x(13) ;
C2      = x(14);
C3      = x(15);
C4      = x(16);
I1      = x(17);
I2      = x(18);
I3      = x(19);
% -------RyR------
PO     = x(20);
PRI     = x(21);
PI     = x(22);
PRyR    = x(23);
% --------IK------
atof   = x(24);
itof   = x(25);
atos   = x(26);
itos   = x(27);
nKs     = x(28);
% ------ Ca buffers -------
LTRPNCa = x(29);
HTRPNCa = x(30);
% ---------Cellular concentrations and membrane voltage -----
Nai   = x(31);
Ki    = x(32);
Cai   = x(33);
CaSS  = x(34);
CaJSR = x(35);
CaNSR = x(36);
V       = x(37);
% --------IK------
aur     = x(38);
iur     = x(39);
aKss    = x(40);
CK2    = x(41);
CK1    = x(42);
OK     = x(43);
IK     = x(44);
Caout  = x(45);
Naout  = x(46);
Kout   = x(47);

dx = zeros(length(x),1);

%----------end init diff eqs. -----------------------------------

%C-------COMPUTE REVERSAL POTENTIALS-------------------------------------------
% Nae=10;Nai=10;ENa=(p.RT_over_F)*log(Nae/Nai);      
ENa = p.RT_over_F*log((0.9*p.Nao+0.1*p.Ko)/(0.9*Nai+0.1*Ki)); 
EK =  p.RT_over_F*log(p.Ko/Ki);   
EKr = p.RT_over_F*log((0.98*p.Ko+0.02*p.Nao)/(0.98*Ki+0.02*Nai));
ECaN = 0.5*p.RT_over_F*log(p.Cao/Cai); %for ICab
ECaL = p.ECaL;
%--------- end reversal potentials --------------------------------------

%-------COMPUTE  Ion currents  ---------------------------

INa = p.GNa*(ONa+LONa)*(V-ENa);
ICa = p.GCa*O*(V-ECaL);
ICab = p.GCab*(V-ECaN);
IpCa = p.IpCamax*Cai^2/(p.KmpCa^2+Cai^2);
Jpca = IpCa*p.Acap/(p.Vmyo*p.Faraday)/2;
INab = p.GNab*(V-ENa);
IKtof = p.Gtof*atof^3*itof*(V-EK);
IKtos = p.Gtos*atos*itos*(V-EK);
IK1 = p.GK1*p.Ko/(p.Ko+210.0)*(V-EK)/(1+exp(0.0896*(V-EK)));
IKs = p.GKs*nKs^2*(V-EK);
IKur = p.GKur*aur*iur*(V-EK);
IKss = p.GKss*aKss*(V-EK);
IKr = p.GKr*OK*(V-EKr);

%---------IClCa-------------------------------------------
OClCa = 0.2/(1+exp(-(V-46.7)/7.8));
a = Cai/(Cai+p.KmCl);

ICl = p.GCl*OClCa*a*(V-p.ECl);

%------- INaK, INaCa,-------------------------------------

VF_over_RT=V/p.RT_over_F;

%--------INaK---------------------------------------------
sigma = (exp(p.Nao/67300)-1.0)/7.0;
a1 = 1+0.1245*exp(-0.1*VF_over_RT);
a2 = 0.0365*sigma*exp(-VF_over_RT);
fNaK = 1/(a1+a2);
a1 = p.Ko/(p.Ko+p.KmKo);
a2 = 1+(p.KmNai/Nai)^1.5;

INaK = (p.INaKmax*fNaK*a1)/a2;

%--------INaCa-------------------

% The following 3 lines are added to calculate the NCX current exposed
% subspace Ca: p.beta specifies the fraction of total NCX current existing
% in the subspace.

SS_a1 = (Nai^3*p.Cao*exp(p.eta*VF_over_RT)-p.Nao^3*CaSS*exp((p.eta-1)*VF_over_RT));%/(1+(0.125/CaSS)^2);
SS_a2 = p.Kmcao*Nai^3+p.Kmnao^3*CaSS+p.Knai^3*p.Cao*(1+CaSS/p.Kmcai)+p.Kmcai*p.Nao^3*(1+(Nai/p.Knai)^3)+Nai^3*p.Cao+p.Nao^3*CaSS;
SS_a3 = 1+p.ksat*exp((p.eta-1)*VF_over_RT);

% The following 3 lines are added to calculate the NCX current exposed cytosolic Ca

i_a1 = (Nai^3*p.Cao*exp(p.eta*VF_over_RT)-p.Nao^3*Cai*exp((p.eta-1)*VF_over_RT));%/(1+(0.125/Cai)^2);
i_a2 = p.Kmcao*Nai^3+p.Kmnao^3*Cai+p.Knai^3*p.Cao*(1+Cai/p.Kmcai)+p.Kmcai*p.Nao^3*(1+(Nai/p.Knai)^3)+Nai^3*p.Cao+p.Nao^3*Cai;
i_a3 = 1+p.ksat*exp((p.eta-1)*VF_over_RT);

%cytosolic INCX
i_INaCa = p.vNCXmax*i_a1*(1-p.beta)/(i_a2*i_a3);

%subspace INCX
SS_INaCa = p.vNCXmax*SS_a1*p.beta/(SS_a2*SS_a3);

INaCa = i_INaCa + SS_INaCa;

i_JNCX = -i_INaCa*p.Acap/(p.Vmyo*p.Faraday);
SS_JNCX = -SS_INaCa*p.Acap/(p.VSS*p.Faraday);

% below is the original version of the INaCa code prior to partitioning
% some NCX to the subspace.

% a1 = (Nai^3*p.Cao*exp(p.eta*VF_over_RT)-p.Nao^3*Cai*exp((p.eta-1)*VF_over_RT))/(1+(0.125/Cai)^2);
% a2 = 1400*Nai^3+88000^3*Cai+12000^3*p.Cao*(1+Cai/3.6)+3.6*p.Nao^3*(1+(Nai/12000)^3)+Nai^3*p.Cao+p.Nao^3*Cai;
% a3 = 1+p.ksat*exp((p.eta-1)*VF_over_RT);
% 
% 
% INaCa = p.vNCXmax*a1/(a2*a3);
% JNCX = -INaCa*p.Acap/(p.Vmyo*p.Faraday);

%-------end Ion currents -----------------------------

%--------------------initial stimulating current Ist-----------------------

if p.caffeine == 1 % if caffeine is being simulated
    dPRyR = 0.95*(1-PRyR);
    Ist = 0;
else  
    dPRyR = -p.dPRyR*PRyR-0.1*(ICa/p.ICaLmax)*exp(-V/30);
    Ist = 0;
    if (t >= p.pacestart && t < p.pacestart+p.pacedur)
      Ist = p.paceamp;
    end
end

%end-----------------initial stimulating current Ist----------------------

%-------------------- Voltage Clamp current -----------------------

IVClamp = 0;
if length(p.VClampAmp) > 0
  Ist = 0;
  ind = find(p.VClampTimes>t);
  if length(ind) > 0
    IVClamp = 1e3*(V-p.VClampAmp(ind(1)))/(p.VClampR*p.Acap);
  end
end

%end----------------- Voltage Clamp current ----------------------


%-------COMPUTE GATING VARIABLE DERIVATIVES-----------------------------------

%-----------------INa--------------------------------------------------

%--------------Transition rates-----------

alphaNa11 = p.P1a1/(p.P2a1*exp(-(V+p.P3a1)/p.P4a11)+ p.P5a11*exp(-(V+p.P3a1)/150));
alphaNa12 = p.P1a1/(p.P2a1*exp(-(V+p.P3a1)/p.P4a12)+ p.P5a12*exp(-(V+p.P3a1)/150));
alphaNa13 = p.P1a1/(p.P2a1*exp(-(V+p.P3a1)/p.P4a13)+ p.P5a13*exp(-(V+p.P3a1)/150));
betaNa11 = p.P1b1*exp(-(V+p.P3a1)/p.P2b1);
betaNa12 = p.P1b2*exp(-(V-p.P2b2)/p.P2b1);
betaNa13 = p.P1b3*exp(-(V-p.P2b3)/p.P2b1);
alphaNa4 = 1/(p.P1a4*exp(-(V+7)/p.P2a4)+p.P3a4);
alphaNa5 = p.P1a5*exp(-(V+7)/p.P2a5);
betaNa5 = p.P1b5 + p.P2b5*(V+7);
betaNa4 = (alphaNa13*alphaNa4*alphaNa5)/(betaNa13*betaNa5);
alphaNa6 = alphaNa4/p.P1a6;
betaNa6 = p.P1b6*exp(-V/p.P2b6);
alphaNa7 = p.P1a7*exp(V/p.P2a7);
betaNa7 = p.P1b7*exp(-V/p.P2b7);
alphaNa8 = p.P1a8;
betaNa8 = p.P1b8;

%--------end transition rates----------------

I2Na = (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3));
dCNa3 = betaNa8*LCNa3+betaNa11*CNa2+alphaNa5*ICNa3-(alphaNa11+betaNa5+alphaNa8)*CNa3;
dCNa2  = betaNa8*LCNa2+alphaNa11*CNa3+betaNa12*CNa1+alphaNa5*ICNa2-(betaNa11+alphaNa12+betaNa5+alphaNa8)*CNa2; 
dCNa1  = betaNa8*LCNa1+alphaNa12*CNa2+betaNa13*ONa+alphaNa5*IFNa-(betaNa12+alphaNa13+betaNa5+alphaNa8)*CNa1;
dONa   = betaNa8*LONa+alphaNa13*CNa1+betaNa4*IFNa-(betaNa13+alphaNa4+alphaNa8)*ONa;
dIFNa  = alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa12*ICNa2-(betaNa4+alphaNa5+alphaNa6+betaNa12)*IFNa;
dI1Na  = alphaNa6*IFNa+betaNa7*I2Na-(betaNa6+alphaNa7)*I1Na;
dICNa2 = alphaNa11*ICNa3+betaNa12*IFNa+betaNa5*CNa2-(betaNa11+alphaNa12+alphaNa5)*ICNa2;
dICNa3 = betaNa11*ICNa2+betaNa5*CNa3-(alphaNa11+alphaNa5)*ICNa3;
dLONa = alphaNa13*LCNa1+alphaNa8*ONa-(betaNa8+betaNa13)*LONa;    
dLCNa1 = alphaNa8*CNa1+alphaNa12*LCNa2+betaNa13*LONa-(betaNa8+betaNa12+alphaNa13)*LCNa1;   
dLCNa2 = betaNa12*LCNa1+alphaNa8*CNa2+alphaNa11*LCNa3-(betaNa8+betaNa11+alphaNa12)*LCNa2;   
dLCNa3 = alphaNa8*CNa3+betaNa11*LCNa2-(betaNa8+alphaNa11)*LCNa3;  

%--------------end INa --------------------------------------------------

%--------------ICa --------------------------------------------------

%------------Transition rates---------------

alpha = 0.4*exp((V+p.alphaV)/15.0); % Bondarenko 2010
beta = 0.13*exp(-(V+p.betaV)/18.0); % Bondarenko 2010

gamma = (p.Kpcmax*CaSS)/(p.Kpchalf+CaSS);

Kpcf = 2.5;

C1 = (1-(O+C2+C3+C4+I1+I2+I3));

dO  = alpha*C4+p.Kpcb*I1+0.001*(alpha*I2-Kpcf*O)-O*(4*beta+gamma);
dC2 = 4*alpha*C1+2*beta*C3-(beta+3*alpha)*C2;
dC3 = 3*alpha*C2+3*beta*C4-(2*beta+2*alpha)*C3;
dC4 = 2*alpha*C3+4*beta*O+0.01*(4*p.Kpcb*beta*I1-alpha*gamma*C4)+0.002*(4*beta*I2-Kpcf*C4)+4*beta*p.Kpcb*I3-(3*beta*C4+alpha*C4+gamma*Kpcf*C4);
dI1 = gamma*O+0.001*(alpha*I3-Kpcf*I1)+0.01*(alpha*gamma*C4-4*beta*p.Kpcb*I1)-p.Kpcb*I1;
dI2 = 0.001*(Kpcf*O-alpha*I2)+p.Kpcb*I3+0.002*(Kpcf*C4-4*beta*I2)-gamma*I2;
dI3 = 0.001*(Kpcf*I1-alpha*I3)+gamma*I2+gamma*Kpcf*C4-(4*beta*p.Kpcb*I3+p.Kpcb*I3);

%------- end ICa ----------------

%-------COMPUTE K state variables--------------

%----------IKs---------------------
alpha_n = 0.00000481333*(V+26.5)/(1-exp(-0.128*(V+26.5)));
beta_n = 0.0000953333*exp(-0.038*(V+26.5));


%----------Steady state variables for IKur, IKtof, IKtos and IKss--------
% % Bondarenko 2004
ass = 1.0/(1+exp(-(V+22.5)/7.7));
iss = 1.0/(1+exp((V+45.2)/5.7)); 
% % Bondarenko 2004

% Li et al 2010
% iss = 1.0/(1+exp((V+45.2)/5.7));
% Li et al 2010

% %----------IKur-------------------
% 
% 
% % Bondarenko 2004
tau_iur = p.tauiur-170*iss; 
tau_aur = 0.493*exp(-0.0629*V)+2.058;
% % Bondarenko 2004

% Li et al. 2010
% iurss = 1.0/(1+exp((V+42.1)/5.4));
% ass = 1/(1+exp(-(V+6.19)/9.6));
% 
% tau_iur = 643+1000*iurss;
% tau_aur = 0.493*exp(-0.0629*V)+2.058;
% Li et al. 2010

%34 = aur
%35 = iur

%---------IKtof------------------

% % Bondarenko 2004
alpha_i = (0.000152*exp(-(V+13.5)/7))/(0.067083*(exp(-(V+33.5)/7))+1); % 
beta_i = (0.00095*exp((V+33.5)/7))/(0.051335*exp((V+33.5)/7)+1);
alpha_a = 0.18064*exp(0.03577*(V+30));
beta_a = 0.3956*exp(-0.06237*(V+30));
% % Bondarenko 2004

% Li et al. 2010
% itoss = 1.0/(1+exp((V+51.4)/5));
% tau_i = 9.66+(10.94/(1+exp((V+54.1)/5)));
% 
% alpha_a = 0.18064*exp(0.03577*(V+45));
% beta_a = 0.3956*exp(-0.06237*(V+45));
% Li et al. 2010

%-----------IKtos-------------------
c = (1+exp((p.tisign*(V+p.tinum2))/p.tidenom));
tau_ti = p.ticoeff*exp(((-(V+45)/240)^2)*p.strain)+p.tinum1/c;
tau_ta = 0.493*exp(-0.0629*V)+2.058;

%--------------IKss-------------------
tau_Kss = 39.3*exp(-0.0862*V)+13.17;

%--------------IKr--------------------
alpha_a0 = 0.022348*exp(0.01176*V);
beta_a0 = 0.047002*exp(-0.0631*V);
% Bondarenko 2004
alpha_a1 = 0.013733*exp(0.038198*V);
beta_a1 = 0.0000689*exp(-0.04178*V);
% Bondarenko 2004

% Li et al 2010
% alpha_a1 = 0.0335*exp(0.0109*V);
% beta_a1 = 0.0703*exp(0.0287*(V+5));
% Lit et al 2010

alpha_iKr = 0.090821*exp(0.023391*(V+5));
beta_iKr = 0.006497*exp(-0.03268*(V+5));

CK0 = (1-(CK1+CK2+OK+IK));

dCK2= ((p.kf*CK1+beta_a1*OK)-(p.kb*CK2+alpha_a1*CK2));
dCK1= ((alpha_a0*CK0+p.kb*CK2)-(beta_a0*CK1+p.kf*CK1));
dOK = ((alpha_a1*CK2+beta_iKr*IK)-(beta_a1*OK+alpha_iKr*OK));
dIK = (alpha_iKr*OK-beta_iKr*IK);


%---------end K state variables--------------------------------------------

%----------COMPUTE RyR state derivatives-----------------------------------

kCaSR=p.MaxSR-(p.MaxSR-p.MinSR)/(1+(p.EC50SR/CaJSR)^2.5);
koSRCa=p.koCa/kCaSR;
kiSRCa=p.kiCa*kCaSR;

PR = (1.0-(PO+PRI+PI));
dPO= (koSRCa*(CaSS)^2*PR-p.kom*PO)-(kiSRCa*CaSS*PO-p.kim*PI);
dPRI= (p.kom*PI-koSRCa*(CaSS)^2*PRI)-(p.kim*PRI-kiSRCa*CaSS*PR);
dPI= (kiSRCa*CaSS*PO-p.kim*PI)-(p.kom*PI-koSRCa*(CaSS)^2*PRI);

%-------COMPUTE Intracellular calcium fluxes Jup, Jleak, Jrel, Jtr, Jtrpn, Jxfer-------

Jup = p.Q10SERCA*((p.SERCAVmax*(Cai/p.kmf)^p.HSERCA)-(p.SERCAVmax*(CaNSR/p.kmr)^p.HSERCA))/(1+((Cai/p.kmf)^p.HSERCA)+((CaNSR/p.kmr)^p.HSERCA));
Jleak = p.kleak*(CaNSR-Cai);
Jrel = p.Grel*(PO)*(CaJSR-CaSS)*PRyR; 
Jtr = (CaNSR - CaJSR)/p.tautr;
Jxfer = (CaSS-Cai)/p.tauxfer;
Jtrpn = p.khtrpn_plus*Cai*(p.HTRPNtot-HTRPNCa)+p.kltrpn_plus*Cai*(p.LTRPNtot-LTRPNCa)-(p.khtrpn_minus*HTRPNCa+p.kltrpn_minus*LTRPNCa);

%-------COMPUTE buffer scale factors--------------------------------

a1 = p.CMDNtot*p.KmCMDN/((CaSS+p.KmCMDN)^2);               
%a2 = p.EGTAtot*p.KmEGTA/((CaSS+p.KmEGTA)^2);                                   
betaSS = 1/(1+a1);
CaCMDNSS = p.CMDNtot*CaSS/(CaSS+p.KmCMDN);
		
a1 = p.CSQNtot*p.KmCSQN/((CaJSR+p.KmCSQN)^2);             
betaJSR = 1/(1+a1);                                   
CaCSQNJSR = p.CSQNtot*CaJSR/(CaJSR+p.KmCSQN);

a1 = p.CMDNtot*p.KmCMDN/((Cai+p.KmCMDN)^2);
%a2 = p.EGTAtot*p.KmEGTA/((Cai+p.KmEGTA)^2);
betai = 1/(1+a1);
CaCMDNi = p.CMDNtot*Cai/(Cai+p.KmCMDN);

%-------COMPUTE CONCENTRATION DERIVATIVES-------------------------
    
fmyo = p.Acap/(p.Vmyo*p.Faraday);
fss = p.Acap/(2*p.VSS*p.Faraday);

% ---------COMPUTE VOLTAGE DERIVATIVES ------------------

a1 = INa+INab+INaK+INaCa;                      %Na-currents
a2 = IKtof+IKtos+IK1+IKs+IKss+IKur+IKr;        %K- currents
a3 = IpCa+ICa+ICab;                            %Ca-currents;

Itot =  a1+a2+a3+ICl+Ist+IVClamp;
Jtotal = Jup + SS_JNCX*p.VSS/p.Vmyo +i_JNCX + Jpca;

%--------END MAIN-----------------------------------------------------------

%----------------Transient Sodium Current-------------------

dx(1) = dCNa2;
dx(2) = dCNa1;
dx(3) = dONa;  
dx(4) = dIFNa; 
dx(5) = dI1Na; 
dx(6) = dCNa3;
dx(7) = dICNa2;
dx(8) = dICNa3;
dx(9) = dLONa;    
dx(10) = dLCNa1;   
dx(11) = dLCNa2;   
dx(12) = dLCNa3;   

%----------------L-type Ca Current-------------------------
dx(13)  = dO; 
dx(14) = dC2; 
dx(15) = dC3; 
dx(16) = dC4; 
dx(17) = dI1; 
dx(18) = dI2; 
dx(19) = dI3; 

%----------------Ryanodine Receptors-------------------------
dx(20) = dPO;
dx(21) = dPRI;
dx(22) = dPI;
dx(23) = dPRyR;

%----------------Misc. Potassium Currents-------------------------
% dx(24) = alpha_a*(1-atof)-(beta_a*atof);  %atof (Li et al. 2010 paper)
% dx(25) = (itoss-itof)/tau_i;            %itof (Li et al. 2010 paper)
dx(24) = alpha_a*(1-atof)-(beta_a*atof); % alternative atof (Bondarenko 2004 and Li et al 2010 cellML code)
dx(25) = alpha_i*(1-itof)-(beta_i*itof); % alternative itof (Bondarenko 2004 and Li et al 2010 cellML code)
dx(26) = (ass-atos)/tau_ta;               %atos
dx(27) = (iss-itos)/tau_ti;               %itos
dx(28) = alpha_n*(1-nKs)-beta_n*nKs;       %nKs  
dx(38) = (ass-aur)/tau_aur;                %aur  (Bondarenko 2004)
dx(39) = (iss-iur)/tau_iur;              %iur  (Bondarenko 2004)
% dx(38) = (ass-aur)/tau_aur;                %aur  (Li et al. 2010 paper)
% dx(39) = (iurss-iur)/tau_iur;                %iur  (Li et al. 2010 paper)
dx(40) = (ass-aKss)/tau_Kss;               %aKss 

%-------Troponin buffer ---------------------------------
dx(29) = p.kltrpn_plus*Cai*(p.LTRPNtot - LTRPNCa) - p.kltrpn_minus * LTRPNCa;        
dx(30) = p.khtrpn_plus*Cai*(p.HTRPNtot - HTRPNCa) - p.khtrpn_minus * HTRPNCa;        

%------Ion concentrations
dx(31) = -(INa+INab+3*INaK+3*INaCa)*fmyo;	                     %[Na]i
dx(32) = -(IKtof+IKtos+IK1+IKs+IKss+IKur+IKr-2*INaK+Ist)*fmyo;           %[K]i
dx(33) = betai*(Jleak+Jxfer-Jup-Jtrpn-(ICab-2*i_INaCa+IpCa)*fmyo/2);  %[Ca]i
dx(34) = betaSS*(Jrel*p.VJSR/p.VSS - Jxfer*p.Vmyo/p.VSS - ICa*fss + SS_INaCa*p.Acap/(p.VSS*p.Faraday)); %[Ca]ss
dx(35) = betaJSR*(Jtr - Jrel);	                                     %[Ca]JSR
dx(36) = (Jup-Jleak)*p.Vmyo/p.VNSR - Jtr*p.VJSR/p.VNSR;	             %[Ca]nsr

%----------Membrane voltage------------------
dx(37) = -Itot;

%----------Delayed rectifier current IKr------------
dx(41) = dCK2;
dx(42) = dCK1;
dx(43) = dOK;
dx(44) = dIK;
dx(45) = (ICab-2*INaCa+IpCa+ICa)*p.Acap/(2*p.Faraday);
dx(46) = (INa+INab+3*INaK+3*INaCa)*p.Acap/p.Faraday;
dx(47) = (IKtof+IKtos+IK1+IKs+IKss+IKur+IKr-2*INaK)*p.Acap/p.Faraday;

% Logging part
nout = max(nargout,1)-1; % Number of output parameters in excess of dx
if nout > 0
  if nargin == 4 && length(loginfo)> 0
    logvar = fieldnames(loginfo);
    nlog = length(logvar);
    argout = cell(1,nlog);
    for i=1:nlog
      switch logvar{i}
       case 'state_der'
        argout{i} = dx;
       case 'currents'
        argout{i} = [INa,INab,IKtof,IKtos,IK1,IKs,IKss,IKur,IKr,IpCa,ICa,...
                     ICl,ICab,INaCa,INaK,IVClamp];
       case 'ionfluxes'
        argout{i} = [Jxfer,Jtr,Jtotal,Jpca,Jup,SS_JNCX,i_JNCX,Jrel];
       case 'buffers'
        argout{i} = [CaCMDNSS, CaCMDNi, CaCSQNJSR, LTRPNCa, HTRPNCa];
       case 'rates'
        argout{i} = [alphaNa11,alphaNa12,alphaNa13,betaNa11,betaNa12,...
                     betaNa13,alphaNa3,betaNa3,alphaNa2,betaNa2,alphaNa4,...
                     betaNa4,alphaNa5,betaNa5,alpha,beta,gamma,Kpcf];
%       case 'hiddenstates'
%        argout{i} = [CNa3,C1,CK0];
      end
    end
    varargout{1} = argout;
  else
    varargout{1} = [];
  end
end

