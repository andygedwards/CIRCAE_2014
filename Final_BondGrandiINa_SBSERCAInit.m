function [p, x0, varargout]  = Final_BondGrandiINa_SBSERCAInit(varargin) %loginfo, names,
                                                  %currents
 
% Initial conditions and parameters 
% to Mouse model from Bondarenko et al. 2004
% Implemented by Andy Edwards and Kevin Vincent. January 2013   

  if nargout < 2
    error('To few output arguments.')
  end
  nxargout = max(nargout-2,0);
  
%---------PHYSICAL CONSTANTS----------------------------------------------------
	
p.Faraday=96.485;
p.Temp=310;
p.Rgas=8.314;
p.RT_over_F=p.Rgas*p.Temp/p.Faraday;
p.Avogadros = 6.022e23;  
%C-------CELL GEOMETRY PARAMETERS----------------------------------------------
  
%--The assumption of 1uF/cm^2 holds for this model
%--therefore Acap in cm^2 is equal to whole cell
%--capacitance in uF, i.e Cm=1

p.Acap=1.534e-4;
p.Vmyo=25.84e-6;
p.VJSR=0.12e-6;
p.VNSR=2.098e-6;
p.VSS=1.485e-9;

%C-------STANDARD IONIC CONCENTRATIONS-----------------------------------------
%--------Concentrations in uM not mM ------------------------------------------
p.Ko=5400;
p.Nao=140000;
p.Cao=1000;

%C-------MEMBRANE CURRENT PARAMETERS-------------------------------------------
% K+ currents
p.GKr=0.078;
p.GKs= 0.00575;
p.Gtos = 0.0;  
p.Gtof = 0.4067; % Bondarenko 2004
p.GKss = 0.05656;   % Bondarenko 2004
p.GKur = 0.16;   % Bondarenko 2004
p.GK1  = 0.2938;

p.kf = 0.023761;
p.kb = 0.036778;
p.tauiur = 1200; % Bondarenko 2004
% end K+ currents

% NCX
p.beta = 0.3;
p.vNCXmax = 2.0; 
p.Kmcao = 1300;
p.Kmcai = 3.59;
p.Kmnao = 87500;
p.Knai = 12290;
p.ksat = 0.27;
p.eta = 0.35;
% NCX

% Na-K ATPase
p.INaKmax = 0.88;
p.KmNai = 16650; % Altered to achieve stable and realistic Nai - See Despa et al. 2005
p.KmKo = 1500;
% Na-K ATPase

% Sarcolemmal Ca2+ ATPase
p.IpCamax = 0.085;   % Drastically reduced to achieve 10% flux contribution in WT model
p.KmpCa = 0.5;
% Sarcolemmal Ca2+ ATPase

% Background currents
p.GCab = 0.000165; 
p.GNab = 0.0018;
% Background currents

% Ca2+ activated Cl current
p.GCl = 10;   
p.KmCl = 10;
% Ca2+ activated Cl current
p.ECl = -40;  %Fixed reversal potential (RT/(-F)*ln([Cl]o/[Cl]i))


%C-------SR PARAMETERS---------------------------------------------------------

p.caffeine = 0;
p.Grel = 25;
p.EC50SR = 600;
p.MaxSR = 15;
p.MinSR = 1;
p.dPRyR = 0.04;
p.koCa = 10;
p.kiCa = 0.5;
p.kom = 0.06;
p.kim = 0.005;
p.RyRCaLdep = 0.1;
p.kleak = 0.1e-5;   %SR Ca2+ leak constant
p.kmr = 2100;
p.kmf = 0.300;
p.Q10SERCA = 2.6;
p.HSERCA = 1.787;
p.SERCAVmax = 0.2;   
p.tautr=1;
p.tauxfer=8;
p.kaplus=0.006075;
p.kaminus=0.07125;
p.kbplus=0.00405;
p.kbminus=0.965;
p.kcplus=0.009;
p.kcminus=0.0008;
p.n=4;
p.m=3;

%-------L-TYPE CALCIUM CHANNEL PARAMETERS-------------------------------------

p.GCa = 0.1729;
p.ECaL = 63;        %Fixed reversal potential (RT/2F*ln([Ca]o/[Ca]i))
p.Kpcmax = 0.23324;  
p.Kpchalf = 10; % Bondarenko 2010
p.Kpcb = 0.0005;
p.ICaLmax = 7.0;
p.alphaV = 15;
p.betaV = 15;

%-------Ito parameters-----------------------------------------------------

% Ito,f

p.icoeff = 40;
p.inum = 54.1;
p.idenom  = 5;

% Ito,s

p.tisign = 1;
p.tinum2 = 35.2;
p.tidenom = 10.7;
p.tinum1 = 2000;
p.ticoeff = 200;
p.strain = 0;

%-------BUFFER PARAMETERS--------------------------------------------------

%--------- Time dependent buffers-----------

p.LTRPNtot=70;
p.HTRPNtot=140;
p.khtrpn_plus=0.00237;
p.khtrpn_minus=3.2e-5;
p.kltrpn_plus=0.0327;
p.kltrpn_minus=0.0196;

%-------Fast buffers------------

p.CMDNtot=50.0;
p.CSQNtot=15000;
p.KmCMDN=0.238;
p.KmCSQN=800;

p.EGTAtot=0.0;
p.KmEGTA= 0.0;

%--------- Pace current parameters
p.pacestart = 0.0;
p.pacedur = 5.0;
p.paceamp = -8.0;

%--------- Voltage Clamp current

p.VClampAmp = [];
p.VClampTimes = [];
p.VClampR = 50; % In Mohm

%------------ INa
p.GNa = 13;
p.P1a1 = 3.802; p.P2a1 = 0.1027;p.P3a1 = 2.5;p.P4a11 = 17;p.P4a12 = 15;p.P4a13 = 12;p.P5a11 = 0.2;p.P5a12 = 0.23;p.P5a13 = 0.25;
p.P1b1 = 0.1917; p.P2b1 = 20.3;
p.P1b2 = 0.2; p.P2b2 = 2.5;
p.P1b3 = 0.22; p.P2b3 = 7.5;
p.P1a4 = 0.188495; p.P2a4 = 16.6; p.P3a4 = 0.393956;
p.P1a5 = 7e-07; p.P2a5 = 7.2; 
p.P1b5 = 0.0084/1.9; p.P2b5 = 2e-05; % Grandi's value for P1b5 is Bondarenko/1.9
p.P1a6 = 100; 
p.P1b6 = 8.9554e-07; p.P2b6 = 11.3944; 
p.P1a7 = 4.8696e-05; p.P2a7 = 23.2696; 
p.P1b7 = 2.868e-04; p.P2b7 = 35.9898;
p.P1a8 = 1e-08/1000;
p.P1b8 = 9.8e-03/1000;


%--------- end parameters ------------------------------

%---------- Initial conditions for sodium channel --------------------

u = -78.3418;

if nargin > 0 
  argin = cell2cell(varargin);
  for i = 1:length(argin) 
    if isstruct(argin{i})
      fields = fieldnames(argin{i})';
      % Copy the fields in the struct to inputstruct
      for field =fields
        inputstruct.(field{1}) = argin{i}.(field{1});
      end
      if length(inputstruct)>0
        inputfields = fieldnames(inputstruct)';
        for field=inputfields
          field = field{1};
          if isfield(p,field)
            p.(field) = inputstruct.(field);
          end
        end
      end
    end
  end 
end

alphaNa11 = p.P1a1/(p.P2a1*exp(-(u+p.P3a1)/p.P4a11)+ p.P5a11*exp(-(u+p.P3a1)/150));
alphaNa12 = p.P1a1/(p.P2a1*exp(-(u+p.P3a1)/p.P4a12)+ p.P5a12*exp(-(u+p.P3a1)/150));
alphaNa13 = p.P1a1/(p.P2a1*exp(-(u+p.P3a1)/p.P4a13)+ p.P5a13*exp(-(u+p.P3a1)/150));
betaNa11 = p.P1b1*exp(-(u+p.P3a1)/p.P2b1);
betaNa12 = p.P1b2*exp(-(u-p.P2b2)/p.P2b1);
betaNa13 = p.P1b3*exp(-(u-p.P2b3)/p.P2b1);
alphaNa4 = 1/(p.P1a4*exp(-(u+7)/p.P2a4)+p.P3a4);
alphaNa5 = p.P1a5*exp(-(u+7)/p.P2a5);
betaNa5 = p.P1b5 + p.P2b5*(u+7);
betaNa4 = (alphaNa13*alphaNa4*alphaNa5)/(betaNa13*betaNa5);
alphaNa6 = alphaNa4/p.P1a6;
betaNa6 = p.P1b6*exp(-u/p.P2b6);
alphaNa7 = p.P1a7*exp(u/p.P2a7);
betaNa7 = p.P1b7*exp(-u/p.P2b7);
alphaNa8 = p.P1a8;
betaNa8 = p.P1b8;

A=[-(alphaNa11+alphaNa5) betaNa11 0 0 0 betaNa5 0 0 0 0 0 0 0;...
    alphaNa11 -(betaNa11+alphaNa5+alphaNa12) betaNa12 0 0 0 betaNa5 0 0 0 0 0 0;...
        0 alphaNa12 -(betaNa12+betaNa4+alphaNa5+alphaNa6) betaNa6 0 0 0 betaNa5 alphaNa4 0 0 0 0;...
        0 0 alphaNa6 -(betaNa6+alphaNa7) betaNa7 0 0 0 0 0 0 0 0; ...
        0 0 0 alphaNa7 -betaNa7 0 0 0 0 0 0 0 0;...
        alphaNa5 0 0 0 0 -(betaNa5+alphaNa11+alphaNa8) betaNa11 0 0 betaNa8 0 0 0; ...
        0 alphaNa5 0 0 0 alphaNa11 -(betaNa11+alphaNa8+alphaNa12+betaNa5) betaNa12 0 0 betaNa8 0 0;...
        0 0 alphaNa5 0 0 0 alphaNa12 -(betaNa12+alphaNa8+alphaNa13+betaNa5) betaNa13 0 0 betaNa8 0; ...
        0 0 betaNa4 0 0 0 0 alphaNa13 -(alphaNa4+betaNa13+alphaNa8) 0 0 0 betaNa8;...
        0 0 0 0 0 alphaNa8 0 0 0 -(alphaNa11+betaNa8) betaNa11 0 0; ...
        0 0 0 0 0 0 alphaNa8 0 0 alphaNa11 -(betaNa8+betaNa11+alphaNa12) betaNa12 0;...
        0 0 0 0 0 0 0 alphaNa8 0 0 alphaNa12 -(betaNa8+betaNa12+alphaNa13) betaNa13; ...
        1 1 1 1 1 1 1 1 1 1 1 1 1];

b = zeros(13,1);
b(13) = 1;

x = linsolve(A, b);
  
%---------- Initial conditions----------------------------------------

x0.CNa2 = x(7);               % CNa2      (1)
x0.CNa1 = x(8);               % CNa1      (2)
x0.ONa = x(9);                % ONa       (3)
x0.IFNa = x(3);               % IFNa      (4)
x0.I1Na = x(4);               % I1Na      (5)
x0.CNa3 = x(6);               % CNa3      (6)
x0.ICNa2 = x(2);              % ICNa2     (7)
x0.ICNa3 = x(1);              % ICNa3     (8)
x0.LONa = x(13);              % LONa      (9)
x0.LCNa1 = x(12);             % LCNa1     (10)
x0.LCNa2 = x(11);             % LCNa2     (11)
x0.LCNa3 = x(10);             % LCNa3     (12)
x0.O = 9.30308e-19;           % O         (13)
x0.C2 = 0.124216e-3;          % C2        (14)
x0.C3 = 0.578679e-8;          % C3        (15)
x0.C4 = 0.119816e-12;         % C4        (16)
x0.I1 = 0.497923e-18;         % I1        (17)
x0.I2 = 0.345847e-13;         % I2        (18)
x0.I3 = 0.185106e-13;         % I3        (19)
x0.PO = 0.149102e-4;          % PO1(RyR)  (20)
x0.PRI = 0.951726e-10;        % PO2(RyR)  (21)
x0.PI = 0.167740e-3;          % PC2(RyR)  (22)
x0.PRyR = 0.0;                % PRyR      (23)
x0.atof = 0.265563e-2;        % atof      (24)
x0.itof = 0.999977;           % itof      (25)
x0.atos = 0.417069e-3;        % atos      (26)
x0.itos = 0.998543;           % itos      (27)
x0.nks = 0.262753e-3;         % nKs       (28)
x0.LTRPNCa = 11.2684;         % LTRPNCa   (29)
x0.HTRPNCa = 125.290;         % HTRPNCa   (30)
x0.Nai = 13000;               % Nai       (31)
x0.Ki  = 143000;              % Ki        (32)
x0.Cai = 0.115001;            % Cai       (33)
x0.CaSS = 0.115001;           % CaSS      (34)
x0.CaJSR = 1200;              % CaJSR     (35)
x0.CaNSR = 1200;              % CaNSR     (36)
x0.V = -78.3418;              % V         (37)
x0.aur = 0.417069e-3;         % aur       (38)
x0.iur = 0.998543;            % iur       (39)
x0.aKss = 0.417069e-3;        % aKss      (40)
x0.CK2 = 0.641229e-3;         % CK2       (41)
x0.CK1 = 0.992513e-3;         % CK1       (42)
x0.OK = 0.175298e-3;          % OK        (43)
x0.IK = 0.319129e-4;          % IK        (44)
x0.Caout = 0.0;               % Caout     (45)
x0.Naout = 0.0;               % Naout     (46)
x0.Kout = 0.0;                % Kout      (47)

if exist('inputstruct','var')
    if length(inputstruct)>0
        for field=inputfields
          field = field{1};
          if isfield(x0,field)
            x0.(field) = inputstruct.(field);
          end
        end
    end
end

%---------- End Initial conditions----------------------------------------

%-------------- Logging parameters--------------------------------
names.modelname = 'Bondarenko mouse model, with SB SERCA, RyR, and NCX';
%names.states = {'CNa2','CNa1','ONa','IFNa','I1Na','I2Na','ICNa2', ...
%                    'ICNa3','O','C2','C3','C4','I1','I2','I3','PO1','PO2', ...
%                    'PC2','PRyR','atof','itof','atos','itos','nKs', ...
%                    'LTRPNCa','HTRPNCa','Nai','Ki','Cai','CaSS', ...
%                    'CaJSR','CaNSR','V','aur','iur','aKss','CK2','CK1',...
%                    'OK','IK','Caout','Naout','Kout'};

names.states = fieldnames(x0)';
% Place the ionexchangers last, ie 'INaCa','INaK'. (For logging purposes)
names.currents = {'INa','INab','IKtof','IKtos','IK1', ...
                 'IKs','IKss','IKur','IKr','IpCa','ICa',...
                 'ICl','ICab','INaCa','INaK','IVClamp'};

names.ionfluxes = {'Jxfer','Jtr','Jtotal','Jpca','Jup','SS_JNCX','i_JNCX','Jrel'};
names.buffers = {'CaCMDNSS','CaCMDNi','CaCSQNJSR','LTRPNCa','HTRPNCa'};
names.rates = {'alphaNa11','alphaNa12','alphaNa13','betaNa11','betaNa12',...
               'betaNa13','alphaNa4','betaNa4','alphaNa5','betaNa5','alphaNa6',...
               'betaNa6','alphaNa7','betaNa7','alphaNa8','betaNa8','alpha','beta','gamma', ...
               'Kpcf'};
names.hiddenstates = {'CNa3','C1','CK0'};

% The possible logging variables. 
possible_logvar = {'state_der','currents','ionfluxes','buffers','rates'};%,'hiddenstates'};



loginfo = [];
if nargin > 0 
  % If ioninfo to be logged then add both buffers and currents as logging
  % parameters
  if ismemb('ioninfo',argin)
    if nxargout < 3
      error('If logging ioninfo, then 5 or more output variables are needed.')
    end
    if ~ismemb('buffers',argin)
      argin = [argin {'buffers'}];
    end
    if ~ismemb('currents',argin)
      argin = [argin {'currents'}];
    end
    if ~ismemb('ionfluxes',argin)
      argin = [argin {'ionfluxes'}];
    end
  end%}
  
  [tmp,ind] = ismemb(argin,possible_logvar);
  logvar = possible_logvar(sort(ind(ind>0)));
  logging = length(logvar) > 0;
  for i=1:length(logvar)
    if strcmp(logvar{i},'state_der')
      field = 'states';
    else
      field = logvar{i};
    end
    loginfo.(logvar{i}) = length(names.(field));
  end
else
  logging = 0;
  loginfo = [];
end

if ~logging
end


% Explicit set the different ion compartments
ions.types = {'Ca','Na','K'};
ionvalence = {2 1 1};

statenames.Ca = {'Cai','CaSS','CaJSR','CaNSR'};
statenames.Na = {'Nai'};
statenames.K = {'Ki'};

buffers.Ca = {'CaCMDNSS','CaCMDNi','CaCSQNJSR','LTRPNCa','HTRPNCa'};
buffers.Na = {};
buffers.K = {};

currents.Ca = {'IpCa','ICa','ICab'};
currents.Na = {'INa','INab'};
currents.K =  {'IKtof','IKtos','IK1','IKs','IKss','IKur','IKr'};

% Ion contributions from special currents (i.e. exchangers and pumps)
currents.special.Ca.INaCa = -2;
currents.special.Na.INaCa = 3;
currents.special.Na.INaK = 3;
currents.special.K.INaK = -2;

% Iterate through the ions and set the differents parameters
for i = 1:length(ions.types)
  ion = ions.types{i};
  % Setting the indexes to ion compartments
%  ions.ind.(ion).states = indfind(names.states,statenames.(ion));
  % Setting volumes to the corresponding ion compartments
  for n=1:length(statenames.(ion))
    statename = statenames.(ion){n};
    ions.conc.states.(ion).(statename) = p.(volumestr(statename));
  end
  for n=1:length(buffers.(ion))
    buffername = buffers.(ion){n};
    ions.conc.buffers.(ion).(buffername) = p.(volumestr(buffername));
  end
  % Setting ioncurrents together with their corresponding conversion
  % factor
  for n=1:length(currents.(ion))
    curname = currents.(ion){n};
    ions.currents.(ion).(curname) = p.Acap/(ionvalence{i}*p.Faraday);
  end
  special_names = fieldnames(currents.special.(ion));
  for n=1:length(special_names)
    curname = special_names{n};
    ions.currents.(ion).(curname) =currents.special.(ion).(curname)*...
        p.Acap/(ionvalence{i}*p.Faraday);
  end
end

% Index of currents that are super threshold and subthreshold
suprth = [1 5 6 9 10 14];
currents.suprth = suprth;
currents.subth = setdiff([1:1:15],suprth );

% Assigning out variables
if nxargout >= 1
 varargout{1} = loginfo;
end
if nxargout >=2
  varargout{2} = names;
end
if nxargout >= 3 
  varargout{3} = ions;
end
if nxargout >=4
  varargout{4} = currents;
end

%******************************************************
function out = volumestr(instr)
% VOLUMESTR Return a volume corresponding to the compartment name instr.
  substrings = {'SS','JSR','NSR'};

  for i=1:length(substrings)
    ind = findstr(substrings{i},instr);
    found = length(ind)>0;
    if found
      break;
    end
  end
  
  if found
    out = ['V' substrings{i}];
  else
    out = 'Vmyo';
  end

