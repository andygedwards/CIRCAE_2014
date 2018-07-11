function [T,S,varargout] = paceone(cellmodel,pacetime,varargin)
% PACEONE Pace a givven cellmodel ones
%    
%  x0,p,loginfo)
  nxargout = max(nargout-2,0);
  nxargin = max(nargin-2,0);
  if nargin < 2
    error('Too few input arguments.')
  end
  if nargout < 2
    error('Too few output arguments.')
  end
  
  % The model needs to be initialized
  if nxargin <= 1
    if nxargin == 1
      arg = varargin{1};
    else
      arg = {};
    end
    % Call the init function
    cellmodel_init = str2func([func2str(cellmodel) 'Init']);
    [p, x0,loginfo] = feval(cellmodel_init,arg);
    x0 = struct2array(x0)';
  end
  
  
  % The initparameters are supplied
  if nxargin == 2
    x0 = varargin{1};
    p = varargin{2};
    loginfo = [];
  end
  
  if nxargin == 3
    x0 = varargin{1};
    p = varargin{2};
    loginfo = varargin{3};
  end
  options = [];
  if length(p.VClampAmp)>0
    options = odeset('RelTol',5e-6,'AbsTol',5e-8);
  end
%  options = odeset('RelTol',1e-12,'AbsTol',1e-15);%, ...
%                   'OutputFcn','odeplot','OutputSel',33);
%  options = odeset('MaxStep',0.05);%,'OutputFcn','odeplot','OutputSel',33);
  
  % Use different ode solver depending on if logging or not
  if length(loginfo)>0
    if p.pacestart == 0 && length(p.VClampAmp)==0
      [T,S,logres] = ode15slog(cellmodel,[0 pacetime],x0,options,p,loginfo);
      varargout{1} = logres;
    elseif ~(p.pacestart == 0)
      [T,S,logres] = ode15slog(cellmodel,[0 p.pacestart],x0,options,p, ...
                               loginfo);
      [tmpT,tmpS,tmplogres] = ode15slog(cellmodel,[p.pacestart pacetime],S(end,:),...
                                        options,p,loginfo);
      T = [T;tmpT(2:end)];
      S = [S;tmpS(2:end,:)];
      logres.currents = [logres.currents;tmplogres.currents(2:end,:)];
      logres.buffers = [logres.buffers;tmplogres.buffers(2:end,:)];
      varargout{1} = [logres ;tmplogres(2:end,:)];
    else
      timevec = [0 p.VClampTimes];
      T = [];
      S = [];
      tmpS = x0;
      for i = 1:length(timevec)-1
        [tmpT,tmpS,tmplogres] = ode15slog(cellmodel,[timevec(i) timevec(i+1)],...
                                          tmpS(end,:),options,p,loginfo);
        T = [T;tmpT(2:end)];
        S = [S;tmpS(2:end,:)];
        if i == 1
            logres.currents = tmplogres.currents(2:end,:);
            logres.ionfluxes = tmplogres.ionfluxes(2:end,:);
            logres.buffers = tmplogres.buffers(2:end,:);   
        else
            logres.currents = [logres.currents;tmplogres.currents(2:end,:)];
            logres.ionfluxes = [logres.ionfluxes;tmplogres.ionfluxes(2:end,:)];
            logres.buffers = [logres.buffers;tmplogres.buffers(2:end,:)];
        end
      end
      varargout{1} = logres;
    end
  else
    if p.pacestart == 0 && length(p.VClampAmp)==0
      [T,S] = ode15s(cellmodel,[0 pacetime],x0,options,p);
    elseif ~(p.pacestart == 0)
      [T,S] = ode15s(cellmodel,[0 p.pacestart],x0,options,p);
      [tmpT,tmpS] = ode15s(cellmodel,[p.pacestart pacetime],S(end,:),...
                           options,p);
      T = [T;tmpT(2:end)];
      S = [S;tmpS(2:end,:)];
    else
      timevec = [0 p.VClampTimes];
      T = [];
      S = [];
      tmpS = x0;
      for i = 1:length(timevec)-1
        [tmpT,tmpS] = ode15s(cellmodel,[timevec(i) timevec(i+1)],tmpS(end,:),...
                             options,p);
        T = [T;tmpT(2:end)];
        S = [S;tmpS(2:end,:)];
      end
    end
    for i=1:nxargout
      varargout{i} = [];
    end
  end
  
