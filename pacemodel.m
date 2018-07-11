function [Tout,Sout,varargout] = pacemodel(cellmodel,paceintervall,numpace,varargin)
%PACEMODEL  Pace a given cell model
%   [T,S] = PACEMODEL(CELLMODEL,PACEINTERVALL,NUMPACE) paces a the given
%   CELLMODEL NUMPACE times with intervall PACEINTERVALL. The states and
%   time vectors are returned in T, and S as NUMPACE numbers of matrixes
%   stored in cellarrays. CELLMODEL is a function handler to a function
%   that delivers the derivative of the states.
%
%   Possible parameters are:
%   'prepace'   - The model are prepaced with no stimuli for
%                 paceintervall number of ms.
%
%   misc log    - Variables that are not state varibles can be logged by 
%   parameters    defining them in the corresponding Init file, and to
%                 actually log them in the cellmodel function.
%
%   'ioninfo'   - Calculate ionic informations during a paceing. This
%                 concist of specific ioncurrents, concentrations of ions
%                 and buffers and the actuall numbers of ions. This
%                 option assume that both 'currents' and 'ionbuffers' are
%                 valid logparameters.
%  
%   model       - Model parameters and initial values that are set in the
%   parameters    init file to the cellmodel can be set. This is done by
%                 passing a struct that contains the name of the parameter 
%   initial       or state variable as fieldname.
%   statevalues   If the name is not present as a parameter or state 
%                 in the init file, it will be ignored. You will get
%                 a confirmation of wich parameters or initial states
%                 that are successfully set by you when simulation is
%                 started.
%  
%   dynamic     - If 'dynamicstop' is passed then in addition to reaching 
%   stop          numpace number of paces then the simulation can be
%                 stopped if the states do not change more than 1 %.
%                 If an array  is passed after 'dynamicstop', then only 
%                 these states are checked for dynamicstop. A second
%                 numerical value can be passed after the array, setting
%                 the stop criteria to something else than 1 %.
%  
%  EXAMPLE:      
%  >> ind_V = 33
%  >> [T,S] = pacemodel(@Bond,100,3)
%  >> plot(T{1},S{1}(:,ind_V),T{2},S{2}(:,ind_V),T{3},S{3}(:,ind_V))
%  >> param.v3 = 0.55
%  >> x0.V     = -80
%  >> [T,S] = pacemodel(@Bond,100,3,param,x0)
%  >> [T,S,loginfo,x1,ioninfo] = pacemodel(@Bond,100,3,'ioninfo')
%    x
%    
  nxargin = max(nargin-3,0);
  nxargout = max(nargout-2,0);
  if nargin < 3
    error('Not enough input parameters')
  end
  if nargout < 2
    error('Too few output arguments')
  end
  if numpace < 0
    error('Number of paces have to be positive.')
  end
  if paceintervall < 0
    error('Pace intervall have to be positive.')
  end
  
  inputstruct = {};
  argin = {};
  prepace = 0;
  dynamicstop = 0;
  
  % If there are more than three inputparameters
  if nxargin > 0
    % Flatten the inputparameters to a one dimensional cellarray
    argin = cell2cell(varargin);
  end

  % Call the init function
  cellmodel_init = str2func([func2str(cellmodel) 'Init']);
  [p, structx0, loginfo, names, ioninfo, currents] = ...
      feval(cellmodel_init,argin);

  % Collects any structs in the input, i.e. if the user wants to set
  % any initial parameters or values
  for i = 1:length(argin) 
    if isstruct(argin{i})
      fields = fieldnames(argin{i})';
      % Copy the fields in the struct to inputstruct
      for field =fields
        inputstruct.(field{1}) = argin{i}.(field{1});
      end
    end
  end
  
  % If ioninfo are to be calculated
  if ismemb('ioninfo',argin)
    makeioninfo = 1;
    if nxargout < 3
      error('If logging ioninfo then 5 or more output variables are needed.')
    end
  else
    makeioninfo = 0;
  end
 
  % Check for 'prepace' arguments
  if ismemb('prepace',argin)
    prepace = 1;
    ind = indfind(argin,'prepace');
    if length(argin)>ind && isnumeric(argin{ind+1})
      prepace_time = argin{ind+1};
      if prepace_time < 0
        error('Prepace time have to be positive!')
      end
    else
      prepace_time = paceintervall;
    end
  end

  % Check for 'dynamicstop' arguments
  if ismemb('dynamicstop',argin)
    dynamicstop = 1;
    ind = indfind(argin,'dynamicstop');
    dynamicstop_states = [1:1:length(names.states)];
    dynamicstop_value = 1e-2; % A stop value of 1% as default
    if length(argin)>ind && isnumeric(argin{ind+1})
      dynamicstop_states = argin{ind+1};
        if length(argin)>ind +1 && isnumeric(argin{ind+2})
          dynamicstop_value = argin{ind+2};
        end
    end
  end
  
  % If any values is collected in the inputstruct then change these for
  % the correspondning fields in x0 or p.
  pchange = {};
  x0change = {};
  if length(inputstruct)>0
    inputfields = fieldnames(inputstruct)';
    for field=inputfields
      field = field{1};
      if isfield(p,field)
        p.(field) = inputstruct.(field);
        pchange.(field) = inputstruct.(field);
      elseif isfield(structx0,field)
        structx0.(field) = inputstruct.(field);
        x0change.(field) = inputstruct.(field);
      end
    end
  end
  x0 = struct2array(structx0);

  disp('Pacemodel')
  sum(x0(1:12))+x0(43)
  
  % Initilize the output arguments
  Tout = cell(1,prepace + numpace);
  Sout = cell(1,prepace + numpace);
  
  % If logging
  if length(loginfo)>0;
    logging = 1;
    logvar = fieldnames(loginfo);
    nlogvar = length(logvar);
    for n=1:nlogvar
      logout.(logvar{n}) = cell(1,prepace + numpace);
    end
  else
    logging = 0;
  end
  disp('    **********************************************************')
  disp(['    Paceing ' names.modelname ' ' num2str(numpace) ' times,'])
  disp(['    with intervall of ' num2str(paceintervall) ' ms, '...
       num2str(round(1000/paceintervall*60)) ' beats/min.'])

  if (prepace)
    disp(' ')
    disp(['    Prepace with zero stimuli for ' num2str(prepace_time) ' ms.'])
  end
  
  if (logging)
    out = '    Logging for: ';
    for log = logvar'
      out = [out log{1} ', '];
    end
    disp(' ')
    disp([out(1:end-2) '.'])
    disp('    The result is returned in third output variable.')
  end
  
  %makeioninfo=1; % hack
  if makeioninfo
    disp(' ')
    disp('    Ioninfo are calculated. Returned in fifth output variable')
  end
  
  if dynamicstop
    disp(' ')
    disp(['    Simulation can be stopped by dynamic stop, with stop criteria '...
          num2str(dynamicstop_value*100) ' %.'] )
    if ~(length(dynamicstop_states) == length(names.states))
      disp('    The following states are checked:')
      disp(names.states(dynamicstop_states))
    end
    disp('')
   
  end
  if length(pchange)>0
    disp(' ')
    disp('    Following parameters are set by user:')
    disp(pchange)
  end
  
  if length(x0change)>0
    disp(' ')
    disp('    Following initial values are set by user:')
    disp(x0change)
  end
  disp('    **********************************************************')
  disp(' ')
  
  nout = 1;
  modeltime = 0.0;

  % Take the realtime of the whole simulation
  % Do the prepace, if any.
  if prepace
%    disp(' ')
%    disp(['*** Pre-pace in modeltime: ' num2str(prepace_time) ' ms. ****'])
    tmp = p.pacedur;
    p.pacedur = 0.0;
    [T,S,logres] = paceone(cellmodel,prepace_time,x0,p,loginfo);
    p.pacedur = tmp;
    x0 = S(end,:);
    Tout{nout} = T;
    Sout{nout} = S;
    if logging
      for n=1:nlogvar
        logout.(logvar{n}){nout} = logres.(logvar{n});
      end
    end
    nout = nout+1;
    modeltime = prepace_time;
  end
  
  acc_time = 0;
  stoped_by_dynamicstop = 0;
  % The main paceing
  for i = 1:numpace
    tic
      modeltime = modeltime + paceintervall;
      if ~dynamicstop
        disp(['    Pace # ' num2str(i) ' out of ' num2str(numpace) '. In ' ...
                          'modeltime: ' num2str(modeltime) ' ms.'])
      end
      [T,S,logres] = paceone(cellmodel,paceintervall,x0,p,loginfo);
      Tout{nout} = T;
      Sout{nout} = S;
      x0 = S(end,:);
      if logging
        for n=1:nlogvar
          logout.(logvar{n}){nout} = logres.(logvar{n});
        end
      end
      acc_time = acc_time+toc;
      if dynamicstop && i > 1 && i ~= numpace
        [diff ind] = statecompare(Sout{i-1},S,dynamicstop_states);
        
        if diff < dynamicstop_value 
          stoped_by_dynamicstop = 1;
          Tout = Tout(1:nout);
          Sout = Sout(1:nout);
          if logging
            for n=1:nlogvar
              logout.(logvar{n}) = logout.(logvar{n})(1:nout);
            end
          end
          disp(' ')
          disp('    Stopped by dynamic stop criteria.')
          disp(['    State ' names.states{ind} ...
                ' finaly did not change more than ' ...
                num2str(dynamicstop_value*100) ' %.'])
          break;
        end
        disp(['    Pace # ' num2str(i) ' of maximal ' num2str(numpace)]) 
        disp(['    ' names.states{ind} ' changed by ' num2str(ceil(diff*1000)/10) ' %.'])
        disp(' ')
      end
      nout = nout+1;
      time_estimate = max(acc_time/i*(numpace-i),0.0);
      if ~dynamicstop
        disp(' ')
        disp(['    Estimated time left: ' formattime(time_estimate) '.'])
      end
    end
    
  % Print the realtime of the whole simulation
  
  disp(' ')
  disp(['    Total simulation time ' formattime(acc_time) '.'])
  disp(' ')
  
  if nxargout >= 1
    logout.names = names;
    varargout{1} = logout;
  end

  % If ioninfo are to be calculated
  ioninfoout = [];
  
  if makeioninfo
    ioninfoout = ionstat(Tout,Sout,logout,ioninfo);
  end
  
  if nxargout >=2
    varargout{2} = initialcond(Sout,names.states);
  end

  if nxargout >=3
    varargout{3} = ioninfoout;
  end
  
  if nxargout >=4
    varargout{4} = currents;
  end
  
%************************************************************************
function strout = formattime(secin)
  if secin < 0
    error('Only positive arguments!')
  end
  strout = '';
  minutes = floor(secin/60);
  if minutes > 0
    strout = [ num2str(minutes) ' min '];
  end
  seconds = round(mod(secin,minutes*60));
  strout = [strout num2str(seconds) 's'];
