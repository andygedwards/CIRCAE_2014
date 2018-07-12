function [massType, massFcn, massArgs, massM, Mvs, Mfac] = ...
    odemass(FcnHandlesUsed,ode,t0,y0,options,extras);
%ODEMASS  Helper function for the mass matrix function in ODE solvers
%    ODEMASS determines the type of the mass matrix, initializes massFcn to
%    the mass matrix function and creates a cell-array of extra input
%    arguments. ODEMASS evaluates the mass matrix at(t0,y0).  
%
%   See also ODE113, ODE15S, ODE23, ODE23S, ODE23T, ODE23TB, ODE45.

%   Jacek Kierzenka
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.8 $  $Date: 2002/04/08 20:26:55 $

massType = 0;  
massFcn = [];
massArgs = {};
massM = speye(length(y0));  
Mvs = [];
Mfac = [];

if FcnHandlesUsed     % function handles used    
  Moption = odeget(options,'Mass',[],'fast');
  if isempty(Moption)
    return    % massType = 0
  elseif isnumeric(Moption)
    massType = 1;
    massM = Moption;            
  else % try feval
    massFcn = Moption;
    massArgs = extras;  
    Mstdep = odeget(options,'MStateDependence','weak','fast');
    switch lower(Mstdep)
      case 'none'
        massType = 2;
      case 'weak'
        massType = 3;
      case 'strong'
        massType = 4;
        Mvs = odeget(options,'MvPattern',[],'fast');
      otherwise
        error(['The ODESET property ''MStateDependence'' must be set to' ...
               ' ''none'' | ''weak'' | ''strong'' ']);    
    end      
    if massType > 2   % state-dependent
      massM = feval(massFcn,t0,y0,massArgs{:});
    else   % time-dependent only
      massM = feval(massFcn,t0,massArgs{:});
    end
  end  
  
else % ode-file
  mass = lower(odeget(options,'Mass','none','fast'));
  % Obsolete code -- should be removed in a subsequent release.
  Mconstant = odeget(options,'MassConstant',[],'fast');
  if strcmp(mass,'on') | strcmp(mass,'off') | ...
      strcmp(Mconstant,'on') | strcmp(Mconstant,'off')
    if strcmp(mass,'on') | strcmp(mass,'off')
      warning('MATLAB:odemass:ObsoleteMassValue',['Mass property values ' ...
              'are ''none'', ''M'', ''M(t)'', and ''M(t,y)'' (see ODESET).' ...
              '\n         The ''on'' and ''off'' values are obsolete ' ...
              'and may not be supported in a\n         future release.']);
    end
    if strcmp(Mconstant,'on') | strcmp(Mconstant,'off')
      warning('MATLAB:odemass:ObsoleteMassConstantProperty', ...
            ['The MassConstant property is obsolete and will disappear ' ...
            '\n         in a future release (see the Mass property of ' ...
            'ODESET).']);
    end
    if strcmp(Mconstant,'on')
      mass = 'm';
      warning('MATLAB:odemass:AssumeConstant', ...
          'Assuming Mass value ''M'', a constant mass matrix.');
    elseif strcmp(mass,'on')
      mass = 'm(t)';
      warning('MATLAB:odemass:AssumeTimeDependent', ...
          'Assuming Mass value ''M(t)'', a time-dependent mass matrix.');
    else
      mass = 'none';
      warning('MATLAB:odemass:AssumeNone', ...
          'Assuming Mass value ''none'', no mass matrix.');
    end
  end

  switch(mass)
    case 'none', return;  % massType = 0
    case 'm', massType = 1;
    case 'm(t)', massType = 2;
    case 'm(t,y)', massType = 3;
    otherwise
      error(['Unrecognized property ''Mass'': ' mass ]);
  end
  massFcn = ode;  
  massArgs = [{'mass'}, extras];
  massM = feval(massFcn,t0,y0,massArgs{:});    

end



