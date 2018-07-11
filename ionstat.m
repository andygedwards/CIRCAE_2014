function infoout = ionstat(T,S,logres,ioninfo)
%IONSTAT Return info about all ion handling from a model
%   INFOOUT = IONSTAT(T,S,LOGRES,IONINFO)
%
%   INFOOUT is on the following form
%    COMPARTMENTS: [1x1 struct]
%           SINKS: [1x1 struct]
%          SYSTEM: [1x1 struct]
%  
%   COMPARTMENTS holds the different ionic concentrations and number
%   of ions in each compartments in the model for each iontype.
%  
%   SINKS holds the ionderivatives for each ion current and the
%   accumulated sum of all currents of a specific ion current.
%  
%   SYSTEM holds the total number of ions in the cell and the fractional
%   error.

  allstates = logres.names.states;
  allbuffers = logres.names.buffers;
  allcurrents = logres.names.currents;
  buffers = logres.buffers;
  currents = logres.currents;
  ncell = length(S);
  
  for ion = ioninfo.types
    ion = ion{1};
    currentind = [];
    bufferind = [];
    stateind = [];
    % ********** Collecting info about the input structure *************
    % Index and names of the states that includes any ion concentration
    % info, excluding buffers.
    if isfield(ioninfo.conc.states,ion)
      statenames = fieldnames(ioninfo.conc.states.(ion));
      stateind = indfind(allstates,statenames);
    end
    % Index and names of the buffers that contain ions in the model
    if isfield(ioninfo.conc.buffers,ion)
      buffernames = fieldnames(ioninfo.conc.buffers.(ion));
      bufferind = indfind(allbuffers,buffernames);
    end
    % Index and names of the membrane currents that are of specific ions
    if isfield(ioninfo.currents,ion)
      currentnames = fieldnames(ioninfo.currents.(ion));
      currentind = indfind(allcurrents,currentnames);
    end
    % Index for the collected ion currents going out of the cell
    ionoutind = indfind(allstates,[ion 'out']);
    
    % Assigning variables for output
    compartments.names =  [allstates(stateind) allbuffers(bufferind)];
    compartments.conc = cell(1,ncell);
    compartments.numions = cell(1,ncell);
    compartments.sumions = cell(1,ncell);
    
    sinks.names = currentnames;
    sinks.ionsder = cell(1,ncell);
    sinks.cum = cell(1,ncell);
  
    system.sumions = cell(1,ncell);
    system.fracerror = cell(1,ncell);
    
    % Collecting information for each pace in the simulation
    for n=1:ncell
      % The number of time steps in the present pace
      numt = size(S{n},1);
      % Assigning space for the current sinks
      sinks.ionsder{n} = zeros(numt,length(currentind));
      % Collecting currents sinks
      for j=1:length(currentind)
         sinks.ionsder{n}(:,j)= currents{n}(:,currentind(j))*...
             ioninfo.currents.(ion).(currentnames{j});
      end
      % Collecting the cumulated ions going out of the cell 
      sinks.cum{n} = S{n}(:,ionoutind);
      
      % Collecting the actuall concentrations of the ion compartments
      compartments.conc{n} = [S{n}(:,stateind) buffers{n}(:,bufferind) ];
      % Assigning space for the ion numbers in the corrersponding compartments
      compartments.numions{n} = zeros(numt,length(bufferind)+length(stateind));
      % Collecting the ion numbers in the corrersponding compartments...
      for j=1:length(stateind)
        compartments.numions{n}(:,j) = S{n}(:,stateind(j))*...
            ioninfo.conc.states.(ion).(statenames{j});
      end
      % ... and buffer
      for x=1:length(bufferind)
        compartments.numions{n}(:,x+j) =  buffers{n}(:,bufferind(x))*...
            ioninfo.conc.buffers.(ion).(buffernames{x});
      end
      % Sum all ions in the different compartments
      compartments.sumions{n} = sum(compartments.numions{n},2);
      % Sum all ions in compartments and the accumulated ions going out
      % of the cell
      system.sumions{n} = sinks.cum{n} + compartments.sumions{n};
      % Compute the fractional error, compared with the total amount of
      % ions in the first time step of the first pace
      system.fracerror{n} = system.sumions{n}/system.sumions{1}(1)-1;
    end
  % Assigning out variables
  infoout.compartments.(ion) = compartments;
  infoout.sinks.(ion) = sinks;
  infoout.system.(ion) = system;
  end
  
