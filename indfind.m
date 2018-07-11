function ind = indfind(names,varargin)
%INDFIND Index find
%   IND_FIND(CELL_ARRAY,ARG) return the indexes from the cell array, 
%   (containing strings) that is equal to the strings that is past as
%   arg. An cell array with strings can also be past as one of the arguments.
  
  if nargin < 2
    error('Not enough parameters')
  end
  
  test_names = cell2strcell(varargin);
  
  %if ~iscellstr(test_names)
  %  for i = 1:length(test_names)
  %    if ~ischar(test_names{i})
  %      test_names{i} = '\\\\\\\';
  %    end
  %  end
  %end

  if ~iscellstr(names)
    for i = 1:length(names)
      if ~ischar(names{i})
        names{i} = '///////';
      end
    end
  end
  
  [tmp,ind] = ismember(test_names,names);
  ind = sort(ind(find(ind>0)));
