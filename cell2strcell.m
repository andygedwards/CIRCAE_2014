function strcell = cell2strcell(varargin)
%CELL2STRCELL Convert all arguments to a one dimensional cellarray of
%             strings. Sort out all other arguments that is not string,
%             or character array.
  
  if nargin < 1
    error('To few arguments')
  end
  
  strcell = {};
  for i = 1:nargin
    if iscell(varargin{i})
      tmpcell = varargin{i};
      for n = 1:length(tmpcell)
        if iscell(tmpcell{n})
          tmp = cell2strcell(tmpcell{n});
          strcell = [strcell tmp];
        elseif ischar(tmpcell{n})
          strcell = [strcell tmpcell(n)];
        end
      end
    elseif ischar(varargin{i})
      strcell = [strcell varargin(i)];
    end
  end
  
