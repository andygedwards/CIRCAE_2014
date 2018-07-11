function outcell = cell2cell(varargin)
%CELL2STRCELL Convert all arguments to a one dimensional cellarray of
%             strings. Sort out all other arguments that is not string,
%             or character array.
  
  if nargin < 1
    error('To few arguments')
  end
  
  outcell = {};
  for i = 1:nargin
    if iscell(varargin{i})
      tmpcell = varargin{i};
      for n = 1:length(tmpcell)
        if iscell(tmpcell{n})
          tmp = cell2cell(tmpcell{n});
          outcell = [outcell tmp];
        else
          outcell = [outcell tmpcell(n)];
        end
      end
    elseif ischar(varargin{i})
      outcell = [outcell varargin(i)];
    end
  end
  
