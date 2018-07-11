function [out, varargout]= ismemb(A,S)
% ISMEMB True for set member.
%    The same function as the built in ismember, but can handle cells
%    that are not str, i.e. converts all none char in S to a the dummy
%    string '///////'.
  
% Remove all none char ellement
  if ~iscellstr(S)
    for i = 1:length(S)
      if ~ischar(S{i})
        S{i} = '///////';
      end
    end
  end

  if ~iscellstr(A) && ~ischar(A)
    for i = 1:length(A)
      if ~ischar(A{i})
        A{i} = '\\\\\\\';
      end
    end
  end
  
  [out ind]= ismember(A,S);
  if nargout > 1
    varargout{1} = ind;
  end
