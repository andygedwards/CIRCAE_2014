function s = funstring(fun)
% Yield a string representing fun.

%   Mike Karr, Jacek Kierzenka, 11-19-99
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2002/04/08 20:26:54 $

if isa(fun, 'function_handle')
  s = upper(func2str(fun));
elseif isstr(fun)
  s = upper(fun);
elseif isa(fun, 'inline')
  s = formula(fun);
else
  s = 'unknown';
end
