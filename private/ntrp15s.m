function yinterp = ntrp15s(tinterp,t,y,tnew,ynew,h,dif,k)
%NTRP15S  Interpolation helper function for ODE15S.
%   YINTERP = NTRP15S(TINTERP,T,Y,TNEW,YNEW,H,DIF,K) uses data computed in
%   ODE15S to approximate the solution at time TINTERP.
%   
%   See also ODE15S, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.12 $  $Date: 2002/04/08 20:26:54 $

s = ((tinterp - tnew) / h)';        % may be a row vector

if k == 1
  yinterp = ynew(:,ones(length(tinterp),1)) + dif(:,1) * s;
else                    % cumprod collapses vectors
  K = (1:k)';
  kI = K(:,ones(length(tinterp),1));
  yinterp = ynew(:,ones(length(tinterp),1)) + ...
      dif(:,K) * cumprod((s(ones(k,1),:)+kI-1)./kI);
end
