function yinterp = ntrp23s(tinterp,t,y,tnew,ynew,h,k1,k2)
%NTRP23S  Interpolation helper function for ODE23S.
%   YINTERP = NTRP23S(TINTERP,T,Y,TNEW,YNEW,H,K1,K2) uses data computed in
%   ODE23S to approximate the solution at time TINTERP.
%   
%   See also ODE23S, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.12 $  $Date: 2002/04/08 20:26:54 $

s = ((tinterp - t) / h)';       % may be a row vector

d = 1 / (2 + sqrt(2));
e = h / (1 - 2*d);

p1 = (s .* (1 - s)) * e;
p2 = (s .* (s - 2*d)) * e;

yinterp = y(:,ones(length(tinterp),1)) + k1*p1 + k2*p2;
