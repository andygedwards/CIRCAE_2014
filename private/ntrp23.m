function yinterp = ntrp23(tinterp,t,y,tnew,ynew,h,f)
%NTRP23  Interpolation helper function for ODE23.
%   YINTERP = NTRP23(TINTERP,T,Y,TNEW,YNEW,H,F) uses data computed in ODE23
%   to approximate the solution at time TINTERP.
%   
%   See also ODE23, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.12 $  $Date: 2002/04/08 20:26:54 $

BI = [
    1       -4/3        5/9
    0       1       -2/3
    0       4/3     -8/9
    0       -1      1
    ];

s = ((tinterp - t) / h)';       % may be a row vector

yinterp = y(:,ones(length(tinterp),1)) + f*(h*BI)*cumprod(s(ones(3,1),:));
