function yinterp = ntrp45(tinterp,t,y,tnew,ynew,h,f)
%NTRP45  Interpolation helper function for ODE45.
%   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F) uses data computed in ODE45
%   to approximate the solution at time TINTERP.
%   
%   See also ODE45, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.13 $  $Date: 2002/04/08 20:26:55 $

BI = [
    1       -183/64     37/12       -145/128
    0       0       0       0
    0       1500/371    -1000/159   1000/371
    0       -125/32     125/12      -375/64 
    0       9477/3392   -729/106    25515/6784
    0       -11/7       11/3        -55/28
    0       3/2     -4      5/2
    ];

s = ((tinterp - t) / h)';       % may be a row vector

yinterp = y(:,ones(length(tinterp),1)) + f*(h*BI)*cumprod(s(ones(4,1),:));
