function yint = ntrp23tb(tint,t,y,tnew,ynew,t2,y2)
%NTRP23TB  Interpolation helper function for ODE23TB.
%   YINT = NTRP23TB(TINT,T,Y,TNEW,YNEW,T2,Y2) approximates 
%   the solution at TINT by quadratic interpolation.
%   
%   See also ODE23TB, DEVAL.

%   Mark W. Reichelt, Lawrence F. Shampine, and Yanyuan Ma, 7-1-97
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2002/04/08 20:26:54 $

a1 = (((tint - tnew) .* (tint - t2)) ./ ((t - tnew) .* (t - t2)))';
a2 = (((tint - t) .* (tint - tnew)) ./ ((t2 - t) .* (t2 - tnew)))';
a3 = (((tint - t) .* (tint - t2)) ./ ((tnew - t) .* (tnew - t2)))';
yint = y*a1 + y2*a2 + ynew*a3;
