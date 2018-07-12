function yint = ntrp23t(tint,t,y,tnew,ynew,h,z,znew)
%NTRP23T  Interpolation helper function for ODE23T.
%   YINT = NTRP23T(TINT,T,Y,TNEW,YNEW,H,Z,ZNEW) approximates
%   the solution at TINT by cubic Hermite interpolation.
%   
%   See also ODE23T, DEVAL.

%   Mark W. Reichelt, Lawrence F. Shampine, and Yanyuan Ma, 7-1-97
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.9 $  $Date: 2002/04/08 20:26:54 $

s = ((tint - t)/h)';
s2 = s .* s;
s3 = s .* s2;
v1 = ynew - y - z;
v2 = znew - z;
yint = y(:,ones(length(tint),1)) + z*s + (3*v1 - v2)*s2 + (v2 - 2*v1)*s3;

