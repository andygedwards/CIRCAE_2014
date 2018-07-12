function out = odemxv(t,y,Mfun,Margs,v)
%ODEMXV  Helper function -- evaluates M(t,y)*v
%   Used to get d(M(t,y)*v)/dy when the property MStateDependence is 'strong'  
%
%   See also IC3DAE, ODE15S, ODE23T, ODE23TB.

%   Jacek Kierzenka, Lawrence Shampine
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.5 $  $Date: 2002/04/08 20:26:55 $

out = feval(Mfun,t,y,Margs{:})*v;

