function [f,finv,a]=uself(c,r)
% USELF a self-map of the unit disk 
%    F=USELF(C,R) is based on F(Z)=LAM*(Z-A)/(1-Z*CONJ(A)) with
%    |A|<1, F(A)=0, and |LAM|=1. Here A is calculated so that the circle
%    with center C and radius R is mapped to be the inner circle of an
%    annulus.
%
%
% This file is part of MCSC.
% 
% MCSC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MCSC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MCSC.  If not, see <http://www.gnu.org/licenses/>.

% Copyright Everett Kropf, 2008, 2013,
% except as may be noted below.

% put c on the positive real axis
k = exp(1i*angle(c));
c = abs(c);

% figure out a
b = 1+(c^2-r^2);
discr = sqrt(c^4+r^4-2*(c^2*r^2+c^2+r^2)+1);
if b >= discr
    a = (b-discr)/(2*c);
else
    a = (b+discr)/(2*c);
end

% rerotate the circle
a = a*k;

% viola
f = @(z) (z-a)./(1-conj(a)*z);

if nargout>1
    finv = @(w) (w+a)./(1+conj(a)*w);
end
