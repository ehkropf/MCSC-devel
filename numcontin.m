% Numercial continuation, an MCSC solver
%   Based on Algower & Georg, algorithm 3, appropriate credit, etc.
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

% Copyright Everett Kropf, 2013,
% except as may be noted below.

classdef numcontin < nesolver

  properties
    Ofun; % objective function object
    X0;
    f0;
    
    lambda = 1;
    hlast = inf;
    soltime = inf;
    ofunnorm = 1;
    tolerance = 1e-15;
    
    hmin = 1e-14; % min stepsize
    hmax = 1; % max stepsize
    resmin = 1e-8;
    contrmax = 0.5;
    perturb = 0; %1e-12; %sqrt(eps);
    anglemax = pi/4; % max tangent turning angle
  end

  methods
    function S = numcontin(O,X0)
      if nargin > 0
        
        S.X0 = X0;
        S.f0 = O(X0);
        S.Ofun = O;
      end
    end % constructor
    
    function run_solver(S)
      tic
      continuation(S);
      secs = toc;
      S.soltime = secs;
      
%       if abs(1 - S.lambda) > 1e-6
%         warning('MCSC:checkoutput',...
%           'Solver output may not be trustworthy.')
%       end
      
      fprintf('\nContinuation took %.3f seconds.\n',secs)
      fprintf('Final 1-lambda value was %.6g.\n',1-S.lambda)
      fprintf('Inf norm of objective function is %.6g\n\n',...
        S.ofunnorm); %norm(S.Ofun(S.lastx(2:end)),inf))
    end % run_solver
    
    function continuation(S)
      % continuation(curvefun,x,refdir,monitor,target)
      % Euler-Newton method of continuation for equations F(x)=0, with Broyden
      % update of a numerical jacobian.
      %
      % F maps from R^{n+1} to R^n, describing a curve to follow.
      % Euler predictor, chord-Newton corrector with stepsize control.
      % See Program 3 in Allgower & Georg.
      %
      % Input or set initial point, initial step, min step, corrector tol, max
      % contract, max distance from curve.
      
      h = 1/16;
      monitor = 1;
      target = 1;
      x = [0;S.X0];
      refdir = [1;0*S.X0];
      N = length(x)-1;
      Fx = curvefun(S,x);

      fprintf('Computing initial FD Jacobian...\n')
      A = numjac(S,x,Fx);
      
      % Find a tangent direction (null(A)).
      t = S.tangent(A,1);
      orient = sign(t'*refdir);
      t = orient*t;
      
      adaptnewton = false;
      refinecount = 0;
      
      fprintf('\n  Stepsize     1-lambda      ||F(x)||\n')
      fprintf(  ' ---------    ---------    ----------\n')
      
      while (abs(h) > S.hmin) && (S.ofunnorm > S.tolerance)
        % After several consecutive refinements, re-initialize the Jacobian.
        if refinecount > 3
          fprintf('Recomputing FD Jacobian...\n')
          fprintf('\n  Stepsize     1-lambda      ||F(x)||\n')
          fprintf(  ' ---------    ---------    ----------\n')

          A = numjac(S,x,Fx);
          t = S.tangent(A,orient);
          refinecount = 0;
        end
        
        u = x + h*t;
        Fu = curvefun(S,u);
        
        % Predictor update.
        s = t;
        A = A + (Fu-Fx)*t'/h;                 % Broyden
        [t,Q,R] = S.tangent(A,orient);          % tangent
        if acos(s'*t) > S.anglemax
          % Change in tangent was too drastic.
          h = h/2;
          refinecount = refinecount+1;
          continue
        end
        
        % Perturb to avoid cancellation. (?)
        pf = -S.perturb*sign(Fu) .* double(abs(Fu) < S.perturb);
        
        % Corrector step.
        du = Q(:,1:N) * (R(1:N,1:N)' \ (Fu-pf));
        v = u - du;
        Fv = curvefun(S,v);
        
        % Corrector update.
        A = A + (Fv-pf)*(-du)' / norm(du)^2;
        t = S.tangent(A,orient);            % tangent
        
        % Contraction test.
        contract = norm(Fv-pf) / (norm(Fu-pf) + S.resmin);
        if contract > S.contrmax
          h = h/2;
          refinecount = refinecount+1;
        else
          S.ofunnorm = norm(Fv-(v(1)-1)*S.f0,inf);
%           fprintf(['Stepsize = % 9.3g;  lambda =% 11.5g;  ',...
%             '||F(x)|| =% 10.4g\n'],...
%             h,v(monitor),S.ofunnorm);
          fprintf(' % 9.3g    % 9.3g    % 10.4g\n',...
            h,1-v(monitor),S.ofunnorm);

          x = v;
          Fx = Fv;
          refinecount = 0;
          
          % Switch to newton for the endgame?
          if v(monitor) > target
            if ~adaptnewton
              fprintf('Entering final Newton phase.\n');
              fprintf('\n  Stepsize     1-lambda      ||F(x)||\n')
              fprintf(  ' ---------    ---------    ----------\n')
            end
            adaptnewton = true;
          end
          
          % Adapt step size.
          if adaptnewton
            h = -(v(monitor)-target) / t(monitor);
          else
            h = min( S.hmax, 2*h );
          end
        end
        
      end

      S.hlast = h;
    end % continuation
    
    function H = curvefun(S,X)
      % the homotopy function
      S.lambda = X(1);
      H = S.Ofun(X(2:end)) + (X(1)-1)*S.f0;
    end % curvefun
    
    function J = numjac(S,x,Fx)
      N = length(x)-1;
      J = zeros(N,N+1);
      for k = 1:N+1
        xk = x;  xk(k) = xk(k)+1e-6;
        J(:,k) = 1e6 * (curvefun(S,xk) - Fx);
      end
    end
    
  end % methods

  methods (Static)    
    function [t Q R] = tangent(A,orient)
      N = size(A,1);
      [Q,R] = qr(A');
      t = Q(:,N+1);
      t = orient * sign(det(Q)) * prod(sign(diag(R))) * t;
    end
    
  end % static methods
  
end
