% RSTOBJFUN provides objective function for (rectangle) resistor problems.

% E. Kropf, 2014

classdef rstobjfun < objectivefun

properties
  Prect_
  Cr_
  gencirc_
  N_
  quad_
end

methods
  function F = rstobjfun(Prect, Cr, gencirc, opts)
    if nargin
      args{1} = Prect;
      
      if nargin < 4 && isempty(opts)
        opts = intmapopts;
      end
      args{2:3} = {opts.monitor, opts.fignum};
    else
      args = {};
    end
    F = F@objectivefun(args{:});
    
    if nargin
      % Save data to self.
      F.Prect_ = Prect;
      F.Cr_ = Cr;
      F.gencirc_ = gencirc;
      F.N_ = opts.N;
      
      % Create quadrature object.
      fprime = opts.fpclassh(Prect, Cr, opts.N);
      F.quad_ = gjquad(Prect, fprime, opts.ngj, opts.halfrule);
    end
  end
  
  function Xc = constrain(F, Xu)
    % Convert unconstrained circle unknowns back to human levels.
    
    m = F.Cr_.m;
    Xc = zeros(2*(m-1), 1);
    for j = 1:(m-1)
      Xu(2*j-1) = mod(Xc(2*j-1), 2*pi);
      Xu(2*j) = Xc(2*j-1) + 2*pi/(1 + exp(Xc(2*j)));
    end
    
    if any(isinf(Xu))
      warning('mcsc:unconstr', ...
        'At least one of the unconstrained variables is at infinity.')
    end
  end
  
  function f = fun_eval(F, Xu)
    % Calculate objective function.
    C = F.Cr_;
    c = C.c;
    r = C.r;
    tkj = C.t;
    [m, ~, vl] = polydat(F.Prect_);
    quad = F.quad_;
    
    % Convert unconstrained and store in C.
    Xc = unconstrain(F, Xu);
    for j = 1:(m-1)
      C.t(j+1,1:2) = Xc(2*j-1:2*j);
    end
    
    % Compute integration constant.
    A = arcq(quad, tkj(1,1), 1, tkj(2,1), 2, 0, 1);
    A = (vl(2,1) - vl(1,1))/A;
    
    % Integrate along both sides of inner slits.
    SL = nan(2, m-1);
    for j = 1:m-1
      off = 4 + 2*j;
      lpvn = off + (1:2)';
      left = tkj(1:2,j+1);
      rpvn = flipud(lpvn);
      right = tkj(2:1,j+1);
      right(left > right) = right(left > right) + 2*pi;
      cv(1:2) = c(j+1);
      rv(1:2) = r(j+1);
      
      SL(1:2,j) = A*arcq(quad, left,lpvn, right,rpvn, cv,rv);    
    end
    
    % Components of objective function.
    f = zeros(2*(m-1), 1);
    for j = 1:m-1
      f(2*j-1) = sin(angle(SL(1,j+1)));
      f(2*j) = abs(SL(1,j+1)) - abs(SL(2,j+1));
    end
  end
  
  function plot(F)
    % Draw thyself!
  end
  
  function out = subsref(F, S)
    % Provide behaviour:
    %   f = F(X) -- as function,
    %   revert to default otherwise.    
    if numel(S) == 1 && strcmp(S.type, '()')
      out = fun_eval(F, S.subs{:});
    else
      out = builtin('subsref', F, S);
    end
  end
  
  function Xu = unconstrain(F, C)
    % Give unconstrained list of unknown C variables (the interior slit tip
    % pre-images).
    
    if nargin < 2
      C = F.Cr_;
    end
    m = C.m;
    t = C.t;
    Xu = zeros(2*(m-1), 1);
    for j = 1:(m-1)
      phi = diff([t(j+1,1:2); t(j+1,1) + 2*pi]);
      phi = mod(phi, 2*pi); % enforce sum(phi) = 2*pi
      Xu(2*j-1) = t(j+1,1);
      Xu(2*j) = log(phi(2)/phi(1));
    end
  end
end

end
