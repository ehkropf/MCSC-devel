% RSTMAP creates the rectangular resistor map.

% E. Kropf, 2014

classdef rstmap < mcscmap

properties
  
end

methods
  function F = rstmap(C, genrect, opts)
    % Here we need to find a rectangle with horizontal slits, and the slit tip
    % preimages on a given circle domain.
    % INPUTS:
    %   C = Circle domain, which in general is a solution to a parameter problem
    %   for another polygonal domain.
    %   genrect = List of vertices on C_1 in the circle domain that define the
    %   generalized rectangle.
    
    if nargin == 0
      args = {};
    else
      if ~isa(C, 'circdomain')
        error('First argument must be a circle domain.')
      end
      if ~isequal(size(genrect(:)), [4, 1])
        error('Second argument is a 4 element vector.')
      end
      if nargin < 3 || isempty(opts)
        opts = intmapopts;
      end
      
      m = C.m;
      
      % Fake rectangle with slits to give to the solver.
      Prect = cell(m, 1);
      Prect{1} = polygon([1i, 0, 1, 1+1i]);
      for j = 2:m
        Prect{j} = polygon([1 2]); % Only the angles matter for now.
      end
      Prect = intpolys(Prect{:});
      
      % Circle domain for generalized rectangle.
      Cr = cell(m, 1);
      Cr{1} = {C.c(1), C.r(1), unwrap(C.t(1,genrect))};
      st = C.t(1,genrect(1)) + [0, pi];
      for j = 2:m
        Cr{j} = {C.c(j), C.r(j), st};
      end
      Cr = circdomain(Cr{:});
      
      % Args for superclass constructor are:
      %  1. Polygons object (here we'll give a fake, since that's what
      %    we're solving for).
      %  2. The discretization N variable (number of reflections).
      %  3. Circle domain.
      %  4. Objective function
      %  5. Nonlinear equation solver.
      
      args{1} = Prect;
      args{2} = opts.N;
      args{3} = Cr;
      
      % Instantiate objective function.
      args{4} = rstobjfun(Prect, Cr, genrect, opts);
      
      % Instantiate numerical solver.
      args{5} = numcontin(args{4}, C.Xu);
      args{5}.tolerance = opts.tolerance;
    end
     
    % Superclass constructor calls solver if args is not empty.
    F = F@mcscmap(args{:});
    
    if nargin
      % Post solver operations.
    end
  end
  
  function [arc, line, k, j] = integpath(F, z)
    % Integral path for map evaluations.
    
    % Basically a copy of what's done in intmap.
    vc = F.P.vc;
    tkj = F.C.t;
    c = F.C.c;
    r = F.C.r;
    
    [j, t] = nearestcircle(F, z);
    
    % distance (and direction to nearest prevertex)
    dt = tkj(1:vc(j),j) - t;
    dt = mod(dt+pi, 2*pi) - pi; % dt in [-pi,pi)
    [mindt, k] = min(abs(dt));
    
    if mindt < 100*eps
      % no arc needed (z is the prevertex or at the same angle on a
      % line to the circle center)
      arc = [];
      line.lpvn = sum(vc(1:j-1)) + k;
    else
      arc.left = tkj(k,j);
      arc.lpvn = sum(vc(1:j-1)) + k;
      arc.right = tkj(k,j) - dt(k);
      line.lpvn = 0;
    end
    
    if (j > 1 && (abs(z-c(j)) < r(j) + 100*eps(r(j)))) ...
        || (j == 1 && abs(z-c(j)) > r(j) + 100*eps(r(j)))
      % call it on the circle
      line = [];
    else
      line.left = c(j) + r(j)*exp(1i*t);
    end
  end
end

end