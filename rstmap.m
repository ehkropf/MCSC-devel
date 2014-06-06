% RSTMAP creates the rectangular resistor map.

% E. Kropf, 2014

classdef rstmap < mcscmap

properties
  
end

methods
  function F = rstmap(C)
    % Here we need to find a rectangle with horizontal slits, and the slit tip
    % preimages on a given circle domain.
    % INPUTS:
    %   C = Circle domain, which in general is a solution to a parameter problem
    %   for another polygonal domain.
    %   genrect = List of vertices on C_1 in the circle domain that define the
    %   generalized rectangle.
    
    % Instantiate objective function.
    % Instantiate numerical solver.
    
    % Superclass constructor calls solver.
    
    % Post solver operations.
  end
end

end