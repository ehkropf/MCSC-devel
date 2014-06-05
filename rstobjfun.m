% RSTOBJFUN provides objective function for (rectangle) resistor problems.

% E. Kropf, 2014

classdef rstobjfun < objectivefun

properties
  PR
  C
  quad
end

methods
  function F = rstobjfun(C)
    
  end
  
  function f = fun_eval(F, Xu)
    % Calculate objective function.
    
    
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
end

end
