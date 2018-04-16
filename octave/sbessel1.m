function [y,dy] = sbessel1(n,x) 
% MATLAB function to solve the spherical bessel function 
% of the first kind. 
% Martin Anderson 
% Created: 3-25-96      Last modified: 3-26-96 
% Accepts the argument(s)                       (x) (may be a vector) 
% and the required order of the solution.       (n) 
% 
% Syntax: 
%  y = jn(n,x); 

y = sqrt(pi ./ (2 * x)) .* besselj(n+0.5, x); 

if ( nargout > 1 ) 
  h = .001; 
  h2 = 2 * h; 
  len = length(x); 
  dy = zeros(size(x)); 
  for c = 1:len 
    ord = x(c)-h2:h:x(c)+h2; 
    f = sqrt(pi ./ (2 * ord)) .* besselj(n+0.5, ord); 
    F = (f(4) - f(2)) / h2; 
    F2 = (f(5) - f(1))/ (2 * h2); 
    dy(c) = F + (F - F2)/3; 
  end 
end

  

