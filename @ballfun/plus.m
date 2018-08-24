function h = plus(f, g)
%+   BALLFUN plus.
%   F + G adds BALLFUNs F and G, or a scalar to a BALLFUN if either F or G is a
%   scalar.
%
%   H = PLUS(F, G) is called for the syntax 'F + G'.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');

if (fIsBallfun && gIsBallfun)
    [mf,nf,pf] = size(f.coeffs); 
    [mg,ng,pg] = size(g.coeffs); 
    m = max(mf,mg); 
    n = max(nf,ng);
    p = max(pf,pg); 
    X = zeros(m,n,p); 
    X(1:mf,1:nf,1:pf) = f.coeffs; 
    X(1:mg,1:ng,1:pg) = X(1:mg,1:ng,1:pg) + g.coeffs;
    h = ballfun(X);
elseif (fIsBallfun && isnumeric(g))
    S = size(f.coeffs);
    X = f.coeffs;
    % Add the constant g
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + g;
    h = ballfun(X);
elseif (isnumeric(f) && gIsBallfun)
    S = size(g.coeffs);
    X = g.coeffs;
    % Add the constant f
    X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) = X(1,floor(S(2)/2)+1,floor(S(3)/2)+1) + f;
    h = ballfun(X);
else
    error('BALLFUN:mtimes:unknown', ...
          ['Undefined function ''plus'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
