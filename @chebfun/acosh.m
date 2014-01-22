function F = acosh(F, varargin)
%ACOSH   Inverse hyperbolic cosine of a CHEBFUN.
%   ACOSH(F) computes the inverse hypoerbolic cosine of the CHEBFUN F.
%
%   ACOSH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Call the compose method:
F = compose(F, @acosh, varargin{:});

end
