classdef fun % (Abstract) 
%FUN   Approximate functions on arbitrary domains.
%   Abstract (interface) class for approximating functions on the arbitrary 
%   intervals.
%
% See also DELTAFUN, CLASSICFUN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUN Class Description:
%  [TODO]
%
% Class diagram: [<<CHEBFUN>>] <>-- [<<FUN>>] <----[<<classicfun>>]
%                                             <----[    deltafun  ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

    methods (Static)
        function obj = constructor(op, data, pref)
            % Parse inputs.
            if ( nargin < 1 )
                % We can't return an empty FUN, so pass an empty OP down.
                op = [];
            end

            if ( (nargin < 2) || isempty(data) )
                    data = struct();
            end

            if ( (nargin < 3) || isempty(pref) )
                pref = chebfunpref();
            else
                pref = chebfunpref(pref);
            end

            % TODO:  Allow construction of DELTAFUNs as well?
            if ( isa(op, 'fun') )
                % OP is already a FUN!
                obj = op;
            else
                obj = classicfun.constructor(op, data, pref);
            end
        end
    end

    %% ABSTRACT STATIC METHODS REQUIRED BY THIS CLASS.
    methods ( Abstract = true, Static = true )

        % Map from [-1, 1] to the domain of the FUN.
        m = createMap(domain);  
        
        % Make a FUN. (Constructor shortcut)
        f = make(varargin);
    end
    
    methods ( Abstract = true, Static = false )
        
    end
    
end
