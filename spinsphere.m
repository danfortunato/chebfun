function [uout, tout] = spinsphere(varargin)
%SPINSPHERE  Solve a time-dependent PDE on the unit sphere, using a Fourier
%spectral method and an exponential integrator time-stepping scheme.
%
%   UOUT = SPIN2SPHERE(PDECHAR) solves the PDE specified by the string PDECHAR, 
%   and plots a movie of the solution as it computes it. The space and time 
%   intervals and the initial condition are chosen to produce beautiful movies. 
%   Strings available include include 'AC2s' for Allen-Cahn equation and 'GL2s' 
%   for Ginzburg-Landau equation. The output UOUT is a SPHEREFUN corresponding 
%   to the solution at the final time (a CHEBMATRIX for systems of equations, 
%   each row representing one variable). 
%
%   UOUT = SPINSPHERE(PDECHAR, TSPAN) solves the PDE from TPSAN(1) to TSPAN(END)
%   where TSPAN=[0 T1 T2 ... TF] is a vector of time chunks. The output UOUT is 
%   a CHEBMATRIX, each row corresponding to one variable and each column to one 
%   time chunk (unless TSPAN=[0 TF] and there is only one variable, in which 
%   case the output is a CHEBFUN2 at TF).
%
%   UOUT = SPINSPHERE(PDECHAR, TSPAN, U0) solves the PDE with initial condition 
%   a SPHEREFUN U0 (one variable) or a CHEBMATRIX U0 (systems). 
%
%   UOUT = SPINSPHERE(S) solves the PDE specified by the SPINOPSPHERE S and
%   plots a movie of the solution as it computes it. See HELP/SPINOPSPHERE.
%
%   UOUT = SPINSPHERE(..., PREF) allows one to use the preferences specified by
%   the SPINPREF2S object PREF. See HELP/SPINPREFSPHERE.
% 
%   [UOUT, TOUT] = SPINSPHERE(...) also returns the times chunks TOUT at which
%   UOUT was computed.
%
% Remark 1: Available (case-insensitive) strings PDECHAR are
%
%    - 'ACsphere' for Allen-Cahn equation,
%    - 'GLsphere' for Ginzburg-Landau equation.
%
% See also SPINOP2S, SPINPREF2S, SPINSCHEME, SPIN, SPIN2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

pref = [];
j = 1;
while ( j <= nargin )
    item =  varargin{j};
    if ( isa(item, 'spinoperator') == 1 )
        if ( isa(item, 'spinop') == 1 )
            error('CHEBFUN:SPIN2S', 'Use SPIN for PDEs in one space dimension.')
        elseif ( isa(item, 'spinop2') == 1 )
            error('CHEBFUN:SPINS2', ['Use SPIN2 for PDEs in two space ', ...
                'dimensions.'])
        elseif ( isa(item, 'spinop3') == 1 )
            error('CHEBFUN:SPIN2S', ['Use SPIN3 for PDEs in three space ', ...
                'dimensions.'])
        end
    elseif ( isa(item, 'char') == 1 )        
        isDemo = spinoperator.isDemoCheck(item);
        % This is a char for a demo, e.g., 'AC2s' or 'GL2s':
        if ( isDemo == 1 )
            is2D = ~isempty(strfind(item, '2'));
            is3D = ~isempty(strfind(item, '3'));
            isSphere = ~isempty(strfind(item, 'sphere'));
            is1D = ~is2D && ~is3D && ~isSphere;
            if ( is1D == 1 )
                error('CHEBFUN:SPINSPHERE', ['Use SPIN for PDEs in one ', ...
                    'space dimension.'])
            elseif ( is2D == 1 )
                error('CHEBFUN:SPINSPHERE', ['Use SPIN2 for PDEs in two ', ...
                    'space dimensions.'])
            elseif ( is3D == 1 )
                error('CHEBFUN:SPINSPHERE', ['Use SPIN3 for PDEs in three ', ...
                    'space dimensions.'])
            end
        % This is a preference, e.g., 'N' or 'dt':
        else
            if ( isempty(pref) == 1 )
                pref = spinprefsphere();
            end
            pref.(item) = varargin{j+1};
            varargin{j} = [];
            varargin{j + 1} = [];
            j = j + 2;
            continue
        end
    end
    j = j + 1;
end

% Add the preferences:
if ( isempty(pref) == 0 )
   varargin{end + 1} = pref;
end

% Get rid of the deleted entries:
varargin = varargin(~cellfun(@isempty, varargin));

% SPIN2S is a wrapper for SOLVPDE:
[uout, tout] = spinoperator.solvepde(varargin{:});

end