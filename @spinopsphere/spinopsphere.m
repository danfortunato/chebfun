classdef spinopsphere < spinoperator
%SPINOPSPHERE   Class for representing the spatial part of differential 
%operators for time-dependent PDEs on the unit sphere.
%   SPINOPSPHERE is a class for representing the spatial part S of a time-
%   dependent PDE of the form u_t = S(u) = Lu + N(u) on the unit sphere, where 
%   L is a linear operator and N is a nonlinear operator. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function S = spinopsphere(varargin)
            
            % Zero input argument:
            if ( nargin == 0 )
                return
            
            % One input argument:
            elseif ( nargin == 1 )
                if ( isa(varargin{1}, 'char') == 1 )
                    [L, N, tspan, u0] = parseInputs(varargin{1});
                    S.linearPart = L;
                    S.nonlinearPart = N;
                    S.tspan = tspan;
                    S.init = u0;  
                else
                    error('SPINOPSPHERE:constructor', ['When constructing ', ...
                        'a SPINOPSPHERE with one input argument, this ', ...
                        'argument should be a STRING or a DOUBLE.'])
                end
                
                % The domain is the unit sphere:
                S.domain = 'unit sphere';
            
            % More than one input argument:    
            else
                error('SPINOPSPHERE:constructor', 'Too many input arguments.')
            end
            
        end
        
    end

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    function [L, N, tspan, u0] = parseInputs(pdechar)
        %PARSEINPUTS   Parse inputs when using a STRING.
        %   [L, N, TSPAN, U0] = PARSEINPUTS(PDECHAR), where PDECHAR is a string,
        %   outputs two function handles L and N which represent the linear and
        %   the nonlinear parts of the PDE represented by PDECHAR, the time
        %   interval TSPAN and the initial condition U0.
        
        % Allen-Cahn equation:
        if ( strcmpi(pdechar, 'ACsphere') == 1 )
            L = @(u) 1e-2*lap(u);
            N = @(u) u - u.^3;
            tspan = [0 40];
            u0 = spherefun.sphharm(5,4);
            
        % Ginzburg-Landau equation:
        elseif ( strcmpi(pdechar, 'GLsphere') == 1 )
            L = @(u) 5e-2*lap(u);
            N = @(u) u - (1 + 1.5i)*u.*(abs(u).^2);
            tspan = [0 40];
            u0 = spherefun.sphharm(12, 7);
            
        else
            error('SPINOPSPHERE:parseInputs', 'Unrecognized PDE.')
        end

    end