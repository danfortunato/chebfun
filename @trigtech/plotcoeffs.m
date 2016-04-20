function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display Trigonometric coefficients graphically.
%   PLOTCOEFFS(F) plots the Trigonometric coefficients of a TRIGTECH F on a
%   semilogy scale.  If F is an array-valued TRIGTECH then a curve is plotted
%   for each component (column) of F.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc, in the standard MATLAB manner. If S contains a string
%   'LOGLOG', the coefficients will be displayed on a log-log scale.
%
%   H = PLOTCOEFFS(F) returns a column vector of handles to lineseries
%   objects.
%
%   Note: to make the COEFPLOT easier to read, zero coefficients have a small
%   value added to them (typically EPS*VSCALE(F)) so that the curve displayed
%   is continuous.
%
% See also PLOT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% get linestyle defaults
curDALSO = get(groot,'DefaultAxesLineStyleOrder');
curDLMS = get(groot,'DefaultLineMarkerSize');

% set the linestyle defaults
set(groot,'DefaultAxesLineStyleOrder','.');
set(groot,'DefaultLineMarkerSize',10);

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        varargout{1} = plot([]);
    end
    return
end

% Set defaults:
loglogPlot = false;

% Copy input arguments:
args = varargin;

% Check inputs for 'loglog'.
j = 1;
while ( j <= length(args) )
    if ( strcmpi(args{j}, 'loglog') )
        loglogPlot = true; 
        args(j) = [];
    else
        j = j + 1;
    end
end

% Store the hold state of the current axis:
holdState = ishold;

% The coefficients:
absCoeffs = abs(f.coeffs);

% Get the size:
[n, m] = size(absCoeffs);

% Need to handle odd/even cases separately
isEven = ~mod(n, 2);
if ( isEven )
    coeffIndex = -n/2:n/2-1;
else
    coeffIndex = -(n-1)/2:(n-1)/2;
end

if ( ~loglogPlot )
    % Plot the coefficients:
    h = semilogy(coeffIndex, absCoeffs, args{:});
    if ( ~holdState )
        xlim([min(coeffIndex(1),-1) -min(coeffIndex(1),-1)]);
    end
    % Set the string for the x-axis label.
    xlabelStr = 'Wave number';
else
    if ( isEven )
        % In this case the negative cofficients have an additional term
        % corresponding to the cos(N/2*x) coefficient. We will store
        % the constant mode coefficient in both vectors.
        cPos = absCoeffs(n/2+1:n,:);
        cNeg = absCoeffs(n/2+1:-1:1,:);
    else
        cPos = absCoeffs((n+1)/2:n,:);
        cNeg = absCoeffs((n+1)/2:-1:1,:);
    end
    coeffIndexPos = 1:size(cPos,1);
    coeffIndexNeg = 1:size(cNeg,1);

    % Plot the coefficients for the positive and negative fourier modes
    % separately.
    h = loglog([coeffIndexPos nan coeffIndexNeg], [cPos ; nan(1, m) ; cNeg], args{:});

    % Set the string for the x-axis label.  In this case we will be
    % plotting the absolute value of the wave number + 1 (since we can't
    % represent a wave number <= 0 on a logarithmic scale.
    xlabelStr = '|Wave number|+1';
end

% For constant functions, plot a dot:
if ( n == 1 )
    set(h, 'marker', 'o');
end

% Add title and labels
title(gca, 'Fourier coefficients')
xlabel(gca, xlabelStr)
ylabel(gca, 'Magnitude of coefficient')

% By default, set grid on
grid(gca, 'on')

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h;
end

% set the default linestyle to factory
set(groot,'DefaultAxesLineStyleOrder',curDALSO);
set(groot,'DefaultLineMarkerSize',curDLMS);

end
