function out = feval(f, x, y, z)
%FEVAL  Evaluate a CHEBFUN3 at one or more points.
%   FEVAL(F,X,Y,Z) evaluates the CHEBFUN3 F at the point(s) in (X,Y,Z), 
%   where X, Y, and Z are scalars, vectors or tensors of doubles.
%
%   See also SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    out = [];
    return
end

% Case where all 3 inputs are chebfuns
if ( isa(x, 'chebfun') && isa(y, 'chebfun') && isa(z, 'chebfun') )
    if ( norm(x.domain-y.domain) == -norm(y.domain-z.domain) )
    % Extract chebfun3 along the path x = x(t), y = y(t), z = z(t), i.e.,
    % to construct the 1D chebfun f(x(t), y(t), z(t)).
        out = chebfun(@(t) feval(f, x(t), y(t), z(t)), x.domain);
        return
    else
    % Throw error for domain mismatch
        error('CHEBFUN:CHEBFUN3:feval:domains',...
              'For chebfun inputs the domains of each chebfun must match.');
    end
end

% Determine the number of semicolons
nsc = 0;
if ( strcmp(x, ':') )
    nsc = nsc + 1;
end
if ( strcmp(y, ':') )
    nsc = nsc + 1;
end
if ( strcmp(z, ':') )
    nsc = nsc + 1;
end

if ( nsc == 3 )
% return F if NSC == 3
    out = f;  
    return
elseif ( nsc > 0 )
% case where output is either a chebfun or chebfun2
% in this case the numeric inputs must be scalars
    cols = f.cols;
    rows = f.rows;
    tubes = f.tubes;
    core = f.core;
    if ( ~strcmp(x, ':') && isscalar(x) )
    % x is a scalar
        core = chebfun3.txm(core,cols(x),1);
    elseif ( ~strcmp(x, ':') )
    % throw error if x is not : and not scalar
        error('CHEBFUN:CHEBFUN3:feval:input',...
              'With semicolon inputs the numeric inputs must be scalars.');
    end
    if ( ~strcmp(y, ':') && isscalar(y) )
    % y is a scalar
        core = chebfun3.txm(core,rows(y),2);
    elseif ( ~strcmp(y, ':') )
    % throw error if y is not : and not scalar
        error('CHEBFUN:CHEBFUN3:feval:input',...
              'With semicolon inputs the numeric inputs must be scalars.');
    end
    if ( ~strcmp(z, ':') && isscalar(z) )
    % z is a scalar
        core = chebfun3.txm(core,tubes(z),3);
    elseif ( ~strcmp(z, ':') )
    % throw error if z is not : and not scalar
        error('CHEBFUN:CHEBFUN3:feval:input',...
              'With semicolon inputs the numeric inputs must be scalars.');
    end

    % After we squeeze the core the result is either a matrix or a vector.
    core = squeeze(core);

    if ( nsc == 1)
    % chebfun output
        if ( strcmp(x, ':') )
        % chebfun in x
            out = cols*core;
            return    
        elseif ( strcmp(y, ':') )
        % chebfun in y
            out = rows*core;
            return    
        else
        % chebfun in z
            out = tubes*core;
            return
        end
    else
    % chebfun2 output
error('No chebfun2 support')
    end
else
% case where all 3 inputs are numeric
    % get sizes
    sx = size(x);
    sy = size(y);
    sz = size(z);
    
    if ( norm(sx-sy) ~= -norm(sy-sz) )
    % throw error if inputs are different sizes
        error('CHEBFUN:CHEBFUN3:feval:input',...
              'When x, y and z are numeric they must be the same size.');
    end

    % extract cols, rows, tubes and core from f
    cols = f.cols;
    rows = f.rows;
    tubes = f.tubes;
    core = f.core;

    % determine if values come from meshgrid or ndgrid 
    % when x, y and z are 3rd order tensors
    isndgrid = false;
    ismeshgrid = false;
    if ( length(sx) == 3 && min(sx) > 1 )
        % check symmetries along different dimensions
        D = zeros(3);
        D(1,1) = norm(reshape(diff(x,1,1),[],1));
        D(1,2) = norm(reshape(diff(x,1,2),[],1));
        D(1,3) = norm(reshape(diff(x,1,3),[],1));
        D(2,1) = norm(reshape(diff(y,1,1),[],1));
        D(2,2) = norm(reshape(diff(y,1,2),[],1));
        D(2,3) = norm(reshape(diff(y,1,3),[],1));
        D(3,1) = norm(reshape(diff(z,1,1),[],1));
        D(3,2) = norm(reshape(diff(z,1,2),[],1));
        D(3,3) = norm(reshape(diff(z,1,3),[],1));

        % check for ndgrd
        if ( norm(D-diag(diag(D))) == 0 )
            isndgrid = true;
        end

        % check for meshgrid
        D = fliplr([D(:,3),D(:,1:2)]);
        if ( norm(D-diag(diag(D))) == 0 )
            ismeshgrid = true;
        end
    end
       
    if ( isndgrid ) 
    % case when data comes from ndgrid
        x = reshape(squeeze(x(:, 1, 1)),[],1);
        y = reshape(squeeze(y(1, :, 1)),[],1);
        z = reshape(squeeze(z(1, 1, :)),[],1);
        cols = transpose(cols(x));
        rows = transpose(rows(y));
        tubes = transpose(tubes(z));
        out = chebfun3.txm(chebfun3.txm(chebfun3.txm(core,cols,1),rows,2),tubes,3);
    elseif ( ismeshgrid )
    % case when data comes from meshgrid
        x = squeeze(x(1, :, 1));
        y = squeeze(y(:, 1, 1));
        z = squeeze(z(1, 1, :));
        cols = transpose(cols(x));
        rows = transpose(rows(y));
        tubes = transpose(tubes(z));
        core = permute(core, [2 1 3]);
        out = chebfun3.txm(chebfun3.txm(chebfun3.txm(core,rows,1),cols,2),tubes,3);
    else
    % generic vector case
        x = x(:);
        y = y(:);
        z = z(:);
        out = 0*x;
        for ii = 1:length(out)
            out(ii) = squeeze(chebfun3.txm(chebfun3.txm(chebfun3.txm(...
                      core,cols(x(ii)),1),rows(y(ii)),2),tubes(z(ii)),3));
        end
        out = reshape(out,sx);
    end
end

end
