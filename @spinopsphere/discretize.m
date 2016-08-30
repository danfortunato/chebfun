function [L, Nc] = discretize(S, N)
%DISCRETIZE   Discretize a SPINOPSPHERE.
%   [L, NC] = DISCRETIZE(S, N) uses a Fourier spectral method in coefficient 
%   space to discretize the SPINOPSPHERE S with N grid points in each direction. 
%   L is the linear part, a N^2xN^2 diagonal matrix and NC is the 
%   differentiation term of the nonlinear part (and hence is linear).
%
% Remark: DISCRETIZE will fail to discretize SPINOPSPHERE objects which are not 
%         of the right form. See HELP/SPINOPSPHERE.
%
% See also SPINOPSPHERE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Set-up:
 
% Get the linear part LFUN, the nonlinear part NFUN, and the number of variables 
% NVARS from S:
funcL = S.linearPart;
nVars = nargin(funcL);

% Get the variables of the workspace:
func = functions(funcL);
wrk = func.workspace{1};
names = fieldnames(wrk);
if ( isempty(names) == 0 )
    lengthNames = size(names, 1);
    for k = 1:lengthNames
        eval(sprintf('%s = wrk.(names{k});', names{k}));
    end
end

%% Discretize the linear part:

% Identity:
I = sparse(eye(N));

% Fourier differentiation matrices:
D = trigspec.diffmat(N, 1);
D2 = trigspec.diffmat(N, 2);

% Fourier multiplication matrices:
P = eye(m+1); 
P = P(:, 1:m); 
P(1,1) = .5; 
P(m+1,1) = .5;
Q = eye(m+1+5); 
Q = Q(4:m+3,:); Q(1,4) = .5; 
Q(1,m+4) = .5;
sin2 = chebfun(@(x) sin(x).^2, [-pi, pi], 'trig');
Msin2 = full(trigspec.multmat(m+1+5, sin2));
Msin2 = Msin2(:, 4:m+4);
Msin2 = Q*Msin2*P;
cossin = chebfun(@(x) cos(x).*sin(x), [-pi, pi], 'trig');
Mcossin = full(trigspec.multmat(m+1+5, cossin));
Mcossin = Mcossin(:, 4:m+4);
Mcossin = Q*Mcossin*P;

% Look for 'laplacian'/'lap':
strL = func2str(funcL);
isLap = isempty(strfind(strL,'laplacian')) && isempty(strfind(strL,'lap'));
isLap = ~isLap;

% NxN identity matrix for the Kronecker products:
I = eye(N);

% Construct the Laplacian operator -- needed for all the operators:
if ( isLap == 1 )
    lapmat = full((kron(I, D2 + Msin2\(Mcossin*D)) + kron(D2, Msin2\I)));
else
    lapmat = 0;
end

% Convert to a string and initialize L:
strL = func2str(funcL);
L = [];

% Get the constants A in front of the Laplacians:
str = strrep(strL, 'laplacian', '');
str = strrep(str, 'lap', '');
func = eval(str);
inputs = cell(1, nVars);
for k = 1:nVars
   inputs{k} = 1; 
end
A = feval(func, inputs{:}); 

% Compute L:
for k = 1:nVars
    L = [L; A(k)*lapmat];  %#ok<*AGROW>
end

%% Disretize the differentiation term of the nonlinear part:

% We only support no differentiation, i.e., NC = 1.
% [TODO]: Support diff order > 1.
Nc = 1;

end