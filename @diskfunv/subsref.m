function varargout = subsref( F, ref )
%SUBSREF   DISKFUNV subsref.
%
%   F(t,r, 'polar') returns the values of the DISKFUNV F evaluated on the 
%   array (theta,r) in polar coords.
%   F(X,Y) or F(X,Y, 'cart')
%   returns values of the DISKFUNV F evaluated on the array
%   (X,Y) in Cartesian coords. 
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3. 
%
%  
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').


% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty DISKFUNV object. 
if ( isempty( F ) )
   varargout = {[]};
   return
end 


indx = ref(1).subs;


switch ( ref(1).type )
    
    case '.'
        if ( numel( ref ) == 1 )
            % This is a get call to get a property. 
            varargout = { get(F, indx) };
        else
            % Probably .^ or maybe .* 
            %t2 = index(2).type;
            t2 = ref(2).type; 
            if ( strcmp(t2,'.') )
                out = get(F, indx, ref(2).subs{:});
            else
                out = get(F, indx);
                out = out( ref(2).subs{:} );
            end
            if ( numel(ref) > 2 )
                varargout = {subsref(out, ref(3:end))};
            else
                varargout = { out };
            end
        end
        
    case '()'
        if ( length(indx) > 1 )
            x = indx{1}; 
            y = indx{2}; % where to evaluate
            vals = feval(F, x, y, ref(1).subs{:}); 
            varargout = { vals }; 
        else
            if ( isa(indx{1},'double') )
                if all( indx{1} == 1  )
                    varargout = F.components(1);
                elseif ( all( indx{1} == 2 ) )
                    varargout = F.components(2);
                else
                    error('CHEBFUN:DISKFUNV:subsref:index', ...
                        'DISKFUNV only contains two components');
                end
            end
        end
        
    otherwise
        error('CHEBFUN:DISKFUNV:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

% Recurse down: 
if ( numel( ref ) > 1 )
   ref(1) = []; 
   varargout = { subsref( varargout{ : }, ref ) }; 
end

end