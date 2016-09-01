function f = getVals2CoeffsTransform
%GETVALS2COEFFSTRANSFORM   Returns the values to coeffs transform on the sphere.

    f = @(u) reshape(trigtech.vals2coeffs( ...
        trigtech.vals2coeffs(u).').', length(u)^2, 1);
    
end