function f = getCoeffs2ValsTransform
%GETCOEFFS2VALSTRANSFORM   Returns the coeffs to values transform on the sphere.

    f = @(u) trigtech.coeffs2vals(trigtech.coeffs2vals( ...
        reshape(u, sqrt(length(u)), sqrt(length(u))).').');
    
end