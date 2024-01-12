function [ areas ] = calc_tri_areas( M )
%CALC_TRI_AREAS Calculates the area of each triangle in M.
%   M needs to be struct of M.VERT and M.TRIV.
%   areas - n x 1
% written by Zorah Lähner (laehner@in.tum.de)

    a = sqrt(sum((M.VERT(M.TRIV(:,1), :) - M.VERT(M.TRIV(:,2), :)).^2, 2));
    b = sqrt(sum((M.VERT(M.TRIV(:,2), :) - M.VERT(M.TRIV(:,3), :)).^2, 2));
    c = sqrt(sum((M.VERT(M.TRIV(:,3), :) - M.VERT(M.TRIV(:,1), :)).^2, 2));
    p = (a + b + c) / 2;
    
    areas = sqrt(p.*(p-a).*(p-b).*(p-c));
end

