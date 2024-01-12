% Full constrained quadratic programming
% A and B are the matrices to be matched 
% Return permutation matrix P so that P*A*P' is matched to B 

function [P_qp] = matching_full_qp(A, B)
    n = size(A, 1);
    [z, ~] = quadprog_admm(A, B, [], ones(2*n,1), zeros(n^2,1), ones(n^2,1), 40, 1.5, 1e-5, 1e-3);
    %Xhat = vec2mat(z, n)';
    %Xhat = reshape(z, [n,n])';
    mat = reshape(z.', n, []).';% replace vec2mat with reshape
    Xhat=mat';
    M = matchpairs(Xhat', -99999, 'max');
    P_qp = full(sparse(M(:, 1), M(:, 2), 1, n, n));