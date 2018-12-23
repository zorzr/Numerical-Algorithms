function [x,res,k] = gauss_seidel(A,b,tol,maxits,x0)
%  Gauss-Seidel method: iterative algotithm for solving linear systems.
%  Splitting methods represent the matrix A as:     A = M - N
%   then at each iteration they solve    Mx = b + Nx
%  Here M is a matrix with the lower triangular values of A.
%  Computational cost:  O(2 n^2) for each iteration
    x = x0;
    r = b - A*x0;
    res0 = norm(r);
    M_inv = inv(tril(A));
    
    res = {};
    for k = 1:maxits
        x = x + M_inv*r;
        r = b - A*x;
        
        res{k} = norm(r);
        if (res{k}/res0 < tol)
            break
        end
    end
    res = cell2mat(res);
end
