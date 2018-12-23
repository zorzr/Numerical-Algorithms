function B = inverse(A)
%  Best way to do the inverse of any matrix: first check if has
%   useful properties, then decide which method to use.
    [C,p] = chol(A);
    
    if (p > 0)
        B = inv_lu(A);
    else
        B = inv_chol(C);
    end 
end

function E = e(i, n)
    E = zeros(n, 1);
    E(i) = 1;
end

function B = inv_lu(A)
%  Matrix inversion using the result of the LU factorization.
%  Passing a matrix with range = 0 will produce an indefinite output.
%  Computational cost:  O(8/3 n^3)
    [L,U,P] = lu(A);
    n = length(A(1,:));
    B = zeros(n,n);
    
    for i = 1:n
        y = lower_triangle(L, P(:,i));
        B(:,i) = upper_triangle(U, y);
    end
end

function B = inv_chol(C)
%  This technique speeds up the inverse computation using the matrix symmetry.
%  Computational cost:  O(1/3 n^3)
    n = length(C(1,:));
    B = zeros(n,n);

    for i = 1:n
        y = lower_triangle(C', e(i,n));
        B(:,i) = upper_triangle(C, y);
    end
end