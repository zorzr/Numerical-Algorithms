function x = gauss_elimination(A, v)
%  Gaussian Elimination: direct method for solving linear systems.
%  Partial pivoting is used. Similar to LU factorization.
%  Computational cost:  O(2/3 n^3)
    U = A;
    b = v;
    n = length(b);
    
    for k = 1:n-1
        u = U(k:n,k);
        j = find(abs(u) == max(abs(u))) + k - 1;
        
        tmp = U(j, k:n);
        U(j, k:n) = U(k, k:n);
        U(k, k:n) = tmp;
        
        tmp = b(j);
        b(j) = b(k);
        b(k) = tmp;
        
        for i = k+1:n
            I = U(i,k) / U(k,k);
            U(i,k:n) = U(i,k:n) - I*U(k,k:n);
            b(i) = b(i) - I*b(k);
        end
    end
    
    x = upper_triangle(U, b);
end
