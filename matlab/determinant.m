function d = determinant(A)
%  Back Substitution method to obtain a upper triangular matrix
%  The determinant is just the product of the elements on the diagonal
    U = A;
    n = length(A(1,:));
    
    for k = 1:n-1
        u = U(k:n,k);
        j = find(abs(u) == max(abs(u))) + k - 1;
        
        tmp = U(j, k:n);
        U(j, k:n) = U(k, k:n);
        U(k, k:n) = tmp;
        
        for i = k+1:n
            I = U(i,k) / U(k,k);
            U(i,k:n) = U(i,k:n) - I*U(k,k:n);
        end
    end
    
    d = U(1,1);
    for i = 2:n
        d = d * U(i,i);
    end
end
