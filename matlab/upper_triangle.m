function x = upper_triangle(U, b)
%  Back substitution to solve linear sistems
    n = length(b);
    x = zeros(n, 1);
    
    x(n) = b(n) / U(n,n);
    for i = n-1:-1:1
        x(i) = (b(i) - U(i, i+1:n)*x(i+1:n)) / U(i,i);
    end
end
