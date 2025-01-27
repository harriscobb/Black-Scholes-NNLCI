function x = thomas_algorithm(a,b,c,d)
    n = length(d);
    c_prime = zeros(n-1,1);
    d_prime = zeros(n,1);

    %Forward sweep
    c_prime(1) = c(1) / b(1);
    d_prime(1) = d(1) / b(1);

    for i = 2:n-1
        denom = b(i) - a(i - 1) * c_prime(i - 1);
        c_prime(i) = c(i) / denom;
        d_prime(i) = (d(i) - a(i - 1) * d_prime(i - 1)) / denom;
    end
    denom = b(n) - a(n - 1) * c_prime(n - 1);
    d_prime(n) = (d(n) - a(n - 1) * d_prime(n - 1)) / denom;

    %Back substitution
    x(n) = d_prime(n);
    for i = n - 1:-1:1
        x(i) = d_prime(i) - c_prime(i) * x(i + 1);
    end
end
