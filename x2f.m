function U = x2f(x,N,L)
    % transform from real space (x coordinates) to Fourier space (indexed by -N:N)
    NU = 2*N + 1;
    Nx = length(x);
    U = zeros(Nx,NU);
    for nn = -N:N
        U(:,nn+N+1) = 1/sqrt(L) * exp(1j*2*pi*nn/L*x);
    end
end