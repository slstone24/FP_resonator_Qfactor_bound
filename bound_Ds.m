function [cvx_optval,popt] = bound_Ds(S, D, A, epsr, U)
    % optimization problem:
    % max. p'*Im(w'*G0)*p + p'*Im(xi)p
    % s.t. p'*Re{D*S}*p = 0
    %      p'*p = Re(chi)
    % S = G0 + xi_Mat
    % A = Im(w'*G0) + Im(xi)
    % init
    D = D(~cellfun('isempty',D)); % remove empty D entry
    ND = length(D);
    
    % construct matrices for SDP
    C = cell(1,ND);
    for i = 1:ND
        C{i} = U'*Mat_real(D{i}*(S))*U;
    end
    
    % solve SDP by cvx
    A = U'*A*U;
    C{ND+1} = U'*U;
    n = size(A,1);
    
    cvx_begin quiet
        variable X(n,n) hermitian
        maximize( real(A(:)'*X(:)) );
        subject to
            for i = 1:ND
                C{i}(:)'*X(:) == 0;
            end
            C{ND+1}(:)'*X(:) == real(epsr-1);
            X == hermitian_semidefinite(n);
    cvx_end
    
    popt = extract_p_opt(X);
end
