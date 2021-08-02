function Dopt = get_Dopt(p0,S,einc)
    
    if isnan(p0) % if optimal current is nan, return a random diagonal D matrix
        Dopt = 1j*diag(rand(size(einc)));
        fprintf('add a random D matrix as Dopt ... \n')
        return
    end
        
    % diagonal D matrix
    q = S*p0 + einc;
    d = p0.*conj(q);
    d = d/max(d);
    Dopt = diag(d); % maximally violated constraint
end
