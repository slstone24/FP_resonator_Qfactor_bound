function p_opt = extract_p_opt(X)
    %popt is eigenvector corresponding to largest eigenvalue
        
    if any(isnan(X(:))) % load output with nan if X is nan
        p_opt = nan;
        return
    end

    % return largest eigenvector of X
    [p_opt,~] = eigs(X,1);
end
