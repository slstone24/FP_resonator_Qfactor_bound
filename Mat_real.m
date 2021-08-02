function [ReX, ImX] = Mat_real(X)
    % notice: symmetric and unsymmetric part of matrix X
    ReX = 1/2 * (X + X');
    ImX = 1/(2j) * (X - X');
end
