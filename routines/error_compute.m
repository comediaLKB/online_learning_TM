function error = error_compute(x_true, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error = error_compute(x_true, x)
% Computes the RMS error metric between the ground truth x_true and the
% estimate x.
% It allows for multiplication by a constant and for a constant phase
% offset.
%
% Inputs:       x_true      - ground truth signal
%                    x      - current esimtate of x_true
%
% Outputs:       error      - normalized error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gamma = sum(x_true .* conj(x)) / sum((abs(x)).^2);
    error = sum((abs(x_true - gamma * x)).^2) / sum((abs(x_true)).^2);
end