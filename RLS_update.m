function [T, p12] = RLS_update(T, lambda, x, p12, y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T, p12] = RLS_update(T, lambda, x, p12, y)
% Implements one iteration of the Inverse QR-decomposition-based recursive
% least-squares (QRD-RLS) algorithm, see [Haykin, Adaptive filter theory,
% 5th ed., chapter 10].
%
% Inputs:            T      - current estimate of the transfer matrix
%               lambda      - forgetting factor
%                    x      - input of current data point (x, y)
%                  p12      - current estimate of the square root of the
%                             inverse covariance matrix
%                    y      - output of current data point (x, y)
%
% Outputs:           T      - new estimate of the transfer matrix
%                  p12      - new estimate of the square root of the
%                             inverse covariance matrix
%
% Propagating the square root of the inverse covariance matrix prevents
% numerical instabilities arising from the inversion and update of the
% covariance matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A11 = 1;
    A21 = zeros(numel(x), 1);
    A12 = 1 / sqrt(lambda) * x' * p12;
    A22 = 1 / sqrt(lambda) * p12;
    A = [A11, A12; ...
         A21, A22];
    [~, R] = qr(A');
    R = R';
    R11 = R(1, 1);
    R21 = R(2:end, 1);
    p12 = R(2:end, 2:end);

    % Update
    T = T + (R21 / R11 * (y - x' * T'))';
end