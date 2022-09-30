function [merit_avg] = RLS_TM(N, M, n_meas, n_rep, noise_coeff, persistence, N_phase_steps, mode, lambda, reg_constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [merit_avg] = RLS_TM(N, M, n_meas, n_rep, noise_coeff, persistence, N_phase_steps, mode, lambda, reg_constant)
% Simulates light control through a scattering medium via a the online
% recursive least-squares (RLS) implementation described in our
% accompanying paper.
% The scattering medium is allowed to be dynamic and the measurement system
% affected by additive white random gaussian noise.
%
% Inputs:            N      - number of independently controllable elements
%                             to shape the beam incident onto the
%                             scattering medium
%                    M      - number of independent observable output modes
%               n_meas      - number of intensity measurements
%                n_rep      - number of repetitions
%          noise_coeff      - noise coefficient: it corresponds to
%                             1/sqrt(SNR), where SNR is the signal-to-noise
%                             ratio
%          persistence      - persistence time of the scattering medium,
%                             expressed in number of intensity measurements
%        N_phase_steps      - number of phase-stepped intensity images
%                             collected for 1 field measurement
%                 mode      - light control mode
%                             {'focusing', 'energy_transmission', ...
%                              'psf_engineering'}
%               lambda      - forgetting factor of the RLS algorithm
%         reg_constant      - regularization constant of the loss function
%                             minimized by the RLS algorithm
%
% Outputs:   merit_avg      - figure of merit for the reconstruction of the
%                             transfer matrix, averaged over all
%                             repetitions
%
% We hereby recover the transfer matrix from (n_meas / N_phase_steps) field
% measurements, parallelized over the M output pixels, yielding a total
% number of Mx(n_meas / N_phase_steps) measurements.
% We expect the figure of merit to increase (or remain constant at worst)
% with increasing number of measurements, increasing the oversampling ratio
% of the linear system.
% Note, that past measurements do not need to be stored in memory, unlike
% in the algorithm in TM.m, which makes this implementation advantageous
% when dealing with long (>> N) streams of data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    set(0, 'DefaultFigureRenderer', 'painters');
    % Total incident energy
    I_tot = 1;

    mu = 0;
    sigma = 0.8;
    % Hadamard patterns
    H = hadamard(N);
    % Standard deviation of dynamic part of the transfer matrix
    dyn_std = sqrt(N_phase_steps / persistence);

    % Fourier mask for psf engineering
    mask = zeros(sqrt(M), sqrt(M));
    [uu, vv] = meshgrid(linspace(-sqrt(M)/2, sqrt(M)/2, sqrt(M)), ...
                        linspace(-sqrt(M)/2, sqrt(M)/2, sqrt(M)));
    radius_big = max(uu(:)) * 0.9;
    m = 1;
    mask(uu.^2 + vv.^2 < radius_big^2) = 1;
    theta = atan2(vv, uu);
    mask = mask .* exp(1i * m * theta);

    merit = nan(n_rep, round(n_meas / N_phase_steps));
    for i_rep = 1 : n_rep
        p12 = 1 / sqrt(reg_constant) * eye(N);

        % Ground truth Transfer Matrix
        T = normrnd(mu, sigma, M, N) + 1i * normrnd(mu, sigma, M, N);

        switch mode
            case 'focusing'
                back = 0;
        end

        % First estimate of t
        Tc = zeros(M, N) + 1i * zeros(M, N);

        % Loop: it takes N_phase to get one field measurement
        for i_meas = 1 : floor(n_meas / N_phase_steps)
            switch mode
                case 'energy_transmission'
                    % Mean transmittance
                    x_plane = ones(N, 1);
                    % Prefactor for energy conservation
                    prefactor = sqrt(I_tot / sum(abs(x_plane).^2));
                    y_mean = prefactor * x_plane' * T';
                    % Add noise
                    noise = normrnd(0, sqrt(noise_coeff * (abs(y_mean).^2) / 2), size(y_mean, 1), size(y_mean, 2)) + 1i * ...
                            normrnd(0, sqrt(noise_coeff * (abs(y_mean).^2) / 2), size(y_mean, 1), size(y_mean, 2));
                    y_mean = y_mean + noise;
                    % Transmittance
                    trans_mean = sum(abs(y_mean).^2);
            end            
            % Measurement
            x = exp(1i * (H(:, mod(i_meas-1, N)+1)+1) * pi/2);
            % Prefactor for energy conservation
            prefactor = sqrt(I_tot / sum(abs(x).^2));
            y = prefactor * x' * T';

            % Add noise
            sigma_noise = sqrt(noise_coeff * mean(abs(y).^2) / 2);

            % Add noise
            noise = normrnd(0, sigma_noise, size(y, 1), size(y, 2)) + 1i * ...
                    normrnd(0, sigma_noise, size(y, 1), size(y, 2));
            y = y + noise;
            switch mode
                case 'focusing'
                    back = back + abs(y(round(M/2))).^2 / n_meas * N_phase_steps;
            end

            % RLS update
            [Tc, p12] = RLS_update(Tc, lambda, x, p12, y);

            % Input to check output pattern
            switch mode
                case 'focusing'
                    x_opt_test = Tc(round(M/2), :)';
                    x_opt_test = exp(1i * angle(x_opt_test));
                    % Prefactor for energy conservation
                    prefactor = sqrt(I_tot / sum(abs(x_opt_test).^2));
                    y_test = prefactor * x_opt_test' * T';
                    % Add noise
                    noise = normrnd(0, sigma_noise, size(y_test, 1), size(y_test, 2)) + 1i * ...
                            normrnd(0, sigma_noise, size(y_test, 1), size(y_test, 2));
                    y_test = y_test + noise;
                    % Enhancement
                    merit(i_rep, i_meas) = abs(y_test(round(M/2))).^2;
                case 'energy_transmission'
                    [~, ~, V] = svd(Tc);
                    x_opt_test = V(:, 1);
                    x_opt_test = exp(1i * angle(x_opt_test));
                    % Prefactor for energy conservation
                    prefactor = sqrt(I_tot / sum(abs(x_opt_test).^2));
                    y_test = prefactor * x_opt_test' * T';
                    % Add noise
                    noise = normrnd(0, sigma_noise, size(y_test, 1), size(y_test, 2)) + 1i * ...
                            normrnd(0, sigma_noise, size(y_test, 1), size(y_test, 2));
                    y_test = y_test + noise;
                    % Transmittance
                    merit(i_rep, i_meas) = sum(abs(y_test).^2) / trans_mean;
                case 'psf_engineering'
                    Tc_2d = reshape(Tc, [sqrt(M), sqrt(M), N]);
                    Tc_2d_gt = reshape(T, [sqrt(M), sqrt(M), N]);
                    Tc_filt = zeros(size(Tc_2d));
                    Tc_filt_gt = zeros(size(Tc_2d));
                    for i_n = 1 : N
                        Tc_filt(:, :, i_n) = ifftshift(fft2(ifftshift(...
                                               fftshift(fft2(fftshift(...
                                               Tc_2d(:, :, i_n)))) .* mask)));
                        Tc_filt_gt(:, :, i_n) = ifftshift(fft2(ifftshift(...
                                               fftshift(fft2(fftshift(...
                                               Tc_2d_gt(:, :, i_n)))) .* mask)));
                    end
                    Tc_filt = reshape(Tc_filt, [M, N]);
                    x_opt_test = (Tc_filt(round(M/2), :))';
                    x_opt_test = exp(1i * angle(x_opt_test));
                    % Prefactor for energy conservation
                    prefactor = sqrt(I_tot / sum(abs(x_opt_test).^2));
                    y_output = prefactor * x_opt_test' * T';
                    merit(i_rep, i_meas) = error_compute(T(:), Tc(:));
            end

            % Dynamic medium
            T = T + normrnd(0, sqrt(dyn_std^2 / 2), M, N) + 1i * ...
                    normrnd(0, sqrt(dyn_std^2 / 2), M, N);
            T = T / sqrt(1 + (dyn_std / sqrt(2) / sigma)^2);
        end
        switch mode
            case 'focusing'
                merit(i_rep, :) = merit(i_rep, :) / back;
        end
    end
    merit_avg = mean(merit, 1);
end