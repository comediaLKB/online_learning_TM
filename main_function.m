function main_function()

    set(0, 'defaultTextInterpreter', 'latex');
    set(0, 'defaultAxesTickLabelInterpreter', 'latex');  
    set(0, 'defaultLegendInterpreter', 'latex');
    set(0, 'DefaultFigureRenderer', 'painters');
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'defaultFigureUnits', 'normalized', 'defaultfigureposition', [0.1, 0.1, 0.6, 0.4])

    close all
    
    N = 8*8;
    limit_phase_control_no_noise_static = (N-1) * pi/4 + 1;
    n_rep = 50;
    noise_coeff = [1e-2, 0.6, 1];
    N_phase_steps = 4;
    n_meas = round(N_phase_steps * N * 5);
    persistence = N * N_phase_steps * [1e6, 3, 1];
    M = 7*7;
    M_psf = 11*11;
    reg_constant = 1;
    i_curr = 0;

    for i_pers = 1 : numel(persistence)
        lambda = 7^(-N_phase_steps / persistence(i_pers));
        for i_noise = 1 : numel(noise_coeff)
            i_curr = i_curr + 1;
            % Focusing
            tic()
            [sbr_tm] = TM(N, M, n_meas, n_rep, noise_coeff(i_noise), ...
                       persistence(i_pers), N_phase_steps, 'focusing');
            [sbr_rls_tm] = RLS_TM(N, M, n_meas, n_rep, noise_coeff(i_noise), ...
                           persistence(i_pers), N_phase_steps, 'focusing', lambda, reg_constant);
            fprintf(['Focusing mode ' num2str(i_curr) '/' ...
                    num2str(numel(persistence)*numel(noise_coeff)) ' done.\n'])

            % Max energy transmission
            [energy_tm] = TM(N, M, n_meas, n_rep, noise_coeff(i_noise), ...
                          persistence(i_pers), N_phase_steps, 'energy_transmission');
            [energy_rls_tm] = RLS_TM(N, M, n_meas, n_rep, noise_coeff(i_noise), ...
                              persistence(i_pers), N_phase_steps, 'energy_transmission', lambda, reg_constant);
            fprintf(['Max energy mode ' num2str(i_curr) '/' ...
                    num2str(numel(persistence)*numel(noise_coeff)) ' done.\n'])

            % Point-spread-function engineering
            [err_T_tm] = TM(N, M_psf, n_meas, n_rep, noise_coeff(i_noise), ...
                         persistence(i_pers), N_phase_steps, 'psf_engineering');
            [err_T_rls_tm] = RLS_TM(N, M_psf, n_meas, n_rep, noise_coeff(i_noise), ...
                             persistence(i_pers), N_phase_steps, 'psf_engineering', lambda, reg_constant);
            fprintf(['PSF engineering mode ' num2str(i_curr) '/' ...
                    num2str(numel(persistence)*numel(noise_coeff)) ' done.\n'])

            % PLOTS - FOCUSING
            figure(1)
            curr = (i_noise-1)*numel(persistence) + i_pers;
            subplot(numel(noise_coeff), numel(persistence), curr)
            plot(linspace(0, n_meas / (N * N_phase_steps), size(sbr_tm, 2)), ...
                        sbr_tm / limit_phase_control_no_noise_static);
            hold on;
            plot(linspace(0, n_meas / (N * N_phase_steps), size(sbr_rls_tm, 2)), ...
                 sbr_rls_tm / limit_phase_control_no_noise_static);
            xlim([1 / (N * N_phase_steps), n_meas / (N * N_phase_steps)]);
            y_limits = ylim;
            if i_pers == 1
                ylim([-0.005, y_limits(2)*1.005]);
            else
                ylim([y_limits(1), y_limits(2)*1.005]);
            end
            area([N * N_phase_steps, min(N * N_phase_steps + persistence(i_pers), n_meas)] / (N * N_phase_steps), ...
                 y_limits(2)*1.1*ones(1, 2), -100, ...
                 'FaceColor', 'yellow', ...
                 'FaceAlpha', 0.2, ...
                 'EdgeColor', 'none');
            title(['$SNR$ = ' num2str(sqrt(1/noise_coeff(i_noise)), '%.2f') ...
                   ', $T_p / T_{TM}$ = ' sprintf('%0.1e', persistence(i_pers) / (N * N_phase_steps))]);

            h = gca; h.LineWidth = 1; h.FontSize = 14;
            if i_noise == numel(noise_coeff)
                xlabel('$t / T_{TM}$', 'FontSize', 14)
            else
                set(gca,'xticklabel', [])
            end
            if i_pers == 1
                ylabel('$\eta / \eta_{max}$', 'FontSize', 14)         
            end
            if curr == 1
                legend('TM', ...
                       'RLS TM', ...
                       'Location', 'southeast')
            end
            sgtitle('FOCUSING')

            % PLOTS - ENERGY TRANSMISSION
            figure(2)
            curr = (i_noise-1)*numel(persistence) + i_pers;
            subplot(numel(noise_coeff), numel(persistence), curr)
            plot(linspace(0, n_meas / (N * N_phase_steps), size(energy_tm, 2)), ...
                        energy_tm);
            hold on;
            plot(linspace(0, n_meas / (N * N_phase_steps), size(energy_rls_tm, 2)), ...
                 energy_rls_tm);
            xlim([1 / (N * N_phase_steps), n_meas / (N * N_phase_steps)]);
            y_limits = ylim;
            ylim([0.9, y_limits(2)*1.005]);
            area([N * N_phase_steps, min(N * N_phase_steps + persistence(i_pers), n_meas)] / (N * N_phase_steps), ...
                 y_limits(2)*1.1*ones(1, 2), -100, ...
                 'FaceColor', 'yellow', ...
                 'FaceAlpha', 0.2, ...
                 'EdgeColor', 'none');
            title(['$SNR$ = ' num2str(sqrt(1/noise_coeff(i_noise)), '%.2f') ...
                   ', $T_p / T_{TM}$ = ' sprintf('%0.1e', persistence(i_pers) / (N * N_phase_steps))]);

            h = gca; h.LineWidth = 1; h.FontSize = 14;
            if i_noise == numel(noise_coeff)
                xlabel('$t / T_{TM}$', 'FontSize', 14)
            else
                set(gca,'xticklabel', [])
            end
            if i_pers == 1
                ylabel('$T / \langle T \rangle$', 'FontSize', 14)         
            end
            if curr == 1
                legend('TM', 'RLS TM', ...
                       'Location', 'southeast')
            end
            sgtitle('MAX ENERGY TRANSMISSION')
            
            % PLOTS - PSF ENGINEERING
            figure(3)
            curr = (i_noise-1)*numel(persistence) + i_pers;
            subplot(numel(noise_coeff), numel(persistence), curr)
            plot(linspace(0, n_meas / (N * N_phase_steps), size(err_T_tm, 2)), ...
                        err_T_tm);
            hold on;
            plot(linspace(0, n_meas / (N * N_phase_steps), size(err_T_rls_tm, 2)), ...
                 err_T_rls_tm);
            xlim([1 / (N * N_phase_steps), n_meas / (N * N_phase_steps)]);
            y_limits = ylim;
            ylim([min(err_T_rls_tm)*0.9, 1]);
            area([N * N_phase_steps, min(N * N_phase_steps + persistence(i_pers), n_meas)] / (N * N_phase_steps), ...
                 y_limits(2)*1.1*ones(1, 2), -100, ...
                 'FaceColor', 'yellow', ...
                 'FaceAlpha', 0.2, ...
                 'EdgeColor', 'none');
            title(['$SNR$ = ' num2str(sqrt(1/noise_coeff(i_noise)), '%.2f') ...
                   ', $T_p / T_{TM}$ = ' sprintf('%0.1e', persistence(i_pers) / (N * N_phase_steps))]);

            h = gca; h.LineWidth = 1; h.FontSize = 14;
            if i_noise == numel(noise_coeff)
                xlabel('$t / T_{TM}$', 'FontSize', 14)
            else
                set(gca,'xticklabel', [])
            end
            if i_pers == 1
                ylabel('norm. error', 'FontSize', 14)         
            end
            if curr == 1
                legend('TM', 'RLS TM', '$T_{stab}$', ...
                       'Location', 'southeast')
            end
            sgtitle('PSF ENGINEERING')
            fprintf(['Estimated remaining time: ' ...
                    sprintf('%0.1f', toc()/60*(numel(persistence)*numel(noise_coeff)-i_curr)) ...
                    ' minutes. Hold on ...\n'])
        end
    end
    disp('Simulations complete.')
end
