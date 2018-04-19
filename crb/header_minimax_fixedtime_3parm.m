% Header File -- Min-max CRLB analysis, with TR/flip angle variation
% Worst-case M0s/T1/T2 CRLB over T1/T2 ROI, varying TRs and flips
% Instead of TR-compensation, use fixed total imaging time
%
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014
% Version Two: Here we don't consider estimating kappa, but still consider
% robustness to flip angle variation.

% Set imaging constraints
TR_tot = 11.8*2;%41.9;                      % ms (sum of 2*TRs_min + TRd_min)
TRs_min = 11.8;%12.2;                     % ms (measured: 12.24)
TRd_min = 50;%17.5;                     % ms (measured: 17.47)
TRs_res = 0.1;                      % TR SPGR resolution
TRd_res = 0.1;                      % TR DESS resolution

% Implicit parameters, including a T1/T2/kappa range of interest
% Inner minimax problem: tight bounds
% T1 = linspace(600, 1200, 5);
T1 = linspace(800, 1400, 5);        % WM/GM voxels of interest
% T2 = linspace(50, 100, 5);      
T2 = linspace(50, 120, 5);          % WM/GM voxels of interest
kap = linspace(0.9, 1.1, 5);        % a.u. (10% variation)
wf = 0;                             % rad/ms

% Outer robustness criterion: check minima over loose bounds
% T1_rob  = linspace(600, 1500, 6);   % ms 
T1_rob  = linspace(400, 4000, 9);   % ms
T2_rob  = linspace(40,  200,  9);   % ms
kap_rob = 2 .^ linspace(-1, 1, 5);  % a.u. (factor-of-two variation)

% Rescale measured noise variance to current voxel size
oldVoxelVol = 1 * 1 * 5;                        % mm^3
newVoxelVol = (240/256) * (240/256) * (30/6);   % mm^3
oldvar_im_1coil = 2.62e-7;                      % a.u.; at 1x1x5 resolution 
newvar_im_1coil = oldvar_im_1coil * (oldVoxelVol / newVoxelVol);

% Other constant declarations
TE = 4.67;                          % ms 
noise_var_ssos = newvar_im_1coil/2; % High-SNR regime approximation

a = 0.1;                            % T1 relative importance parameter
b = 1;                              % T2 relative importance parameter

tol = 0.01;                         % Global minima tolerance
time_comp = 0;                      % Toggle time compensation on/off
pr = 0;                             % Toggle print on/off
savedat = 0;                        % Toggle saving .mat files on/off

T1range = [0 100];                  % sigmaT1 plot range
T2range = [0 10];                   % sigmaT2 plot range
Psirange = a*T1range + b*T2range;   % Psi plot range

% Choose a scan profile: (0,2), (1,1), or (2,1)
n_spgr = 2;                         % Number of SPGR scans
n_dess = 0;                         % Number of DESS scans
profile = sprintf('%u SPGR, %u DESS', n_spgr, n_dess);

% Introduce the experiment and start timer.
fprintf('\n\nScan Profile: (%u SPGR, %u DESS)\n', n_spgr, n_dess);
fprintf('Total TR constraint: %1.2f ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f,%1.1f)\n', 0, a, b);
tic; 

switch profile   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Sequence One: 0 SPGR, 2 DESS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '0 SPGR, 2 DESS'
        % Controllable parameter grid: varying TRs and flip angles
        TRd_max = TR_tot - TRd_min;
        TR1_dess = [TRd_min : TRd_res : TRd_max];
        dess1 = linspace(5,90,18)*pi/180;   % rad
        dess2 = linspace(5,90,18)*pi/180;   % rad

        % Grid definitions 
        sig_T1_02   = inf(length(TR1_dess), length(dess1), length(dess2));
        sig_T2_02   = inf(length(TR1_dess), length(dess1), length(dess2));
        sigw_T1_02  = inf(length(TR1_dess), length(dess1), length(dess2));
        sigw_T2_02  = inf(length(TR1_dess), length(dess1), length(dess2));
        Psi_02      = inf(length(TR1_dess), length(dess1), length(dess2));

        % Covariance matrix: scaled identity with M = 4
        Sigma_inv = (1 ./ noise_var_ssos) * speye(4);

        %% Step One: Min-max CRLB analysis over parameter grid
        spgr_02 = []';
        TRs_02 = []';
        for i_d1 = 1:length(dess1)
        for i_d2 = 1:length(dess2)
            dess_02 = [dess1(i_d1) dess2(i_d2)]';

            for i_TRd1 = 1:length(TR1_dess)
                % Second DESS TR is already specified from first + constraint
                TRd_02 = [TR1_dess(i_TRd1) TR_tot-TR1_dess(i_TRd1)];

                % Inner maximization: store max SDs over tight (T1,T2,kap) range
                worst_sigT1  = NaN(length(T1), length(T2), length(kap));
                worst_sigT2  = NaN(length(T1), length(T2), length(kap));

                for i_t1 = 1:length(T1)
                for i_t2 = 1:length(T2)
                for i_kap = 1:length(kap)
                    [~, ~, worst_sigT1(i_t1, i_t2, i_kap),... 
                        worst_sigT2(i_t1, i_t2, i_kap)]...
                        = norm_crlb_dess_3parm(T1(i_t1), T2(i_t2), wf,...
                        kap(i_kap)*spgr_02, kap(i_kap)*dess_02,...
                        TRs_02, TRd_02, TE, Sigma_inv, time_comp);
                end
                end
                end

                % Store the worst-case sig_T1 and sig_T2 values
                sigw_T1_02(i_TRd1, i_d1, i_d2) = max(worst_sigT1(:));
                sigw_T2_02(i_TRd1, i_d1, i_d2) = max(worst_sigT2(:));

                % Store the worst-case Psi and its corresponding (T1, T2) index
                worst_Psi = a*worst_sigT1 + b*worst_sigT2; 
                [Psi_02(i_TRd1, i_d1, i_d2), idx_Psiw] = max(worst_Psi(:));

                % Extract this index and use it to save corresponding sigmas
                [T1_idx, T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw); 
                sig_T1_02(i_TRd1, i_d1, i_d2)  = worst_sigT1(T1_idx, T2_idx, kap_idx);
                sig_T2_02(i_TRd1, i_d1, i_d2)  = worst_sigT2(T1_idx, T2_idx, kap_idx);
            end
        end
        end

        % Find the indices of (multiple) Psi_02 minima, to within a tolerance
        Psimin_idx_02 = find( (Psi_02 - min(Psi_02(:))) ./ min(Psi_02(:)) <= tol);
        num_min_02 = length(Psimin_idx_02);

        %% Step Two: Select one minimum based on robustness over broad (T1, T2, kap) range
        % Indices of each of the "num_min_02" total minima
        TRd1_idx_02 = NaN(num_min_02, 1);
        d1_idx_02  = NaN(num_min_02, 1);
        d2_idx_02  = NaN(num_min_02, 1);

        sigT1_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));
        sigT2_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));

        for i_min = 1:num_min_02
            % Convert the 1D index to ND-grid indices
            [TRd1_idx_02(i_min), d1_idx_02(i_min), d2_idx_02(i_min)]...
                = ind2sub(size(Psi_02), Psimin_idx_02(i_min));

            % Store the minimizing parameters
            TRd1min_02 = TR1_dess(TRd1_idx_02(i_min));
            TRd2min_02 = TR_tot - TRd1min_02;
            d1min_02   = dess1(d1_idx_02(i_min));
            d2min_02   = dess2(d2_idx_02(i_min));

            % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
            for i_t1 = 1:length(T1_rob)
            for i_t2 = 1:length(T2_rob)
            for i_kap = 1:length(kap_rob)
                [~, ~, sigT1_worst_02(i_min, i_t1, i_t2, i_kap),...
                    sigT2_worst_02(i_min, i_t1, i_t2, i_kap)] = ...
                    norm_crlb_dess_3parm(T1_rob(i_t1), T2_rob(i_t2), wf,...
                    kap_rob(i_kap)*[]', kap_rob(i_kap)*[d1min_02 d2min_02]', []',...
                    [TRd1min_02 TRd2min_02]', TE, Sigma_inv, time_comp);
            end
            end
            end
        end

        % Store the worst-case sigT1 and sigT2
        sigT1_min_02 = max(reshape(sigT1_worst_02,...
            [num_min_02 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        sigT2_min_02 = max(reshape(sigT2_worst_02,...
            [num_min_02 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Compute the worst-case Psi_02 value over the wider range
        Psi_min_02 = a*sigT1_worst_02 + b*sigT2_worst_02;
        [Psi_maxmin_02, wide_idx_02] = max(reshape(Psi_min_02,...
            [num_min_02 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

        % Extract the parameter indices that minimize Psi_maxmin_02 over minima
        % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
        [Psi_star_02, i_star_02] = min(Psi_maxmin_02);
        sigT1_star_02 = sigT1_min_02(i_star_02);
        sigT2_star_02 = sigT2_min_02(i_star_02);
        
        TRd1star_02 = TRd1_idx_02(i_star_02);
        d1star_02   = d1_idx_02(i_star_02);
        d2star_02   = d2_idx_02(i_star_02);

        % Display results 
        fprintf('Selected scan design: (%1.1f ms, %1.1f deg); (%1.1f ms, %1.1f deg)\n',...
            TR1_dess(TRd1star_02), dess1(d1star_02)*180/pi,...
            TR_tot-TR1_dess(TRd1star_02), dess2(d2star_02)*180/pi);
        fprintf('\t Worst-case sigT1: %6.2f\n',...
            sigw_T1_02(TRd1star_02, d1star_02, d2star_02));
        fprintf('\t Robust-range sigT1: %6.2f\n', sigT1_star_02);
        fprintf('\t Worst-case sigT2: %6.3f\n',...
            sigw_T2_02(TRd1star_02, d1star_02, d2star_02));
        fprintf('\t Robust-range sigT2: %6.2f\n', sigT2_star_02);
        fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
            a, sig_T1_02(TRd1star_02, d1star_02, d2star_02),...
            b, sig_T2_02(TRd1star_02, d1star_02, d2star_02),...
            Psi_02(TRd1star_02, d1star_02, d2star_02));
        fprintf('\t Robust-range Psi: %6.3f\n', Psi_star_02);
        fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
            num_min_02, 100*tol, toc);

        % Save useful outputs
        suffix = sprintf('_%2.2fms_%1.1f,%d', TR_tot, a, b);
        suffix = strrep(suffix, ',', '_');
        suffix = strrep(suffix, '.', 'p');
        if (savedat)
            cd(strcat('TRtot', suffix)); 
            mkdir('(0,2) Scan Profile'); cd('(0,2) Scan Profile');
            label_02 = strcat('0SPGR,2DESS', suffix, '_3parm', '.mat');
            save(label_02, 'sig_T1_02', 'sig_T2_02', 'sigw_T1_02',...
                'sigw_T2_02', 'Psi_02', 'a', 'b', 'Sigma_inv',...
                'TR1_dess', 'dess1', 'dess2', 'T1', 'T2');
        end

        %% Graphical Output: 2D Heatmaps for 0 SPGR, 2 DESS
        % Minimizer location(s)
        if (TR1_dess(TRd1star_02) == TR_tot/2)
            xmin = [dess1(d1star_02) dess2(d2star_02)] * 180/pi; 
            ymin = [dess2(d2star_02) dess1(d1star_02)] * 180/pi;
        else 
            xmin = [dess1(d1star_02)] * 180/pi;
            ymin = [dess2(d2star_02)] * 180/pi;
        end
        
        % T1 sigma vs. (dess1/dess2);
        figure; hold on; imagesc(180*dess1/pi, 180*dess2/pi, ...
            squeeze(sigw_T1_02(TRd1star_02, :, :))', T1range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(dess1*180/pi)' minmax(dess2*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_vs_d1d2_02', suffix, '.eps')), end;
        
        % T2 sigma vs. (dess1/dess2);
        figure; hold on; imagesc(180*dess1/pi, 180*dess2/pi, ...
            squeeze(sigw_T2_02(TRd1star_02, :, :))', T2range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(dess1*180/pi)' minmax(dess2*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_vs_d1d2_02', suffix, '.eps')), end;

        % Worst-case Psi vs. (dess1/dess2);
        figure; hold on; imagesc(180*dess1/pi, 180*dess2/pi, ...
            squeeze(Psi_02(TRd1star_02, :, :))', Psirange);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(dess1*180/pi)' minmax(dess2*180/pi)']);
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_vs_d1d2_02', suffix, '.eps')), end;
            if (savedat), cd ../.., end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Sequence Two: 1 SPGR, 1 DESS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '1 SPGR, 1 DESS'
        % Controllable parameter grid: varying TRs and flip angles 
        TRs_max = TR_tot - TRd_min;
        TR1_spgr = [TRs_min : TRs_res : TRs_max];
        spgr1 = linspace(5,90,18)*pi/180;   % rad
        dess1 = linspace(5,90,18)*pi/180;   % rad
        
        % Grid definitions
        sig_T1_11   = inf(length(TR1_spgr), length(spgr1), length(dess1));
        sig_T2_11   = inf(length(TR1_spgr), length(spgr1), length(dess1));
        sigw_T1_11  = inf(length(TR1_spgr), length(spgr1), length(dess1));
        sigw_T2_11  = inf(length(TR1_spgr), length(spgr1), length(dess1));
        Psi_11      = inf(length(TR1_spgr), length(spgr1), length(dess1));
        
        % Covariance matrix: scaled identity with M = 3
        Sigma_inv = (1 ./ noise_var_ssos) * speye(3);
        
        %% Step One: Min-max CRLB analysis over parameter grid
        for i_s1 = 1:length(spgr1)
            spgr_11 = [spgr1(i_s1)]';
            
            for i_d1 = 1:length(dess1)
                dess_11 = [dess1(i_d1)]';
                
                for i_TRs1 = 1:length(TR1_spgr)
                    TRs_11 = [TR1_spgr(i_TRs1)]';
                    TRd_11 = TR_tot - sum(TRs_11);
                    
                    % Inner maximization: store max SDs over (T1,T2,kap) range
                    worst_sigT1 = NaN(length(T1), length(T2), length(kap));
                    worst_sigT2 = NaN(length(T1), length(T2), length(kap));
                    
                    for i_t1 = 1:length(T1)
                    for i_t2 = 1:length(T2)
                    for i_kap = 1:length(kap)
                        [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                            worst_sigT2(i_t1, i_t2, i_kap)]...
                            = norm_crlb_dess_3parm(T1(i_t1), T2(i_t2), wf,...
                            kap(i_kap)*spgr_11, kap(i_kap)*dess_11,...
                            TRs_11, TRd_11, TE, Sigma_inv, time_comp);
                    end
                    end
                    end
                    
                    % Store the worst-case sig_T1 and sig_T2 values
                    sigw_T1_11(i_TRs1, i_s1, i_d1) = max(worst_sigT1(:));
                    sigw_T2_11(i_TRs1, i_s1, i_d1) = max(worst_sigT2(:));
                    
                    % Store the worst-case Psi and its corresponding (T1, T2) index
                    worst_Psi = a*worst_sigT1 + b*worst_sigT2;
                    [Psi_11(i_TRs1, i_s1, i_d1), idx_Psiw] = max(worst_Psi(:));
                    
                    % Extract this index and use it to save corresponding sigmas
                    [T1_idx, T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw); 
                    sig_T1_11(i_TRs1, i_s1, i_d1)  = worst_sigT1(T1_idx, T2_idx, kap_idx);
                    sig_T2_11(i_TRs1, i_s1, i_d1)  = worst_sigT2(T1_idx, T2_idx, kap_idx);
                end
            end
        end
        
        % Find the indices of (multiple) Psi_11 minima, to within a tolerance
        Psimin_idx_11 = find( (Psi_11 - min(Psi_11(:))) ./ min(Psi_11(:)) <= tol);
        num_min_11 = length(Psimin_idx_11);
        
        %% Step Two: Select one minimum based on robustness over broad (T1, T2, kap) range
        % Indices of each of the "num_min_11" total minima
        TRs1_idx_11 = NaN(num_min_11, 1);
        s1_idx_11 = NaN(num_min_11, 1);
        d1_idx_11 = NaN(num_min_11, 1);
        
        sigT1_worst_11  = NaN(num_min_11, length(T1_rob), length(T2_rob), length(kap_rob));
        sigT2_worst_11  = NaN(num_min_11, length(T1_rob), length(T2_rob), length(kap_rob));
        
        for i_min = 1:num_min_11
            % Convert the 1D index to ND-grid indices
            [TRs1_idx_11(i_min), s1_idx_11(i_min), d1_idx_11(i_min)]...
                = ind2sub(size(Psi_11), Psimin_idx_11(i_min));

            % Store the minimizing parameters
            TRs1min_11 = TR1_spgr(TRs1_idx_11(i_min));
            TRd1min_11 = TR_tot - TRs1min_11;
            s1min_11 = spgr1(s1_idx_11(i_min));
            d1min_11 = dess1(d1_idx_11(i_min));

            % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
            for i_t1 = 1:length(T1_rob)
            for i_t2 = 1:length(T2_rob)
            for i_kap = 1:length(kap_rob)
                [~, ~, sigT1_worst_11(i_min, i_t1, i_t2, i_kap),...
                    sigT2_worst_11(i_min, i_t1, i_t2, i_kap)] = ...
                    norm_crlb_dess_3parm(T1_rob(i_t1), T2_rob(i_t2),...
                    wf, kap_rob(i_kap)*[s1min_11]', kap_rob(i_kap)*[d1min_11]',...
                    [TRs1min_11]', [TRd1min_11]', TE, Sigma_inv, time_comp);
            end
            end
            end
        end
        
        % Store the worst-case sigT1 and sigT2
        sigT1_min_11 = max(reshape(sigT1_worst_11,...
            [num_min_11 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        sigT2_min_11 = max(reshape(sigT2_worst_11,...
            [num_min_11 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Compute the worst-case Psi_11 value over the wider range
        Psi_min_11 = a*sigT1_worst_11 + b*sigT2_worst_11;
        [Psi_maxmin_11, wide_idx_11] = max(reshape(Psi_min_11,...
            [num_min_11 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Extract the parameter indices that minimize Psi_maxmin_11 over minima
        % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
        [Psi_star_11, i_star_11] = min(Psi_maxmin_11);
        sigT1_star_11 = sigT1_min_11(i_star_11);
        sigT2_star_11 = sigT2_min_11(i_star_11);
        
        TRs1star_11 = TRs1_idx_11(i_star_11);
        s1star_11 = s1_idx_11(i_star_11);
        d1star_11 = d1_idx_11(i_star_11);
        
        % Display results
        fprintf('Selected scan design: (%1.1f ms, %1.1f deg); (%1.1f ms, %1.1f deg)\n',...
            TR1_spgr(TRs1star_11), spgr1(s1star_11)*180/pi,...
            TR_tot - TR1_spgr(TRs1star_11), dess1(d1star_11)*180/pi);
        fprintf('\t Worst-case sigT1: %6.2f\n',...
            sigw_T1_11(TRs1star_11, s1star_11, d1star_11));
        fprintf('\t Robust-range sigT1: %6.2f\n', sigT1_star_11);
        fprintf('\t Worst-case sigT2: %6.3f\n',...
            sigw_T2_11(TRs1star_11, s1star_11, d1star_11));
        fprintf('\t Robust-range sigT2: %6.2f\n', sigT2_star_11);
        fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
            a, sig_T1_11(TRs1star_11, s1star_11, d1star_11),...
            b, sig_T2_11(TRs1star_11, s1star_11, d1star_11),...
            Psi_11(TRs1star_11, s1star_11, d1star_11));
        fprintf('\t Robust-range Psi: %6.3f\n', Psi_star_11);
        fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
            num_min_11, 100*tol, toc);
        
        % Save useful outputs
        suffix = sprintf('_%2.2fms_%1.1f,%d', TR_tot, a, b);
        suffix = strrep(suffix, ',', '_');
        suffix = strrep(suffix, '.', 'p');
        if (savedat)
            cd(strcat('TRtot', suffix)); 
            mkdir('(1,1) Scan Profile'); cd('(1,1) Scan Profile');
            label_11 = strcat('1SPGR,1DESS', suffix, '_3parm', '.mat');
            save(label_11, 'sig_T1_11', 'sig_T2_11', 'sigw_T1_11',...
                'sigw_T2_11', 'Psi_11', 'a', 'b', 'Sigma_inv',...
                'TR1_spgr', 'spgr1', 'dess1', 'T1', 'T2');
        end
        
        %% Graphical Output: 2D Heatmaps for 1 SPGR, 1 DESS
        % Minimizer locations(s)
        xmin = [spgr1(s1star_11)*180/pi];
        ymin = [dess1(d1star_11)*180/pi];
        
        % T1 sigma vs. (spgr1/dess1);
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(sigw_T1_11(TRs1star_11, :, :))', T1range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_vs_s1d1_11', suffix, '.eps')), end;
        
        % T2 sigma vs. (spgr1/dess1);
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(sigw_T2_11(TRs1star_11, :, :))', T2range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_vs_s1d1_11', suffix, '.eps')), end;
        
        % Worst-case Psi vs. (spgr1/dess1);
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(Psi_11(TRs1star_11, :, :))', Psirange);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_vs_s1d1_11', suffix, '.eps')), end;
            if (savedat), cd ../.., end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Scan Sequence Three: 2 SPGR, 1 DESS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '2 SPGR, 1 DESS'
        % Controllable parameter grid: varying TRs and flip angles
        TRs_max = TR_tot - TRs_min - TRd_min;
        TR1_spgr = [TRs_min : TRs_res : TRs_max];
        TR2_spgr = [TRs_min : TRs_res : TRs_max];
        spgr1 = linspace(5,90,18)*pi/180;   % rad
        spgr2 = linspace(5,90,18)*pi/180;   % rad
        dess1 = linspace(5,90,18)*pi/180;   % rad

        % Grid definitions
        sig_T1_21   = inf(length(TR1_spgr), length(TR2_spgr), length(spgr1), length(spgr2), length(dess1));
        sig_T2_21   = inf(length(TR1_spgr), length(TR2_spgr), length(spgr1), length(spgr2), length(dess1));
        sigw_T1_21  = inf(length(TR1_spgr), length(TR2_spgr), length(spgr1), length(spgr2), length(dess1));
        sigw_T2_21  = inf(length(TR1_spgr), length(TR2_spgr), length(spgr1), length(spgr2), length(dess1));
        Psi_21      = inf(length(TR1_spgr), length(TR2_spgr), length(spgr1), length(spgr2), length(dess1));

        % Covariance matrix: scaled identity with M = 4
        Sigma_inv = (1 ./ noise_var_ssos) * speye(4);

        %% Step One: Min-max CRLB analysis over parameter grid
        for i_s1 = 1:length(spgr1)
        for i_s2 = 1:length(spgr2)
            spgr_21 = [spgr1(i_s1) spgr2(i_s2)]';

            for i_d1 = 1:length(dess1)
                dess_21 = [dess1(i_d1)]';

                for i_TRs1 = 1:length(TR1_spgr)
                for i_TRs2 = 1:length(TR2_spgr)
                    TRs_21 = [TR1_spgr(i_TRs1) TR2_spgr(i_TRs2)]';
                    TRd_21 = TR_tot - sum(TRs_21);

                    % Only consider points for which TRd_21 >= TRd_min
                    % Otherwise, do nothing, so Psi_21 remains infinite
                    if (TRd_21 >= TRd_min)
                        % Exploit redundancy: if we have already recorded values
                        % when the SPGR indices are swapped, just copy them!
                        if (isfinite(Psi_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1)))
                            Psi_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...
                                = Psi_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1);
                            sig_T1_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...
                                = sig_T1_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1);
                            sig_T2_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...
                                = sig_T2_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1);
                            sigw_T1_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...
                                = sigw_T1_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1);
                            sigw_T2_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...
                                = sigw_T2_21(i_TRs2, i_TRs1, i_s2, i_s1, i_d1);

                        % Otherwise, we have to do the CRLB
                        else
                            % Inner maximization: store max SDs over (T1,T2,kap) range
                            worst_sigT1  = NaN(length(T1), length(T2), length(kap));
                            worst_sigT2  = NaN(length(T1), length(T2), length(kap));

                            for i_t1 = 1:length(T1)
                            for i_t2 = 1:length(T2)
                            for i_kap = 1:length(kap)
                                [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                                    worst_sigT2(i_t1, i_t2, i_kap)]...
                                    = norm_crlb_dess_3parm(T1(i_t1), T2(i_t2), wf,...
                                    kap(i_kap)*spgr_21, kap(i_kap)*dess_21,...
                                    TRs_21, TRd_21, TE, Sigma_inv, time_comp);
                            
                            end
                            end
                            end

                            % Store the worst-case sig_T1 and sig_T2 values
                            sigw_T1_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1) = max(worst_sigT1(:));
                            sigw_T2_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1) = max(worst_sigT2(:));

                            % Store the worst-case Psi and its corresponding (T1, T2) index
                            worst_Psi = a*worst_sigT1 + b*worst_sigT2;
                            [Psi_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1), idx_Psiw]...
                                = max(worst_Psi(:));

                            % Extract this index and use it to save corresponding sigmas
                            [T1_idx, T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw);
                            sig_T1_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...  
                                = worst_sigT1(T1_idx, T2_idx, kap_idx);
                            sig_T2_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1)...  
                                = worst_sigT2(T1_idx, T2_idx, kap_idx);
                        end
                    end
                end
                end
            end
        end
        end

        % Find the indices of (multiple) Psi_21 minima, to within a tolerance
        Psimin_idx_21 = find( ((Psi_21 - min(Psi_21(:))) ./ min(Psi_21(:))) <= tol);
        num_min_21 = length(Psimin_idx_21);

        %% Step Two: Select one of these minima based on robustness over kappa 
        % Indices of each of the "num_min_21" total minima
        TRs1_idx_21 = NaN(num_min_21, 1);
        TRs2_idx_21 = NaN(num_min_21, 1);
        s1_idx_21   = NaN(num_min_21, 1);
        s2_idx_21   = NaN(num_min_21, 1);
        d1_idx_21   = NaN(num_min_21, 1);

        sigT1_worst_21  = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));
        sigT2_worst_21  = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));

        for i_min = 1:num_min_21
            % Convert the 1D index to ND-grid indices
            [TRs1_idx_21(i_min), TRs2_idx_21(i_min), s1_idx_21(i_min), s2_idx_21(i_min),...
                d1_idx_21(i_min)] = ind2sub(size(Psi_21), Psimin_idx_21(i_min));

            % Store the minimizing parameters
            TRs1min_21 = TR1_spgr(TRs1_idx_21(i_min));
            TRs2min_21 = TR2_spgr(TRs2_idx_21(i_min));
            TRd1min_21 = TR_tot - (TRs1min_21 + TRs2min_21);
            s1min_21   = spgr1(s1_idx_21(i_min));
            s2min_21   = spgr2(s2_idx_21(i_min));
            d1min_21   = dess1(d1_idx_21(i_min));

            % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
            for i_t1 = 1:length(T1_rob)
            for i_t2 = 1:length(T2_rob)
            for i_kap = 1:length(kap_rob)
                [~, ~, sigT1_worst_21(i_min, i_t1, i_t2, i_kap),...
                    sigT2_worst_21(i_min, i_t1, i_t2, i_kap)] = ...
                    norm_crlb_dess_3parm(T1_rob(i_t1), T2_rob(i_t2), wf,...
                    kap_rob(i_kap)*[s1min_21 s2min_21]',...
                    kap_rob(i_kap)*[d1min_21]', [TRs1min_21 TRs2min_21]',...
                    [TRd1min_21]', TE, Sigma_inv, time_comp);
            end
            end
            end
        end

        % Store the worst-case sigT1 and sigT2
        sigT1_min_21 = max(reshape(sigT1_worst_21,...
            [num_min_21 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        sigT2_min_21 = max(reshape(sigT2_worst_21,...
            [num_min_21 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
        
        % Compute the worst-case Psi_21 value over the wider range
        Psi_min_21 = a*sigT1_worst_21 + b*sigT2_worst_21; 
        [Psi_maxmin_21, wide_idx_21] = max(reshape(Psi_min_21,...
            [num_min_21 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

        % Extract the parameters that minimize Psi_maxmin_21 over local minima
        % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
        [Psi_star_21, i_star_21] = min(Psi_maxmin_21);
        sigT1_star_21 = sigT1_min_21(i_star_21);
        sigT2_star_21 = sigT2_min_21(i_star_21);
        
        TRs1star_21 = TRs1_idx_21(i_star_21);
        TRs2star_21 = TRs2_idx_21(i_star_21);
        s1star_21   = s1_idx_21(i_star_21);
        s2star_21   = s2_idx_21(i_star_21);
        d1star_21   = d1_idx_21(i_star_21);

        % Display results
        fprintf('Selected scan design: ');
        fprintf('%1.1f ms, %2.0f deg); (%1.1f ms, %2.0f deg); (%1.1f ms, %2.0f deg)\n',...
            TR1_spgr(TRs1star_21), spgr1(s1star_21)*180/pi,...
            TR2_spgr(TRs2star_21), spgr2(s2star_21)*180/pi,...
            TR_tot - sum([TR1_spgr(TRs1star_21) TR2_spgr(TRs2star_21)]),...
            dess1(d1star_21)*180/pi);
        fprintf('\t Worst-case sigT1: %6.2f\n',...
            sigw_T1_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
        fprintf('\t Robust-range sigT1: %6.2f\n', sigT1_star_21);
        fprintf('\t Worst-case sigT2: %6.3f\n',...
            sigw_T2_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
        fprintf('\t Robust-range sigT2: %6.2f\n', sigT2_star_21);
        fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
            a, sig_T1_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21),...
            b, sig_T2_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21),...
            Psi_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
        fprintf('\t Robust-range Psi: %6.3f\n', Psi_star_21);
        fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
            num_min_21, 100*tol, toc);

        % Save useful outputs 
        suffix = sprintf('_%2.2fms_%1.1f,%d', TR_tot, a, b);
        suffix = strrep(suffix, ',', '_');
        suffix = strrep(suffix, '.', 'p');
        if (savedat)
            cd(strcat('TRtot', suffix)); 
            mkdir('(2,1) Scan Profile'); cd('(2,1) Scan Profile');
            label_21 = strcat('2SPGR,1DESS', suffix, '_3parm', '.mat');
            save(label_21, 'sig_T1_21', 'sig_T2_21', 'sigw_T1_21', 'sigw_T2_21',...
                'Psi_21', 'a', 'b', 'Sigma_inv', 'TR1_spgr', 'TR2_spgr',...
                'spgr1', 'spgr2', 'dess1', 'T1', 'T2');
        end

        %% Graphical Output: 2D Heatmaps for 2 SPGR, 1 DESS
        % Minimizer locations(s)
        if (TR1_spgr(TRs1star_21) == TR2_spgr(TRs2star_21))
            xmin = [spgr1(s1star_21) spgr2(s2star_21)] * 180/pi;
            ymin = [spgr2(s2star_21) spgr1(s1star_21)] * 180/pi;
            zmin = [dess1(d1star_21) dess1(d1star_21)] * 180/pi;
        else
            xmin = [spgr1(s1star_21)] * 180/pi;
            ymin = [spgr2(s2star_21)] * 180/pi;
            zmin = [dess1(d1star_21)] * 180/pi;
        end
        
        % T1 sigma vs. (spgr1/spgr2)
        figure; hold on; imagesc(180*spgr1/pi, 180*spgr2/pi, ...
            squeeze(sigw_T1_21(TRs1star_21, TRs2star_21, :, :, d1star_21))', T1range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(spgr2*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_vs_s1s2_21', suffix, '.eps')), end;
        
        % T2 sigma vs. (spgr1/spgr2);
        figure; hold on; imagesc(180*spgr1/pi, 180*spgr2/pi, ...
            squeeze(sigw_T2_21(TRs1star_21, TRs2star_21, :, :, d1star_21))', T2range);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(spgr2*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_vs_s1s2_21', suffix, '.eps')), end;
        
        % Worst-case Psi vs. (spgr1/spgr2);
        figure; hold on; imagesc(180*spgr1/pi, 180*spgr2/pi, ...
            squeeze(Psi_21(TRs1star_21, TRs2star_21, :, :, d1star_21))', Psirange);...
            scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(spgr2*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_vs_s1s2_21', suffix, '.eps')), end;
        
            
        % T1 sigma vs. (spgr1/dess1)
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(sigw_T1_21(TRs1star_21, TRs2star_21, :, s2star_21, :))', T1range);...
            scatter(xmin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_vs_s1d1_21', suffix, '.eps')), end;
        
        % T2 sigma vs. (spgr1/dess1);
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(sigw_T2_21(TRs1star_21, TRs2star_21, :, s2star_21, :))', T2range);...
            scatter(xmin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_vs_s1d1_21', suffix, '.eps')), end;
        
        % Worst-case Psi vs. (spgr1/dess1);
        figure; hold on; imagesc(180*spgr1/pi, 180*dess1/pi, ...
            squeeze(Psi_21(TRs1star_21, TRs2star_21, :, s2star_21, :))', Psirange);...
            scatter(xmin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr1*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_vs_s1d1_21', suffix, '.eps')), end;
        
        
        % T1 sigma vs. (spgr2/dess1)
        figure; hold on; imagesc(180*spgr2/pi, 180*dess1/pi, ...
            squeeze(sigw_T1_21(TRs1star_21, TRs2star_21, s1star_21, :, :))', T1range);...
            scatter(ymin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T1 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr2*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t1_vs_s2d1_21', suffix, '.eps')), end;
        
        % T2 sigma vs. (spgr2/dess1);
        figure; hold on; imagesc(180*spgr2/pi, 180*dess1/pi, ...
            squeeze(sigw_T2_21(TRs1star_21, TRs2star_21, s1star_21, :, :))', T2range);...
            scatter(ymin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr2*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('sigw_t2_vs_s2d1_21', suffix, '.eps')), end;
        
        % Worst-case Psi vs. (spgr2/dess1);
        figure; hold on; imagesc(180*spgr2/pi, 180*dess1/pi, ...
            squeeze(Psi_21(TRs1star_21, TRs2star_21, s1star_21, :, :))', Psirange);...
            scatter(ymin(1), zmin(1), 250, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
            if(~pr), xlabel('SPGR flip2 (deg)', 'FontSize', 14), end;...
            if(~pr), ylabel('DESS flip1 (deg)', 'FontSize', 14), end;...
            if(~pr), title('Worst-case Psi (ms)', 'FontSize', 14), end;...
            axis xy square; axis([minmax(spgr2*180/pi)', minmax(dess1*180/pi)']);...
            colormap('hot'); colorbar;...
            if (pr) print('-depsc', strcat('psiw_vs_s2d1_21', suffix, '.eps')), end;
            if (savedat), cd ../.., end
    otherwise 
        warning('Unexpected scan profile!');
end