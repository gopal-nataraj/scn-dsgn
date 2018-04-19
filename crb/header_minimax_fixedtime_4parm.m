% Header File -- Min-max CRLB analysis, with TR/flip angle variation
% Worst-case M0s/T1/T2/kappa CRLB over T1/T2 ROI, varying TRs and flips
% Instead of TR-compensation, use fixed total imaging time
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

% Set imaging constraints
TR_tot = 17.28*2;                   % ms
TRs_min = 8;                        % ms
TRd_min = 17.28;                    % ms
TRs_res = 1;                        % TR SPGR resolution
TRd_res = 1;                        % TR DESS resolution

% Implicit parameters, including a T1/T2/kappa range of interest
% Inner minimax problem: tight bounds
T1 = linspace(800, 1400, 5);        % WM/GM voxels of interest
T2 = linspace(50, 120, 5);          % WM/GM voxels of interest
% kap = 2 .^ linspace(-0.5, 0.5, 5);  % a.u. (factor-of-sqrt(2) variation)
kap = linspace(0.9, 1.1, 5);        % a.u. (10% variation)
wf = 0;                             % rad/ms

% Outer robustness criterion: check minima over loose bounds
T1_rob  = linspace(400, 2000, 9);   % ms 
T2_rob  = linspace(40, 200, 9);     % ms
kap_rob = 2 .^ linspace(-1, 1, 5);  % a.u. (factor-of-two variation)

% Other constant declarations
TE = 4.67;                          % ms
noise_var_1coil = 2.62e-07;         % a.u.
noise_var_ssos = noise_var_1coil/2; % High-SNR regime approximation

a = 0;                            % T1 relative importance parameter
b = 1;                              % T2 relative importance parameter
c = 0;                              % kap relative importance parameter

tol = 0.01;                         % Global minima tolerance
time_comp = 0;                      % Toggle time compensation on/off
pr = 1;                             % Toggle print on/off
savedat = 0;                        % Toggle saving .mat files on/off
disp_scale = [-1 4];                % Logarithmic display scaling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence One: 0 SPGR, 2 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controllable parameter grid: varying TRs and flip angles
TRd_max = TR_tot - TRd_min;
TR1_dess = [TRd_min : TRd_res : TRd_max];
dess1 = linspace(5,90,18)*pi/180;   % rad
dess2 = linspace(5,90,18)*pi/180;   % rad

% Introduce the experiment
fprintf('\n\nScan Profile: (0 SPGR, 2 DESS)\n');
fprintf('Total TR constraint: %d ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f,%1.1f,%1.1f)\n', 0, a, b, c);

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
        
        % Inner maximization: store max SDs over (T1,T2,kap) range
        worst_sigT1  = NaN(length(T1), length(T2), length(kap));
        worst_sigT2  = NaN(length(T1), length(T2), length(kap));
        worst_sigkap = NaN(length(T1), length(T2), length(kap));
        
        for i_t1 = 1:length(T1)
        for i_t2 = 1:length(T2)
        for i_kap = 1:length(kap)
            [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                worst_sigT2(i_t1, i_t2, i_kap),...
                worst_sigkap(i_t1, i_t2, i_kap)]...
                = norm_crlb_dess_4parm(T1(i_t1), T2(i_t2), kap(i_kap),...
                wf, spgr_02, dess_02, TRs_02, TRd_02, TE, Sigma_inv, time_comp);
        end
        end
        end
        
        % Store the worst-case sig_T1 and sig_T2 values
        sigw_T1_02(i_TRd1, i_d1, i_d2) = max(worst_sigT1(:));
        sigw_T2_02(i_TRd1, i_d1, i_d2) = max(worst_sigT2(:));
        
        % Store the worst-case Psi and its corresponding (T1, T2, kap) index
        worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
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

%% Step Two: Select one of these minima based on robustness over kappa
% Indices of each of the "num_min_02" total minima
TRd1_idx_02 = NaN(num_min_02, 1);
d1_idx_02  = NaN(num_min_02, 1);
d2_idx_02  = NaN(num_min_02, 1);

sigT1_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_02 = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));

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
            sigT2_worst_02(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_02(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, []', [d1min_02 d2min_02]', []',...
            [TRd1min_02 TRd2min_02]', TE, Sigma_inv, 0);
    end
    end
    end
end

% Compute the worst-case Psi_02 value over the wider range
Psi_min_02 = a*sigT1_worst_02 + b*sigT2_worst_02 + c*sigkap_worst_02;
[Psi_maxmin_02, wide_idx_02] = max(reshape(Psi_min_02,...
    [num_min_02 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameter indices that minimize Psi_maxmin_02 over minima
% Psi_star - min(Psi(:)) is the degradation from narrow->wide range
[Psi_star_02, i_star_02] = min(Psi_maxmin_02);
TRd1star_02 = TRd1_idx_02(i_star_02);
d1star_02   = d1_idx_02(i_star_02);
d2star_02   = d2_idx_02(i_star_02);

% Display results 
fprintf('Selected scan design: (%d ms, %2.0f deg); (%d ms, %2.0f deg)\n',...
    TR1_dess(TRd1star_02), dess1(d1star_02)*180/pi,...
    TR_tot-TR1_dess(TRd1star_02), dess2(d2star_02)*180/pi);
fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
    a, sig_T1_02(TRd1star_02, d1star_02, d2star_02),...
    b, sig_T2_02(TRd1star_02, d1star_02, d2star_02),...
    Psi_02(TRd1star_02, d1star_02, d2star_02));
fprintf('\t Worst-case sigT1: %6.2f\n',...
    sigw_T1_02(TRd1star_02, d1star_02, d2star_02));
fprintf('\t Worst-case sigT2: %6.3f\n',...
    sigw_T2_02(TRd1star_02, d1star_02, d2star_02));
fprintf('\t Worst-case robust Psi: %6.3f\n', Psi_star_02);
fprintf('\t Design selected from %d candidates, at a %d%% tolerance\n\n',...
    num_min_02, 100*tol);

% Save useful outputs
suffix = sprintf('_%2.2fms_%1.1f,%d,%d', TR_tot, a, b, c);
cd(strcat('TRtot', suffix)); mkdir('(0,2) Scan Profile'); cd('(0,2) Scan Profile');
if (savedat)
    label_02 = strcat('0SPGR,2DESS', suffix, '.mat');
    save(label_02, 'sig_T1_02', 'sig_T2_02', 'sigw_T1_02', 'sigw_T2_02',...
        'Psi_02', 'a', 'b', 'c', 'Sigma_inv','TR1_dess', 'dess1',...
        'dess2', 'T1', 'T2', 'kap');
end

%% Graphical Output: 2D Heatmaps for 0 SPGR, 2 DESS
% T1 sigma vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sigw_T1_02(TRd1star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print('-depsc', strcat('sigw_t1_vs_d1d2_02', suffix, '.eps')), end;

% T2 sigma vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sigw_T2_02(TRd1star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print('-depsc', strcat('sigw_t2_vs_d1d2_02', suffix, '.eps')), end;

% T2 sigma vs. (dess1/dess2)
xmin = [dess1(d1star_02)*180/pi dess2(d2star_02)*180/pi];
ymin = [dess2(d2star_02)*180/pi dess1(d1star_02)*180/pi];
figure; hold on; imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(sigw_T2_02(TRd1star_02, :, :))', [0 20]);...
    scatter(xmin, ymin, 150, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');
    xlabel('DESS flip1 (deg)', 'FontSize', 14);...
    ylabel('DESS flip2 (deg)', 'FontSize', 14);...
    title('Worst-case T2 Standard Deviation (ms)', 'FontSize', 14);...
    axis xy square; axis([minmax(dess1*180/pi)' minmax(dess2*180/pi)']);
    colormap('hot'); colorbar;...
    if (pr) print('-depsc', strcat('sigw_t2_vs_d1d2_02', suffix, '.eps')), end;
    
% Psi vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(Psi_02(TRd1star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print('-depsc', strcat('psi_vs_d1d2_02', suffix, '.eps')), end; cd ../..;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Two: 2 SPGR, 1 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controllable parameter grid: varying TRs and flip angles
TRs_max = TR_tot - TRs_min - TRd_min;
TR1_spgr = [TRs_min : TRs_res : TRs_max];
TR2_spgr = [TRs_min : TRs_res : TRs_max];
spgr1 = linspace(5,90,18)*pi/180;   % rad
spgr2 = linspace(5,90,18)*pi/180;   % rad
dess1 = linspace(5,90,18)*pi/180;   % rad

% Introduce the experiment
fprintf('\n\nScan Profile: (2 SPGR, 1 DESS)\n');
fprintf('Total TR constraints: %d ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f,%1.1f,%1.1f)\n', 0, a, b, c);

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
                    worst_sigkap = NaN(length(T1), length(T2), length(kap));

                    for i_t1 = 1:length(T1)
                    for i_t2 = 1:length(T2)
                    for i_kap = 1:length(kap)
                        [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                            worst_sigT2(i_t1, i_t2, i_kap),...
                            worst_sigkap(i_t1, i_t2, i_kap)]...
                            = norm_crlb_dess_4parm(T1(i_t1), T2(i_t2),...
                            kap(i_kap), wf, spgr_21, dess_21, TRs_21,...
                            TRd_21, TE, Sigma_inv, time_comp);
                    end
                    end
                    end
                    
                    % Store the worst-case sig_T1 and sig_T2 values
                    sigw_T1_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1) = max(worst_sigT1(:));
                    sigw_T2_21(i_TRs1, i_TRs2, i_s1, i_s2, i_d1) = max(worst_sigT2(:));

                    % Store the worst-case Psi and its corresponding (T1, T2, kap) index
                    worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap;
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
sigkap_worst_21 = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));

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
    for i_kap = 1:length(kap)
        [~, ~, sigT1_worst_21(i_min, i_t1, i_t2, i_kap),...
            sigT2_worst_21(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_21(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, [s1min_21 s2min_21]', [d1min_21]',...
            [TRs1min_21 TRs2min_21]', [TRd1min_21]', TE, Sigma_inv, 0);
    end
    end
    end
end

% Compute the worst-case Psi_21 value over the wider range
Psi_min_21 = a*sigT1_worst_21 + b*sigT2_worst_21 + c*sigkap_worst_21; 
[Psi_maxmin_21, wide_idx_21] = max(reshape(Psi_min_21,...
    [num_min_21 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameters that minimize Psi_maxmin_21 over local minima
% Psi_star - min(Psi(:)) is the degradation from kappa variation
[Psi_star_21, i_star_21] = min(Psi_maxmin_21); 
TRs1star_21 = TRs1_idx_21(i_star_21);
TRs2star_21 = TRs2_idx_21(i_star_21);
s1star_21   = s1_idx_21(i_star_21);
s2star_21   = s2_idx_21(i_star_21);
d1star_21   = d1_idx_21(i_star_21);

% Display results
fprintf('Selected scan design: ');
fprintf('%d ms, %2.0f deg); (%d ms, %2.0f deg); (%d ms, %2.0f deg)\n',...
    TR1_spgr(TRs1star_21), spgr1(s1star_21)*180/pi,...
    TR2_spgr(TRs2star_21), spgr2(s2star_21)*180/pi,...
    TR_tot - sum([TR1_spgr(TRs1star_21) TR2_spgr(TRs2star_21)]),...
    dess1(d1star_21)*180/pi);
fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
    a, sig_T1_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21),...
    b, sig_T2_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21),...
    Psi_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
fprintf('\t Worst-case sigT1: %6.2f\n',...
    sigw_T1_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
fprintf('\t Worst-case sigT2: %6.3f\n',...
    sigw_T2_21(TRs1star_21, TRs2star_21, s1star_21, s2star_21, d1star_21));
fprintf('\t Worst-case robust Psi: %6.3f\n', Psi_star_21);
fprintf('\t Design selected from %d candidates, at a %d%% tolerance\n\n',...
    num_min_21, 100*tol);

% Save useful outputs 
suffix = sprintf('_%2.2fms_%1.1f,%d,%d', TR_tot, a, b, c);
cd(strcat('TRtot', suffix)); mkdir('(2,1) Scan Profile'); cd('(2,1) Scan Profile');
if (savedat)
    label_21 = strcat('2SPGR,1DESS', suffix, '.mat');
    save(label_21, 'sig_T1_21', 'sig_T2_21', 'sigw_T1_21', 'sigw_T2_21',...
        'Psi_21', 'a', 'b', 'c', 'Sigma_inv', 'TR1_spgr', 'TR2_spgr',...
        'spgr1', 'spgr2', 'dess1', 'T1', 'T2', 'kap');
end

%% Graphical Output: 2D Heatmaps for 2 SPGR, 1 DESS
% T1 sigma vs. (spgr1/spgr2)
figure; imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    sigw_T1_21(TRs1star_21, TRs2star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_s1s2_21', suffix, '.eps')), end;

% T2 sigma vs. (spgr1/spgr2)
figure; imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    sigw_T2_21(TRs1star_21, TRs2star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_s1s2_21', suffix, '.eps')), end;

% Psi vs. (spgr1/spgr2)
figure; imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_s1s2_21', suffix, '.eps')), end;
    

% T1 sigma vs. (spgr1/dess1)
figure; imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    sigw_T1_21(TRs1star_21, TRs2star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_s1d1_21', suffix, '.eps')), end;

% T2 sigma vs. (spgr1/dess1)
figure; imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    sigw_T2_21(TRs1star_21, TRs2star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_s1d1_21', suffix, '.eps')), end;

% Psi vs. (spgr1/dess1)
figure; imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_s1d1_21', suffix, '.eps')), end;


% T1 sigma vs. (spgr2/dess1)
figure; imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    sigw_T1_21(TRs1star_21, TRs2star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_s2d1_21', suffix, '.eps')), end;

% T2 sigma vs. (spgr2/dess1)
figure; imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    sigw_T2_21(TRs1star_21, TRs2star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_s2d1_21', suffix, '.eps')), end;

% Psi vs. (spgr2/dess1)
figure; imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_s2d1_21', suffix, '.eps')), end; cd ../..;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Three: 0 SPGR, 3 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Controllable parameter grid: varying TRs and flip angles
TRd_max = TR_tot - 2*TRd_min;               
TR1_dess = [TRd_min : TRd_res : TRd_max];
TR2_dess = [TRd_min : TRd_res : TRd_max];
TR3_dess = [TRd_min : TRd_res : TRd_max];   % Must be same as TR1, TR2
dess1 = linspace(5,90,18)*pi/180;   % rad
dess2 = linspace(5,90,18)*pi/180;   % rad
dess3 = linspace(5,90,18)*pi/180;   % rad

% Introduce the experiment
fprintf('\n\nScan Profile: (0 SPGR, 3 DESS)\n');
fprintf('Total TR constraints: %d ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f,%1.1f,%1.1f)\n', 0, a, b, c);

% Grid definitions 
sig_T1_03   = inf(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2), length(dess3));
sig_T2_03   = inf(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2), length(dess3));
sigw_T1_03  = inf(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2), length(dess3));
sigw_T2_03  = inf(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2), length(dess3));
Psi_03      = inf(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2), length(dess3));

% Covariance matrix: scaled identity with M = 6
Sigma_inv = (1 ./ noise_var_ssos) * speye(6);

%% Step One: Min-max CRLB analysis over parameter grid
spgr_03 = []';
TRs_03 = []';
min_perm3 = @(var, t1, t2, t3, f1, f2, f3) min([...
    var(t1, t2, f1, f2, f3), var(t1, t3, f1, f3, f2),...
    var(t2, t1, f2, f1, f3), var(t2, t3, f2, f3, f1),...
    var(t3, t1, f3, f1, f2), var(t3, t2, f3, f2, f1)]);
for i_d1 = 1:length(dess1)
for i_d2 = 1:length(dess2)
for i_d3 = 1:length(dess3)
    dess_03 = [dess1(i_d1) dess2(i_d2) dess3(i_d3)]';
    
    for i_TRd1 = 1:length(TR1_dess)
    for i_TRd2 = 1:length(TR2_dess)
        % Locate the "fake" third index from constraint
        TRd3 = TR_tot - sum([TR1_dess(i_TRd1) TR2_dess(i_TRd2)]);
        i_TRd3 = find(TR3_dess == TRd3, 1); 
        TRd_03 = [TR1_dess(i_TRd1) TR2_dess(i_TRd2) TRd3]';
        
        % Only consider points for which TRd3 >= TRd_min
        % Otherwise, do nothing, so Psi_03 remains infinite
        if (min(TRd_03) >= TRd_min)
            % Exploit redundancy: check all permutations of DESS scans to
            % see if we have already computed this combination. 
            if (isfinite(min_perm3(Psi_03, i_TRd1, i_TRd2, i_TRd3,...
                    i_d1, i_d2, i_d3)))
                Psi_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = min_perm3(...
                    Psi_03, i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3);
                sig_T1_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = min_perm3(...
                    sig_T1_03, i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3);
                sig_T2_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = min_perm3(...
                    sig_T2_03, i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3);
                sigw_T1_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = min_perm3(...
                    sigw_T1_03, i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3);
                sigw_T2_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = min_perm3(...
                    sigw_T2_03, i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3);
                
            % Otherwise, we have to do the CRLB
            else
                % Inner maximization: store max SDs over (T1,T2,kap) range
                worst_sigT1  = NaN(length(T1), length(T2), length(kap));
                worst_sigT2  = NaN(length(T1), length(T2), length(kap));
                worst_sigkap = NaN(length(T1), length(T2), length(kap));

                for i_t1 = 1:length(T1)
                for i_t2 = 1:length(T2)
                for i_kap = 1:length(kap)
                    [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                        worst_sigT2(i_t1, i_t2, i_kap),...
                        worst_sigkap(i_t1, i_t2, i_kap)]...
                        = norm_crlb_dess_4parm(T1(i_t1), T2(i_t2),...
                        kap(i_kap), wf, spgr_03, dess_03, TRs_03, TRd_03,...
                        TE, Sigma_inv, time_comp);
                end
                end
                end

                % Store the worst-case sig_T1 and sig_T2 values
                sigw_T1_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = max(worst_sigT1(:));
                sigw_T2_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3) = max(worst_sigT2(:));
                
                % Store the worst-case Psi and its corresponding (T1, T2, kap) index
                worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
                [Psi_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3), idx_Psiw]...
                        = max(worst_Psi(:));

                % Extract this index and use it to save corresponding sigmas
                [T1_idx, T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw);
                sig_T1_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3)...
                    = worst_sigT1(T1_idx, T2_idx, kap_idx);
                sig_T2_03(i_TRd1, i_TRd2, i_d1, i_d2, i_d3)...
                    = worst_sigT2(T1_idx, T2_idx, kap_idx);
            end
        end
    end
    end
end
end
end

% Find the indices of (multiple) Psi_03 minima, to within a tolerance
Psimin_idx_03 = find( ((Psi_03 - min(Psi_03(:))) ./ min(Psi_03(:))) <= tol);
num_min_03 = length(Psimin_idx_03);

%% Step Two: Select one of these minima based on robustness over kappa 
% Indices of each of the "num_min_03" total minima
TRd1_idx_03 = NaN(num_min_03, 1);
TRd2_idx_03 = NaN(num_min_03, 1);
d1_idx_03   = NaN(num_min_03, 1);
d2_idx_03   = NaN(num_min_03, 1);
d3_idx_03   = NaN(num_min_03, 1);

sigT1_worst_03  = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_03  = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_03 = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));

for i_min = 1:num_min_03
    % Convert the 1D index to ND-grid indices
    [TRd1_idx_03(i_min), TRd2_idx_03(i_min), d1_idx_03(i_min), d2_idx_03(i_min),...
        d3_idx_03(i_min)] = ind2sub(size(Psi_03), Psimin_idx_03(i_min));
    
    % Store the minimizing parameters
    TRd1min_03 = TR1_dess(TRd1_idx_03(i_min));
    TRd2min_03 = TR2_dess(TRd2_idx_03(i_min));
    TRd3min_03 = TR_tot - (TRd1min_03 + TRd2min_03);
    d1min_03   = dess1(d1_idx_03(i_min));
    d2min_03   = dess2(d2_idx_03(i_min));
    d3min_03   = dess3(d3_idx_03(i_min));
    
    % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
    % Note now that we have turned on time compensation!
    for i_t1 = 1:length(T1_rob)
    for i_t2 = 1:length(T2_rob)
    for i_kap = 1:length(kap_rob)
        [~, ~, sigT1_worst_03(i_min, i_t1, i_t2, i_kap),...
            sigT2_worst_03(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_03(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, []', [d1min_03 d2min_03 d3min_03]', []',...
            [TRd1min_03 TRd2min_03 TRd3min_03]', TE, Sigma_inv, 1);
    end
    end
    end
end

% Compute the worst-case Psi_tild_03 value over the wider range
Psi_min_03 = a*sigT1_worst_03 + b*sigT2_worst_03 + c*sigkap_worst_03; 
[Psi_maxmin_03, wide_idx_03] = max(reshape(Psi_min_03,...
    [num_min_03 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameter indices that minimize Psi_tild_maxmin_03 over local minima
% Psi_star - min(Psi_tild(:)) is the degradation from narrow->wide range
[Psi_star_03, i_star_03] = min(Psi_maxmin_03); 
TRd1star_03 = TRd1_idx_03(i_star_03);
TRd2star_03 = TRd2_idx_03(i_star_03);
d1star_03   = d1_idx_03(i_star_03); 
d2star_03   = d2_idx_03(i_star_03);
d3star_03   = d3_idx_03(i_star_03);

% Display results 
fprintf('Selected scan design: ');
fprintf('%d ms, %2.0f deg); (%d ms, %2.0f deg); (%d ms, %2.0f deg)\n',...
    TR1_dess(TRd1star_03), dess1(d1star_03)*180/pi,...
    TR2_dess(TRd2star_03), dess2(d2star_03)*180/pi,...
    TR_tot - sum([TR1_dess(TRd1star_03) TR2_dess(TRd2star_03)]),...
    dess3(d3star_03)*180/pi);
fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
    a, sig_T1_03(TRd1star_03, TRd2star_03, d1star_03, d2star_03, d3star_03),...
    b, sig_T2_03(TRd1star_03, TRd2star_03, d1star_03, d2star_03, d3star_03),...
    Psi_03(TRd1star_03, TRd2star_03, d1star_03, d2star_03, d3star_03));
fprintf('\t Worst-case sigT1: %6.2f\n',...
    sigw_T1_03(TRd1star_03, TRd2star_03, d1star_03, d2star_03, d3star_03));
fprintf('\t Worst-case sigT2: %6.3f\n',...
    sigw_T2_03(TRd1star_03, TRd2star_03, d1star_03, d2star_03, d3star_03));
fprintf('\t Worst-case robust Psi: %6.3f\n', Psi_star_03);
fprintf('\t Design selected from %d candidates, at a %d%% tolerance\n\n',...
    num_min_03, 100*tol);

% Save useful outputs
suffix = sprintf('_%2.2fms_%1.1f,%d,%d', TR_tot, a, b, c);
cd(strcat('TRtot', suffix)); mkdir('(0,3) Scan Profile'); cd('(0,3) Scan Profile');
if (savedat)
    label_03 = strcat('0SPGR,3DESS', suffix, '.mat');
    save(label_03, 'sig_T1_03', 'sig_T2_03', 'sigw_T1_03', 'sigw_T2_03',...
        'Psi_03', 'a', 'b', 'c', 'Sigma_inv', 'TR1_dess', 'TR2_dess',...
        'dess1', 'dess2', 'dess3', 'T1', 'T2', 'kap');
end


%% Graphical Output: 2D Heatmaps for 0 SPGR, 3 DESS
% T1 sigma vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(sigw_T1_03(...
    TRd1star_03, TRd2star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_d1d2_03', suffix, '.eps')), end;
    
% T2 sigma vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(sigw_T2_03(...
    TRd1star_03, TRd2star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_d1d2_03', suffix, '.eps')), end;

% Psi vs. (dess1/dess2)
figure; imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_d1d2_03', suffix, '.eps')), end;


% T1 sigma vs. (dess1/dess3)
figure; imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(sigw_T1_03(...
    TRd1star_03, TRd2star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_d1d3_03', suffix, '.eps')), end;
    
% T2 sigma vs. (dess1/dess3)
figure; imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(sigw_T2_03(...
    TRd1star_03, TRd2star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_d1d3_03', suffix, '.eps')), end;

% Psi vs. (dess1/dess3)
figure; imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_d1d3_03', suffix, '.eps')), end;


% T1 sigma vs. (dess2/dess3)
figure; imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(sigw_T1_03(...
    TRd1star_03, TRd2star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t1_vs_d2d3_03', suffix, '.eps')), end;
    
% T2 sigma vs. (dess2/dess3)
figure; imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(sigw_T2_03(...
    TRd1star_03, TRd2star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('sigw_t2_vs_d2d3_03', suffix, '.eps')), end;

% Psi vs. (dess2/dess3)
figure; imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print('-depsc', strcat('psi_vs_d2d3_03', suffix, '.eps')), end; cd ../..;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Four: 0 SPGR, 4 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Controllable parameter grid: varying TRs and flip angles
TRd_max = TR_tot - 3*TRd_min;               
TR1_dess = [TRd_min : TRd_res : TRd_max]; 
TR2_dess = [TRd_min : TRd_res : TRd_max];
TR3_dess = [TRd_min : TRd_res : TRd_max];   % Must be same as TR1, TR2
TR4_dess = [TRd_min : TRd_res : TRd_max];   % Must be same as TR1, TR2
dess1 = linspace(5,90,18)*pi/180;   % rad
dess2 = linspace(5,90,18)*pi/180;   % rad
dess3 = linspace(5,90,18)*pi/180;   % rad
dess4 = linspace(5,90,18)*pi/180;   % rad

% Introduce the experiment
fprintf('\n\nScan Profile: (0 SPGR, 4 DESS)\n');
fprintf('Total TR constraints: %d ms\n', TR_tot);
fprintf('Weightings: (%1.1f,%1.1f,%1.1f,%1.1f)\n', 0, a, b, c);

% Grid definitions 
sig_T1_04   = inf(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3), length(dess4));
sig_T2_04   = inf(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3), length(dess4));
sigw_T1_04  = inf(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3), length(dess4));
sigw_T2_04  = inf(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3), length(dess4));
Psi_04      = inf(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3), length(dess4));

% Covariance matrix: scaled identity with M = 8
Sigma_inv = (1 ./ noise_var_ssos) * speye(8);

%% Step One: Min-max CRLB analysis over parameter grid
spgr_04 = []';
TRs_04 = []';
min_perm3 = @(var, t1, t2, t3, t4, f1, f2, f3, f4) min([...
    var(t1, t2, t3, f1, f2, f3, f4), var(t1, t2, t4, f1, f2, f4, f3),...
    var(t1, t3, t2, f1, f3, f2, f4), var(t1, t3, t4, f1, f3, f4, f2),...
    var(t1, t4, t2, f1, f4, f2, f3), var(t1, t4, t3, f1, f4, f3, f2)]);
min_perm4 = @(var, t1, t2, t3, t4, f1, f2, f3, f4) min([...
    min_perm3(var, t1, t2, t3, t4, f1, f2, f3, f4),...
    min_perm3(var, t2, t1, t3, t4, f2, f1, f3, f4),...
    min_perm3(var, t3, t2, t1, t4, f3, f2, f1, f4),...
    min_perm3(var, t4, t2, t3, t1, f4, f2, f3, f1)]);
for i_d1 = 1:length(dess1)
for i_d2 = 1:length(dess2)
for i_d3 = 1:length(dess3)
for i_d4 = 1:length(dess4)
    dess_04 = [dess1(i_d1) dess2(i_d2) dess3(i_d3) dess4(i_d4)]';
    
    for i_TRd1 = 1:length(TR1_dess)
    for i_TRd2 = 1:length(TR2_dess)
    for i_TRd3 = 1:length(TR3_dess)
        % Locate the "fake" fourth index from constraint
        TRd4 = TR_tot - sum([TR1_dess(i_TRd1) TR2_dess(i_TRd2) TR3_dess(i_TRd3)]);
        i_TRd4 = find(TR4_dess == TRd4, 1); 
        TRd_04 = [TR1_dess(i_TRd1) TR2_dess(i_TRd2) TR3_dess(i_TRd3) TRd4]';
        
        % Only consider points for which TRd4 >= TRd_min
        % Otherwise, do nothing, so Psi_04 remains infinite
        if (min(TRd_04) >= TRd_min)
            % Exploit redundancy: check all permutations of DESS scans to
            % see if we have already computed this combination. 
            if (isfinite(min_perm4(Psi_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4,...
                    i_d1, i_d2, i_d3, i_d4)))
                Psi_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) =...
                    min_perm4(Psi_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4, i_d1, i_d2, i_d3, i_d4);
                sig_T1_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) =...
                    min_perm4(sig_T1_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4, i_d1, i_d2, i_d3, i_d4);
                sig_T2_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) =...
                    min_perm4(sig_T2_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4, i_d1, i_d2, i_d3, i_d4);
                sigw_T1_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) =...
                    min_perm4(sigw_T1_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4, i_d1, i_d2, i_d3, i_d4);
                sigw_T2_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) =...
                    min_perm4(sigw_T2_04, i_TRd1, i_TRd2, i_TRd3, i_TRd4, i_d1, i_d2, i_d3, i_d4);
                
            % Otherwise, we have to do the CRLB
            else
                % Inner maximization: store max SDs over (T1,T2,kap) range
                worst_sigT1  = NaN(length(T1), length(T2), length(kap));
                worst_sigT2  = NaN(length(T1), length(T2), length(kap));
                worst_sigkap = NaN(length(T1), length(T2), length(kap));

                for i_t1 = 1:length(T1)
                for i_t2 = 1:length(T2)
                for i_kap = 1:length(kap)
                    [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
                        worst_sigT2(i_t1, i_t2, i_kap),...
                        worst_sigkap(i_t1, i_t2, i_kap)]...
                        = norm_crlb_dess_4parm(T1(i_t1), T2(i_t2),...
                        kap(i_kap), wf, spgr_04, dess_04, TRs_04, TRd_04,...
                        TE, Sigma_inv, time_comp);
                end
                end
                end
                
                % Store the worst-case sig_T1 and sig_T2 values
                sigw_T1_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) = max(worst_sigT1(:));
                sigw_T2_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4) = max(worst_sigT2(:));
                
                % Store the worst-case Psi and its corresponding (T1, T2, kap) index
                worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
                [Psi_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4),...
                    idx_Psiw] = max(worst_Psi(:));

                % Extract this index and use it to save corresponding sigmas
                [T1_idx, T2_idx, kap_idx] = ind2sub(size(worst_Psi), idx_Psiw);
                sig_T1_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4)...
                    = worst_sigT1(T1_idx, T2_idx, kap_idx);
                sig_T2_04(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3, i_d4)...
                    = worst_sigT2(T1_idx, T2_idx, kap_idx);
            end
        end
    end
    end
    end
end
end
end
end

% Find the indices of (multiple) Psi_04 minima, to within a tolerance
Psimin_idx_04 = find( ((Psi_04 - min(Psi_04(:))) ./ min(Psi_04(:))) <= tol);
num_min_04 = length(Psimin_idx_04);

%% Step Two: Select one of these minima based on robustness over kappa 
% Indices of each of the "num_min_04" total minima
TRd1_idx_04 = NaN(num_min_04, 1);
TRd2_idx_04 = NaN(num_min_04, 1);
TRd3_idx_04 = NaN(num_min_04, 1);
d1_idx_04   = NaN(num_min_04, 1);
d2_idx_04   = NaN(num_min_04, 1);
d3_idx_04   = NaN(num_min_04, 1);
d4_idx_04   = NaN(num_min_04, 1);

sigT1_worst_04  = NaN(num_min_04, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_04  = NaN(num_min_04, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_04 = NaN(num_min_04, length(T1_rob), length(T2_rob), length(kap_rob));

for i_min = 1:num_min_04
    % Convert the 1D index to ND-grid indices
    [TRd1_idx_04(i_min), TRd2_idx_04(i_min), TRd3_idx_04(i_min),...
        d1_idx_04(i_min), d2_idx_04(i_min),...
        d3_idx_04(i_min), d4_idx_04(i_min)]...
        = ind2sub(size(Psi_04), Psimin_idx_04(i_min));
    
    % Store the minimizing parameters
    TRd1min_04 = TR1_dess(TRd1_idx_04(i_min));
    TRd2min_04 = TR2_dess(TRd2_idx_04(i_min));
    TRd3min_04 = TR3_dess(TRd3_idx_04(i_min));
    TRd4min_04 = TR_tot - (TRd1min_04 + TRd2min_04 + TRd3min_04);
    d1min_04   = dess1(d1_idx_04(i_min));
    d2min_04   = dess2(d2_idx_04(i_min));
    d3min_04   = dess3(d3_idx_04(i_min));
    d4min_04   = dess4(d4_idx_04(i_min));
    
    % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
    % Note now that we have turned on time compensation!
    for i_t1 = 1:length(T1_rob)
    for i_t2 = 1:length(T2_rob)
    for i_kap = 1:length(kap_rob)
        [~, ~, sigT1_worst_04(i_min, i_t1, i_t2, i_kap),...
            sigT2_worst_04(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_04(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, []', [d1min_04 d2min_04 d3min_04 d4min_04]', []',...
            [TRd1min_04 TRd2min_04 TRd3min_04 TRd4min_04]', TE, Sigma_inv, 1);
    end
    end
    end
end

% Compute the worst-case Psi_tild_04 value over the wider range
Psi_min_04 = a*sigT1_worst_04 + b*sigT2_worst_04 + c*sigkap_worst_04; 
[Psi_maxmin_04, wide_idx_04] = max(reshape(Psi_min_04,...
    [num_min_04 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameter indices that minimize Psi_tild_maxmin_04 over local minima
% Psi_star - min(Psi_tild(:)) is the degradation from narrow->wide range
[Psi_star_04, i_star_04] = min(Psi_maxmin_04); 
TRd1star_04 = TRd1_idx_04(i_star_04);
TRd2star_04 = TRd2_idx_04(i_star_04);
TRd3star_04 = TRd3_idx_04(i_star_04);
d1star_04   = d1_idx_04(i_star_04); 
d2star_04   = d2_idx_04(i_star_04);
d3star_04   = d3_idx_04(i_star_04);
d4star_04   = d4_idx_04(i_star_04);

% Display results 
fprintf('Selected scan design: ');
fprintf('%d ms, %2.0f deg); (%d ms, %2.0f deg); (%d ms, %2.0f deg); (%d ms, %2.0f deg)\n',...
    TR1_dess(TRd1star_04), dess1(d1star_04)*180/pi,...
    TR2_dess(TRd2star_04), dess2(d2star_04)*180/pi,...
    TR3_dess(TRd3star_04), dess3(d3star_04)*180/pi,...
    TR_tot - sum([TR1_dess(TRd1star_04) TR2_dess(TRd2star_04) TR3_dess(TRd3star_04)]),...
    dess4(d4star_04)*180/pi);
fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
    a, sig_T1_04(TRd1star_04, TRd2star_04, TRd3star_04, d1star_04, d2star_04, d3star_04, d4star_04),...
    b, sig_T2_04(TRd1star_04, TRd2star_04, TRd3star_04, d1star_04, d2star_04, d3star_04, d4star_04),...
    Psi_04(TRd1star_04, TRd2star_04, TRd3star_04, d1star_04, d2star_04, d3star_04, d4star_04));
fprintf('\t Worst-case sigT1: %6.2f\n',...
    sigw_T1_04(TRd1star_04, TRd2star_04, TRd3star_04, d1star_04, d2star_04, d3star_04, d4star_04));
fprintf('\t Worst-case sigT2: %6.3f\n',...
    sigw_T2_04(TRd1star_04, TRd2star_04, TRd3star_04, d1star_04, d2star_04, d3star_04, d4star_04));
fprintf('\t Worst-case robust Psi: %6.3f\n', Psi_star_04);
fprintf('\t Design selected from %d candidates, at a %d%% tolerance\n\n',...
    num_min_04, 100*tol);

% Save useful outputs
suffix = sprintf('_%2.2fms_%1.1f,%d,%d', TR_tot, a, b, c);
cd(strcat('TRtot', suffix)); mkdir('(0,4) Scan Profile'); cd('(0,4) Scan Profile');
if (savedat)
    label_04 = strcat('0SPGR,4DESS', suffix, '.mat');
    save(label_04, 'sig_T1_04', 'sig_T2_04', 'sigw_T1_04', 'sigw_T2_04',...
        'Psi_04', 'a', 'b', 'c', 'Sigma_inv', 'TR1_dess', 'TR2_dess',...
        'TR3_dess', 'dess1', 'dess2', 'dess3', 'dess4', 'T1', 'T2', 'kap');
end
cd ../..;