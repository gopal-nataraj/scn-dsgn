% Header File -- Min-max CRLB analysis, with TR/flip angle variation
% Worst-case M0s/T1/T2/kappa CRLB over T1/T2 ROI, varying TRs and flips 
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

% Controllable parameter grid: varying TRs and flip angles
TR1_spgr = linspace(6,16,6);        % ms
TR2_spgr = linspace(6,16,6);        % ms
TR1_dess = linspace(10,20,6);       % ms
TR2_dess = linspace(10,20,6);       % ms
TR3_dess = linspace(10,20,6);       % ms
spgr1 = linspace(5,85,9) * pi/180;  % rad
spgr2 = linspace(5,85,9) * pi/180;  % rad
dess1 = linspace(5,85,9) * pi/180;  % rad
dess2 = linspace(5,85,9) * pi/180;  % rad
dess3 = linspace(5,85,9) * pi/180;  % rad

% Implicit parameters, including a T1/T2/kappa range of interest
% Inner minimax problem: tight bounds
T1  = linspace(500, 900, 5);        % ms (WM 10% range)
T2  = linspace(50, 90, 5);          % ms (WM 10% range)
kap = 2 .^ linspace(-0.5, 0.5, 5);  % a.u. (factor-of-two variation)
wf = 0;                             % rad/ms

% Outer robustness criterion: check minima over loose bounds
T1_rob  = linspace(400, 1600, 13);  % ms 
T2_rob  = linspace(40, 160, 13);    % ms
kap_rob = 2 .^ linspace(-2, 2, 5);  % a.u.

% Other constant declarations
TE = 4;                             % ms
noise_var_1coil = 2.62e-07;         % a.u.
noise_var_ssos = noise_var_1coil/2; % High-SNR regime approximation

a = 0.1;                            % T1 relative importance parameter
b = 1;                              % T2 relative importance parameter
c = 0;                              % kap relative importance parameter

tol = 0.001;                        % Global minima tolerance
time_comp = 0;                      % Toggle 
pr = 0;                             % Toggle print on/off
savedat = 0;                        % Toggle saving .mat files on/off
disp_scale = [0 5];                 % Logarithmic display scaling

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Scan Sequence One: 1 SPGR, 1 DESS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Grid definitions 
% sig_T1_11   = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% sig_T2_11   = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% sig_kap_11  = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% Psi_11      = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% Psi_idx_11  = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% Psi_tild_11 = NaN(length(TR1_spgr), length(TR1_dess), length(spgr1), length(dess1));
% 
% % Covariance matrix: scaled identity with M = 3
% Sigma_inv = (1 ./ noise_var_ssos) * speye(3);
% 
% %% Step One: Min-max CRLB analysis over parameter grid
% for i_s1 = 1:length(spgr1)
%     spgr_11 = [spgr1(i_s1)]';
%     
%     for i_d1 = 1:length(dess1)
%         dess_11 = [dess1(i_d1)]';
%         
%         for i_TRs1 = 1:length(TR1_spgr)
%             TRs_11 = [TR1_spgr(i_TRs1)]';
%             
%             for i_TRd1 = 1:length(TR1_dess)
%                 TRd_11 = [TR1_dess(i_TRd1)]';
%                 scan_time = sum(TRs_11) + sum(TRd_11); 
%                 
%                 % Inner maximization: store max SDs over (T1,T2,kap) range
%                 worst_sigT1  = NaN(length(T1), length(T2), length(kap));
%                 worst_sigT2  = NaN(length(T1), length(T2), length(kap));
%                 worst_sigkap = NaN(length(T1), length(T2), length(kap));
% 
%                 for i_t1 = 1:length(T1)
%                 for i_t2 = 1:length(T2)
%                 for i_kap = 1:length(kap)
%                     [~, ~, worst_sigT1(i_t1, i_t2, i_kap),...
%                         worst_sigT2(i_t1, i_t2, i_kap),...
%                         worst_sigkap(i_t1, i_t2, i_kap)]...
%                         = norm_crlb_dess_4parm(T1(i_t1), T2(i_t2), kap(i_kap),...
%                         wf, spgr_11, dess_11, TRs_11, TRd_11, TE, Sigma_inv, time_comp);
%                 end
%                 end
%                 end
% 
%                 % Store the worst-case Psi and its corresponding (T1, T2, kap) index
%                 worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
%                 [Psi_11(i_TRs1, i_TRd1, i_s1, i_d1), Psi_idx_11(i_TRs1,...
%                     i_TRd1, i_s1, i_d1)] = max(worst_Psi(:));
%                 Psi_tild_11(i_TRs1, i_TRd1, i_s1, i_d1)...
%                     = Psi_11(i_TRs1, i_TRd1, i_s1, i_d1) * sqrt(scan_time);
% 
%                 % Extract this index and use it to save corresponding sigmas
%                 [T1_idx, T2_idx, kap_idx] = ind2sub([length(T1), length(T2), length(kap)],...
%                     Psi_idx_11(i_TRs1, i_TRd1, i_s1, i_d1));
%                 sig_T1_11(i_TRs1, i_TRd1, i_s1, i_d1)  = worst_sigT1(T1_idx, T2_idx, kap_idx);
%                 sig_T2_11(i_TRs1, i_TRd1, i_s1, i_d1)  = worst_sigT2(T1_idx, T2_idx, kap_idx);
%                 sig_kap_11(i_TRs1, i_TRd1, i_s1, i_d1) = worst_sigkap(T1_idx, T2_idx, kap_idx);
%             end
%         end
%     end
% end
% 
% % Find the indices of (multiple) Psi_tild_11 minima, to within a tolerance
% Psimin_idx_11 = find( ((Psi_tild_11 - min(Psi_tild_11(:))) ./ min(Psi_tild_11(:))) <= tol);
% num_min_11 = length(Psimin_idx_11);
% 
% %% Step Two: Select one of these minima based on robustness over wider ranges 
% % Indices of each of the "num_min_11" total minima
% TRs1_idx_11 = NaN(num_min_11, 1);
% TRd1_idx_11 = NaN(num_min_11, 1);
% s1_idx_11  = NaN(num_min_11, 1);
% d1_idx_11  = NaN(num_min_11, 1);
% 
% sigT1_worst_11  = NaN(num_min_11, length(T1_rob), length(T2_rob), length(kap_rob));
% sigT2_worst_11  = NaN(num_min_11, length(T1_rob), length(T2_rob), length(kap_rob));
% sigkap_worst_11 = NaN(num_min_11, length(T1_rob), length(T2_rob), length(kap_rob));
% 
% for i_min = 1:num_min_11
%     % Convert the 1D index to ND-grid indices
%     [TRs1_idx_11(i_min), TRd1_idx_11(i_min), s1_idx_11(i_min),...
%         d1_idx_11(i_min)] = ind2sub(size(Psi_11), Psimin_idx_11(i_min));
%     
%     % Store the minimizing parameters
%     TRs1min_11 = TR1_spgr(TRs1_idx_11(i_min));
%     TRd1min_11 = TR1_dess(TRd1_idx_11(i_min));
%     s1min_11   = spgr1(s1_idx_11(i_min));
%     d1min_11   = dess1(d1_idx_11(i_min));
%     
% %     % Extract the relevant worst-case (T1, T2, kap) indices 
% %     [T1min_idx_11, T2min_idx_11, kapmin_idx_11] ...
% %         = ind2sub([length(T1), length(T2), length(kap)],...
% %         Psi_idx_11(TRs1_idx_11(i_min), TRd1_idx_11(i_min),...
% %         s1_idx_11(i_min), d1_idx_11(i_min)));
%     
%     % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
%     % Note now that we have turned on time compensation!
%     for i_t1 = 1:length(T1_rob)
%     for i_t2 = 1:length(T2_rob)
%     for i_kap = 1:length(kap_rob)
%         [~, ~, sigT1_worst_11(i_min, i_t1, i_t2, i_kap),...
%             sigT2_worst_11(i_min, i_t1, i_t2, i_kap),...
%             sigkap_worst_11(i_min, i_t1, i_t2, i_kap)] = ...
%             norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
%             wf, [s1min_11]', [d1min_11]', [TRs1min_11]', [TRd1min_11]',...
%             TE, Sigma_inv, 1);
%     end
%     end
%     end
% end
% 
% % Compute the worst-case Psi_tild_11 value over the wider range
% Psi_tild_min_11 = a*sigT1_worst_11 + b*sigT2_worst_11 + c*sigkap_worst_11; 
% [Psi_tild_maxmin_11, wide_idx_11] = max(reshape(Psi_tild_min_11,...
%     [num_min_11 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);
% 
% % Extract the parameters that minimize Psi_tild_maxmin_11 over local minima
% % Psi_tild_star - min(Psi_tild(:)) is the degradation from narrow->wide range
% [Psi_tild_star_11, i_star_11] = min(Psi_tild_maxmin_11); 
% TRs1star_11 = TRs1_idx_11(i_star_11);
% TRd1star_11 = TRd1_idx_11(i_star_11);
% s1star_11   = s1_idx_11(i_star_11);
% d1star_11   = d1_idx_11(i_star_11); 
% 
% % Save useful outputs
% % Remember, Psi = a*sig_T1 + b*sig_T2 + c*sig_kap
% if (savedat)
%     save('1SPGR,1DESS_v3.mat', 'sig_T1_11', 'sig_T2_11', 'sig_kap_11',...
%         'Psi_11', 'Psi_tild_11', 'Psi_idx_11', 'a', 'b', 'c',...
%         'TR1_spgr', 'TR2_spgr', 'TR1_dess', 'TR2_dess', 'spgr1',...
%         'spgr2', 'dess1', 'dess2', 'T1', 'T2', 'kap');
% end
% 
% %% Graphical Output: 2D Heatmaps for 1 SPGR, 1 DESS
% % T1 sigma vs. (spgr1/dess1)
% figure(1); imagesc(180*spgr1/pi, 180*dess1/pi, ...
%     squeeze(log10(sig_T1_11(TRs1star_11, TRd1star_11, :, :)))', disp_scale);...
%     xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
%     title('log(T1 sigma): 1 SPGR, 1 DESS'); axis xy square; colormap('hot'); 
%     if (pr) print -depsc sig_t1_vs_s1d1_11.eps, end;
% 
% % T2 sigma vs. (spgr1/dess1)
% figure(2); imagesc(180*spgr1/pi, 180*dess1/pi, ...
%     squeeze(log10(sig_T2_11(TRs1star_11, TRd1star_11, :, :)))', disp_scale);...
%     xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
%     title('log(T2 sigma): 1 SPGR, 1 DESS'); axis xy square; colormap('hot');
%     if (pr) print -depsc sig_t2_vs_s1d1_11.eps, end;
% 
% % Psi vs. (spgr1/dess1)
% figure(3); imagesc(180*spgr1/pi, 180*dess1/pi, ...
%     squeeze(log10(Psi_11(TRs1star_11, TRd1star_11, :, :)))', disp_scale);...
%     xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
%     title('log(Psi): 1 SPGR, 1 DESS'); axis xy square; colormap('hot');
%     if (pr) print -depsc psi_vs_s1d1_11.eps, end;
%     
% % Psi_tild vs. (spgr1/dess1)
% figure(4); imagesc(180*spgr1/pi, 180*dess1/pi, ...
%     squeeze(log10(Psi_tild_11(TRs1star_11, TRd1star_11, :, :)))', disp_scale);...
%     xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
%     title('log(Psi-tilda): 1 SPGR, 1 DESS'); axis xy square; colormap('hot');
%     if (pr) print -depsc psitild_vs_s1d1_11.eps, end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Two: 0 SPGR, 2 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid definitions 
sig_T1_02   = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));
sig_T2_02   = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));
sig_kap_02  = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));
Psi_02      = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));
Psi_idx_02  = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));
Psi_tild_02 = NaN(length(TR1_dess), length(TR2_dess), length(dess1), length(dess2));

% Covariance matrix: scaled identity with M = 4
Sigma_inv = (1 ./ noise_var_ssos) * speye(4);

%% Step One: Min-max CRLB analysis over parameter grid
spgr_02 = []';
TRs_02 = []';
for i_d1 = 1:length(dess1)
for i_d2 = 1:length(dess2)
    dess_02 = [dess1(i_d1) dess2(i_d2)]';

    for i_TRd1 = 1:length(TR1_dess)
    for i_TRd2 = 1:length(TR2_dess)
        TRd_02 = [TR1_dess(i_TRd1) TR2_dess(i_TRd2)]';
        scan_time = sum(TRs_02) + sum(TRd_02);

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

        % Store the worst-case Psi and its corresponding (T1, T2, kap) index
        worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
        [Psi_02(i_TRd1, i_TRd2, i_d1, i_d2), Psi_idx_02(i_TRd1, i_TRd2,...
            i_d1, i_d2)] = max(worst_Psi(:));
        Psi_tild_02(i_TRd1, i_TRd2, i_d1, i_d2)...
            = Psi_02(i_TRd1, i_TRd2, i_d1, i_d2) * sqrt(scan_time);

        % Extract this index and use it to save corresponding sigmas
        [T1_idx, T2_idx, kap_idx] = ind2sub([length(T1), length(T2), length(kap)],...
            Psi_idx_02(i_TRd1, i_TRd2, i_d1, i_d2));
        sig_T1_02(i_TRd1, i_TRd2, i_d1, i_d2)  = worst_sigT1(T1_idx, T2_idx, kap_idx);
        sig_T2_02(i_TRd1, i_TRd2, i_d1, i_d2)  = worst_sigT2(T1_idx, T2_idx, kap_idx);
        sig_kap_02(i_TRd1, i_TRd2, i_d1, i_d2) = worst_sigkap(T1_idx, T2_idx, kap_idx);
    end
    end
end
end

% Find the indices of (multiple) Psi_tild_02 minima, to within a tolerance
Psimin_idx_02 = find( ((Psi_tild_02 - min(Psi_tild_02(:))) ./ min(Psi_tild_02(:))) <= tol);
num_min_02 = length(Psimin_idx_02);

%% Step Two: Select one of these minima based on robustness over kappa 
% Indices of each of the "num_min_02" total minima
TRd1_idx_02 = NaN(num_min_02, 1);
TRd2_idx_02 = NaN(num_min_02, 1);
d1_idx_02  = NaN(num_min_02, 1);
d2_idx_02  = NaN(num_min_02, 1);

sigT1_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_02  = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_02 = NaN(num_min_02, length(T1_rob), length(T2_rob), length(kap_rob));

for i_min = 1:num_min_02
    % Convert the 1D index to ND-grid indices
    [TRd1_idx_02(i_min), TRd2_idx_02(i_min), d1_idx_02(i_min),...
        d2_idx_02(i_min)] = ind2sub(size(Psi_02), Psimin_idx_02(i_min));
    
    % Store the minimizing parameters
    TRd1min_02 = TR1_dess(TRd1_idx_02(i_min));
    TRd2min_02 = TR2_dess(TRd2_idx_02(i_min));
    d1min_02   = dess1(d1_idx_02(i_min));
    d2min_02   = dess2(d2_idx_02(i_min));
    
%     % Extract the relevant worst-case (T1, T2, kap) indices 
%     [T1min_idx_02, T2min_idx_02, kapmin_idx_02]...
%         = ind2sub([length(T1), length(T2), length(kap)],...
%         Psi_idx_02(TRd1_idx_02(i_min), TRd2_idx_02(i_min),...
%         d1_idx_02(i_min), d2_idx_02(i_min)));
    
    % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
    % Note now that we have turned on time compensation!
    for i_t1 = 1:length(T1_rob)
    for i_t2 = 1:length(T2_rob)
    for i_kap = 1:length(kap_rob)
        [~, ~, sigT1_worst_02(i_min, i_t1, i_t2, i_kap),...
            sigT2_worst_02(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_02(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, []', [d1min_02 d2min_02]', []',...
            [TRd1min_02 TRd2min_02]', TE, Sigma_inv, 1);
    end
    end
    end
end

% Compute the worst-case Psi_tild_02 value over the wider range
Psi_tild_min_02 = a*sigT1_worst_02 + b*sigT2_worst_02 + c*sigkap_worst_02; 
[Psi_tild_maxmin_02, wide_idx_02] = max(reshape(Psi_tild_min_02,...
    [num_min_02 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameter indices that minimize Psi_tild_maxmin_02 over local minima
% Psi_tild_star - min(Psi_tild(:)) is the degradation from narrow->wide range
[Psi_tild_star_02, i_star_02] = min(Psi_tild_maxmin_02); 
TRd1star_02 = TRd1_idx_02(i_star_02);
TRd2star_02 = TRd2_idx_02(i_star_02);
d1star_02   = d1_idx_02(i_star_02); 
d2star_02   = d2_idx_02(i_star_02);

% Save useful outputs
if (savedat)
    save('0SPGR,2DESS_0,1,0.mat', 'sig_T1_02', 'sig_T2_02', 'sig_kap_02',...
        'Psi_02', 'Psi_tild_02', 'Psi_idx_02', 'a', 'b', 'c',...
        'TR1_spgr', 'TR2_spgr', 'TR1_dess', 'TR2_dess', 'spgr1',...
        'spgr2', 'dess1', 'dess2', 'T1', 'T2', 'kap');
end

%% Graphical Output: 2D Heatmaps for 0 SPGR, 2 DESS
% T1 sigma vs. (dess1/dess2)
figure(5); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sig_T1_02(TRd1star_02, TRd2star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t1_vs_d1d2_02_0,1,0.eps', end;

% T2 sigma vs. (dess1/dess2)
figure(6); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sig_T2_02(TRd1star_02, TRd2star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t2_vs_d1d2_02_0,1,0.eps', end;

% Psi vs. (dess1/dess2)
figure(7); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(Psi_02(TRd1star_02, TRd2star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psi_vs_d1d2_02_0,1,0.eps', end;
    
% Psi_tild vs. (dess1/dess2)
figure(8); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(Psi_tild_02(TRd1star_02, TRd2star_02, :, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi-tilda): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psitild_vs_d1d2_02_0,1,0.eps', end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Three: 2 SPGR, 1 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid definitions 
sig_T1_21   = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));
sig_T2_21   = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));
sig_kap_21  = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));
Psi_21      = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));
Psi_idx_21  = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));
Psi_tild_21 = NaN(length(TR1_spgr), length(TR2_spgr), length(TR1_dess),...
    length(spgr1), length(spgr2), length(dess1));

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
        
            for i_TRd1 = 1:length(TR1_dess)
                TRd_21 = [TR1_dess(i_TRd1)]';
                scan_time = sum(TRs_21) + sum(TRd_21);

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
                        wf, spgr_21, dess_21, TRs_21, TRd_21, TE, Sigma_inv, time_comp);
                end
                end
                end
                
                % Store the worst-case Psi and its corresponding (T1, T2, kap) index
                worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap;
                [Psi_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1),...
                    Psi_idx_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1)]...
                    = max(worst_Psi(:));
                Psi_tild_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1)...
                    = Psi_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1) * sqrt(scan_time);
                
                % Extract this index and use it to save corresponding sigmas
                [T1_idx, T2_idx, kap_idx] = ind2sub([length(T1), length(T2), length(kap)],...
                    Psi_idx_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1));
                sig_T1_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1)...  
                    = worst_sigT1(T1_idx, T2_idx, kap_idx);
                sig_T2_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1)...  
                    = worst_sigT2(T1_idx, T2_idx, kap_idx);
                sig_kap_21(i_TRs1, i_TRs2, i_TRd1, i_s1, i_s2, i_d1)...
                    = worst_sigkap(T1_idx, T2_idx, kap_idx);
            end
        end
        end
    end
end
end

% Find the indices of (multiple) Psi_tild_21 minima, to within a tolerance
Psimin_idx_21 = find( ((Psi_tild_21 - min(Psi_tild_21(:))) ./ min(Psi_tild_21(:))) <= tol);
num_min_21 = length(Psimin_idx_21);

%% Step Two: Select one of these minima based on robustness over kappa 
% Indices of each of the "num_min_21" total minima
TRs1_idx_21 = NaN(num_min_21, 1);
TRs2_idx_21 = NaN(num_min_21, 1);
TRd1_idx_21 = NaN(num_min_21, 1);
s1_idx_21   = NaN(num_min_21, 1);
s2_idx_21   = NaN(num_min_21, 1);
d1_idx_21   = NaN(num_min_21, 1);

sigT1_worst_21  = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_21  = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_21 = NaN(num_min_21, length(T1_rob), length(T2_rob), length(kap_rob));

for i_min = 1:num_min_21
    % Convert the 1D index to ND-grid indices
    [TRs1_idx_21(i_min), TRs2_idx_21(i_min), TRd1_idx_21(i_min),...
        s1_idx_21(i_min), s2_idx_21(i_min), d1_idx_21(i_min)]...
        = ind2sub(size(Psi_21), Psimin_idx_21(i_min));
    
    % Store the minimizing parameters
    TRs1min_21 = TR1_spgr(TRs1_idx_21(i_min));
    TRs2min_21 = TR2_spgr(TRs2_idx_21(i_min));
    TRd1min_21 = TR1_dess(TRd1_idx_21(i_min));
    s1min_21   = spgr1(s1_idx_21(i_min));
    s2min_21   = spgr2(s2_idx_21(i_min));
    d1min_21   = dess1(d1_idx_21(i_min));
    
%     % Extract the relevant worst-case (T1, T2, kap) indices 
%     [T1min_idx_21, T2min_idx_21, kapmin_idx_21]...
%         = ind2sub([length(T1), length(T2), length(kap)],...
%         Psi_idx_21(TRs1_idx_21(i_min), TRs2_idx_21(i_min), TRd1_idx_21(i_min),...
%         s1_idx_21(i_min), s2_idx_21(i_min), d1_idx_21(i_min)));
    
    % Evaluate CRLB at minimizers, over (T1_rob, T2_rob, kap_rob) ranges
    % Note now that we have turned on time compensation!
    for i_t1 = 1:length(T1_rob)
    for i_t2 = 1:length(T2_rob)
    for i_kap = 1:length(kap)
        [~, ~, sigT1_worst_21(i_min, i_t1, i_t2, i_kap),...
            sigT2_worst_21(i_min, i_t1, i_t2, i_kap),...
            sigkap_worst_21(i_min, i_t1, i_t2, i_kap)] = ...
            norm_crlb_dess_4parm(T1_rob(i_t1), T2_rob(i_t2), kap_rob(i_kap),...
            wf, [s1min_21 s2min_21]', [d1min_21]',...
            [TRs1min_21 TRs2min_21]', [TRd1min_21]', TE, Sigma_inv, 1);
    end
    end
    end
end

% Compute the worst-case Psi_21 value over the wider range
Psi_tild_min_21 = a*sigT1_worst_21 + b*sigT2_worst_21 + c*sigkap_worst_21; 
[Psi_tild_maxmin_21, wide_idx_21] = max(reshape(Psi_tild_min_21,...
    [num_min_21 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameters that minimize Psi_maxmin_21 over local minima
% Psi_star - min(Psi(:)) is the degradation from kappa variation
[Psi_star_21, i_star_21] = min(Psi_tild_maxmin_21); 
TRs1star_21 = TRs1_idx_21(i_star_21);
TRs2star_21 = TRs2_idx_21(i_star_21);
TRd1star_21 = TRd1_idx_21(i_star_21);
s1star_21   = s1_idx_21(i_star_21);
s2star_21   = s2_idx_21(i_star_21);
d1star_21   = d1_idx_21(i_star_21);

% Save useful outputs
if (savedat)
    save('2SPGR,1DESS_0,1,0.mat', 'sig_T1_21', 'sig_T2_21', 'sig_kap_21',...
        'Psi_21', 'Psi_tild_21', 'Psi_idx_21', 'a', 'b', 'c',...
        'TR1_spgr', 'TR2_spgr', 'TR1_dess', 'TR2_dess',...
        'spgr1', 'spgr2', 'dess1', 'dess2', 'T1', 'T2', 'kap');
end

%% Graphical Output: 2D Heatmaps for 2 SPGR, 1 DESS
% T1 sigma vs. (spgr1/spgr2)
figure(9); imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    sig_T1_21(TRs1star_21, TRs2star_21, TRd1star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t1_vs_s1s2_21_0,1,0.eps', end;

% T2 sigma vs. (spgr1/spgr2)
figure(10); imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    sig_T2_21(TRs1star_21, TRs2star_21, TRd1star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t2_vs_s1s2_21_0,1,0.eps', end;

% Psi vs. (spgr1/spgr2)
figure(11); imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, TRd1star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psi_vs_s1s2_21_0,1,0.eps', end;
    
% Psi_tild vs. (spgr1/spgr2)
figure(12); imagesc(180*spgr1/pi, 180*spgr2/pi, squeeze(log10(...
    Psi_tild_21(TRs1star_21, TRs2star_21, TRd1star_21, :, :, d1star_21)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(Psi-tild): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psitild_vs_s1s2_21_0,1,0.eps', end;
    

% T1 sigma vs. (spgr1/dess1)
figure(13); imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    sig_T1_21(TRs1star_21, TRs2star_21, TRd1star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t1_vs_s1d1_21_0,1,0.eps', end;

% T2 sigma vs. (spgr1/dess1)
figure(14); imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    sig_T2_21(TRs1star_21, TRs2star_21, TRd1star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t2_vs_s1d1_21_0,1,0.eps', end;

% Psi vs. (spgr1/dess1)
figure(15); imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, TRd1star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psi_vs_s1d1_21_0,1,0.eps', end;

% Psi_tild vs. (spgr1/dess1)
figure(16); imagesc(180*spgr1/pi, 180*dess1/pi, squeeze(log10(...
    Psi_tild_21(TRs1star_21, TRs2star_21, TRd1star_21, :, s2star_21, :)))', disp_scale);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi-tild): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psitild_vs_s1d1_21_0,1,0.eps', end;
    
    
% T1 sigma vs. (spgr2/dess1)
figure(17); imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    sig_T1_21(TRs1star_21, TRs2star_21, TRd1star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t1_vs_s2d1_21_0,1,0.eps', end;

% T2 sigma vs. (spgr2/dess1)
figure(18); imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    sig_T2_21(TRs1star_21, TRs2star_21, TRd1star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'sig_t2_vs_s2d1_21_0,1,0.eps', end;

% Psi vs. (spgr2/dess1)
figure(19); imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    Psi_21(TRs1star_21, TRs2star_21, TRd1star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psi_vs_s2d1_21_0,1,0.eps', end;
    
% Psi-tild vs. (spgr2/dess1)
figure(20); imagesc(180*spgr2/pi, 180*dess1/pi, squeeze(log10(...
    Psi_tild_21(TRs1star_21, TRs2star_21, TRd1star_21, s1star_21, :, :)))', disp_scale);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi-tild): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc 'psitild_vs_s2d1_21_0,1,0.eps', end;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan Sequence Four: 0 SPGR, 3 DESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid definitions 
sig_T1_03   = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));
sig_T2_03   = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));
sig_kap_03  = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));
Psi_03      = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));
Psi_idx_03  = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));
Psi_tild_03 = NaN(length(TR1_dess), length(TR2_dess), length(TR3_dess),...
    length(dess1), length(dess2), length(dess3));

% Covariance matrix: scaled identity with M = 6
Sigma_inv = (1 ./ noise_var_ssos) * speye(6);

%% Step One: Min-max CRLB analysis over parameter grid
spgr_03 = []';
TRs_03 = []';
for i_d1 = 1:length(dess1)
for i_d2 = 1:length(dess2)
for i_d3 = 1:length(dess3)
    dess_03 = [dess1(i_d1) dess2(i_d2) dess3(i_d3)]';
    
    for i_TRd1 = 1:length(TR1_dess)
    for i_TRd2 = 1:length(TR2_dess)
    for i_TRd3 = 1:length(TR3_dess)
        TRd_03 = [TR1_dess(i_TRd1) TR2_dess(i_TRd2) TR3_dess(i_TRd3)]';
        scan_time = sum(TRs_03) + sum(TRd_03);
        
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
                wf, spgr_03, dess_03, TRs_03, TRd_03, TE, Sigma_inv, time_comp);
        end
        end
        end
        
        % Store the worst-case Psi and its corresponding (T1, T2, kap) index
        worst_Psi = a*worst_sigT1 + b*worst_sigT2 + c*worst_sigkap; 
        [Psi_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3),...
            Psi_idx_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3)] = max(worst_Psi(:));
        Psi_tild_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3)...
            = Psi_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3) * sqrt(scan_time);
        
        % Extract this index and use it to save corresponding sigmas
        [T1_idx, T2_idx, kap_idx] = ind2sub([length(T1), length(T2), length(kap)],...
            Psi_idx_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3));
        sig_T1_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3)...
            = worst_sigT1(T1_idx, T2_idx, kap_idx);
        sig_T2_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3)...
            = worst_sigT2(T1_idx, T2_idx, kap_idx);
        sig_kap_03(i_TRd1, i_TRd2, i_TRd3, i_d1, i_d2, i_d3)...
            = worst_sigkap(T1_idx, T2_idx, kap_idx);
    end
    end
    end
end
end
end

% Find the indices of (multiple) Psi_tild_03 minima, to within a tolerance
Psimin_idx_03 = find( ((Psi_tild_03 - min(Psi_tild_03(:))) ./ min(Psi_tild_03(:))) <= tol);
num_min_03 = length(Psimin_idx_03);
        
%% Step Two: Select one of these minima based on robustness over kappa 
% Indices of each of the "num_min_03" total minima
TRd1_idx_03 = NaN(num_min_03, 1);
TRd2_idx_03 = NaN(num_min_03, 1);
TRd3_idx_03 = NaN(num_min_03, 1);
d1_idx_03   = NaN(num_min_03, 1);
d2_idx_03   = NaN(num_min_03, 1);
d3_idx_03   = NaN(num_min_03, 1);

sigT1_worst_03  = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));
sigT2_worst_03  = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));
sigkap_worst_03 = NaN(num_min_03, length(T1_rob), length(T2_rob), length(kap_rob));
        
for i_min = 1:num_min_03
    % Convert the 1D index to ND-grid indices
    [TRd1_idx_03(i_min), TRd2_idx_03(i_min), TRd3_idx_03(i_min),...
        d1_idx_03(i_min), d2_idx_03(i_min), d3_idx_03(i_min)]...
        = ind2sub(size(Psi_03), Psimin_idx_03(i_min));
    
    % Store the minimizing parameters
    TRd1min_03 = TR1_dess(TRd1_idx_03(i_min));
    TRd2min_03 = TR2_dess(TRd2_idx_03(i_min));
    TRd3min_03 = TR3_dess(TRd3_idx_03(i_min));
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
Psi_tild_min_03 = a*sigT1_worst_03 + b*sigT2_worst_03 + c*sigkap_worst_03; 
[Psi_tild_maxmin_03, wide_idx_03] = max(reshape(Psi_tild_min_03,...
    [num_min_03 length(T1_rob)*length(T2_rob)*length(kap_rob)]), [], 2);

% Extract the parameter indices that minimize Psi_tild_maxmin_03 over local minima
% Psi_tild_star - min(Psi_tild(:)) is the degradation from narrow->wide range
[Psi_tild_star_03, i_star_03] = min(Psi_tild_maxmin_03); 
TRd1star_03 = TRd1_idx_03(i_star_03);
TRd2star_03 = TRd2_idx_03(i_star_03);
TRd3star_03 = TRd3_idx_03(i_star_03);
d1star_03   = d1_idx_03(i_star_03); 
d2star_03   = d2_idx_03(i_star_03);
d3star_03   = d3_idx_03(i_star_03);

% Save useful outputs
if (savedat)
    save('0SPGR,3DESS_0.1,1,0.mat', 'sig_T1_03', 'sig_T2_03', 'sig_kap_03',...
        'Psi_03', 'Psi_tild_03', 'Psi_idx_03', 'a', 'b', 'c',...
        'TR1_spgr', 'TR2_spgr', 'TR1_dess', 'TR2_dess', 'TR3_dess', 'spgr1',...
        'spgr2', 'dess1', 'dess2', 'dess3', 'T1', 'T2', 'kap');
end

%% Graphical Output: 2D Heatmaps for 0 SPGR, 3 DESS
% T1 sigma vs. (dess1/dess2)
figure(21); imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(sig_T1_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t1_vs_d1d2_03_0.1,1,0.eps', end;
    
% T2 sigma vs. (dess1/dess2)
figure(22); imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(sig_T2_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t2_vs_d1d2_03_0.1,1,0.eps', end;

% Psi vs. (dess1/dess2)
figure(23); imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psi_vs_d1d2_03_0.1,1,0.eps', end;

% Psi_tild vs. (dess1/dess2)
figure(24); imagesc(180*dess1/pi, 180*dess2/pi, squeeze(log10(Psi_tild_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, :, d3star_03)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi-tild): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psitild_vs_d1d2_03_0.1,1,0.eps', end;


% T1 sigma vs. (dess1/dess3)
figure(25); imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(sig_T1_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t1_vs_d1d3_03_0.1,1,0.eps', end;
    
% T2 sigma vs. (dess1/dess3)
figure(26); imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(sig_T2_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t2_vs_d1d3_03_0.1,1,0.eps', end;

% Psi vs. (dess1/dess3)
figure(27); imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psi_vs_d1d3_03_0.1,1,0.eps', end;

% Psi_tild vs. (dess1/dess3)
figure(28); imagesc(180*dess1/pi, 180*dess3/pi, squeeze(log10(Psi_tild_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, :, d2star_03, :)))', disp_scale);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi-tild): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psitild_vs_d1d3_03_0.1,1,0.eps', end;


% T1 sigma vs. (dess2/dess3)
figure(29); imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(sig_T1_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t1_vs_d2d3_03_0.1,1,0.eps', end;
    
% T2 sigma vs. (dess2/dess3)
figure(30); imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(sig_T2_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'sig_t2_vs_d2d3_03_0.1,1,0.eps', end;

% Psi vs. (dess2/dess3)
figure(31); imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(Psi_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psi_vs_d2d3_03_0.1,1,0.eps', end;

% Psi_tild vs. (dess2/dess3)
figure(32); imagesc(180*dess2/pi, 180*dess3/pi, squeeze(log10(Psi_tild_03(...
    TRd1star_03, TRd2star_03, TRd3star_03, d1star_03, :, :)))', disp_scale);...
    xlabel('DESS flip2 (deg)'); ylabel('DESS flip3 (deg)'); colorbar;...
    title('log(Psi-tild): 0 SPGR, 3 DESS'); axis xy square; colormap('hot');...
    if (pr) print -depsc 'psitild_vs_d2d3_03_0.1,1,0.eps', end;


