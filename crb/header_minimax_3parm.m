% Header File -- Min-max CRLB analysis
% Worst-case CRLB over T1/T2 ROI, varying TRs and flip angles w/ M0*
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

%% Define parameter grid: varying flip angles and TR values
M0s = 1;
T1 = linspace(400,1000,4);          % ms (WM/GM nominal range)
T2 = linspace(60,100,4);            % ms (WM/GM nominal range)
wf = 0;                             % rad/ms
TE = 5;                             % ms

TR_spgr = linspace(8,16,5);         % ms
TR_dess = linspace(10,20,6);        % ms

spgr1 = linspace(5,45,41) * pi/180; % rad
spgr2 = linspace(5,45,41) * pi/180; % rad
dess1 = linspace(5,85,41) * pi/180; % rad
dess2 = linspace(5,85,41) * pi/180; % rad

c = 0.1;                            % T2/T1 relative importance parameter
pr = 1;                             % Print on/off

%% Scan Sequence One: 1 SPGR, 1 DESS
% Grid definitions
sig_T1_11 = NaN(length(TR_spgr), length(TR_dess), length(spgr1), length(dess1));
sig_T2_11 = NaN(length(TR_spgr), length(TR_dess), length(spgr1), length(dess1));

% CRLB analysis over parameter grid
for i_s1 = 1:length(spgr1)
    flip_spgr = [spgr1(i_s1)]';
    
    for i_d1 = 1:length(dess1)
        flip_dess = [dess1(i_d1)]';
        
        for i_TRs = 1:length(TR_spgr)
        for i_TRd = 1:length(TR_dess)

            worst_sigT1 = NaN(length(T1), length(T2));
            worst_sigT2 = NaN(length(T1), length(T2));
            for i_t1 = 1:length(T1)
            for i_t2 = 1:length(T2)
                [~, ~, worst_sigT1(i_t1, i_t2), worst_sigT2(i_t1, i_t2)]...
                    = crlb_dess_m0star(M0s, T1(i_t1), T2(i_t2), wf,...
                    flip_spgr, flip_dess, TR_spgr(i_TRs), TR_dess(i_TRd), TE);
            end
            end

            sig_T1_11(i_TRs, i_TRd, i_s1, i_d1) = max(worst_sigT1(:));
            sig_T2_11(i_TRs, i_TRd, i_s1, i_d1) = max(worst_sigT2(:));
        end
        end
    end
end

% Compute objective function, Psi(sig_T1, sig_T2)
Psi_11 = c * sig_T1_11 + sig_T2_11;

% Find Psi_11 minimum-variance indices 
[TRs_11, TRd_11, s1_11, d1_11] = ind2sub(size(Psi_11), ...
    find(Psi_11 == min(Psi_11(:)), 1));

%% Scan Sequence Two: 0 SPGR, 2 DESS
% Grid definitions
sig_T1_02 = NaN(length(TR_dess), length(dess1), length(dess2));
sig_T2_02 = NaN(length(TR_dess), length(dess1), length(dess2));

% CRLB analysis over parameter grid
flip_spgr = []'; 
for i_d1 = 1:length(dess1)
    for i_d2 = 1:length(dess2)
        flip_dess = [dess1(i_d1) dess2(i_d2)]';
        
        for i_TRd = 1:length(TR_dess)
            
            worst_sigT1 = NaN(length(T1), length(T2));
            worst_sigT2 = NaN(length(T1), length(T2));
            for i_t1 = 1:length(T1)
            for i_t2 = 1:length(T2)
                [~, ~, worst_sigT1(i_t1, i_t2), worst_sigT2(i_t1, i_t2)]...
                    = crlb_dess_m0star(M0s, T1(i_t1), T2(i_t2), wf,...
                    flip_spgr, flip_dess, 0, TR_dess(i_TRd), TE);
            end
            end
            
            sig_T1_02(i_TRd, i_d1, i_d2) = max(worst_sigT1(:));
            sig_T2_02(i_TRd, i_d1, i_d2) = max(worst_sigT2(:));
        end
    end
end

% Compute objective function, Psi(sig_T1, sig_T2)
Psi_02 = c * sig_T1_02 + sig_T2_02;

% Find Psi_02 minimum-variance indices (omitting diagonal)
Psi_02_tmp = Psi_02;
for j = 1:size(Psi_02_tmp, 2)
    Psi_02_tmp(:,j,j) = Inf;    % Effectively don't consider diagonal
end
[TRd_02, d1_02, d2_02] = ind2sub(size(Psi_02_tmp), ...
    find(Psi_02_tmp == min(Psi_02_tmp(:)), 1));

%% Scan Sequence Three: 2 SPGR, 1 DESS
% Grid definitions
sig_T1_21 = NaN(length(TR_spgr), length(TR_dess), length(spgr1), ...
    length(spgr2), length(dess1));
sig_T2_21 = NaN(length(TR_spgr), length(TR_dess), length(spgr1), ...
    length(spgr2), length(dess1));

% CRLB analysis over parameter grid
for i_s1 = 1:length(spgr1)
for i_s2 = 1:length(spgr2)
    flip_spgr = [spgr1(i_s1) spgr2(i_s2)]';

    for i_d1 = 1:length(dess1)
        flip_dess = [dess1(i_d1)]';

        for i_TRs = 1:length(TR_spgr)
        for i_TRd = 1:length(TR_dess)

            worst_sigT1 = NaN(length(T1), length(T2));
            worst_sigT2 = NaN(length(T1), length(T2));
            for i_t1 = 1:length(T1)
            for i_t2 = 1:length(T2)
                [~, ~, worst_sigT1(i_t1, i_t2), worst_sigT2(i_t1, i_t2)]...
                    = crlb_dess_m0star(M0s, T1(i_t1), T2(i_t2), wf,...
                    flip_spgr, flip_dess, TR_spgr(i_TRs), TR_dess(i_TRd), TE);
            end
            end
            
            sig_T1_21(i_TRs, i_TRd, i_s1, i_s2, i_d1) = max(worst_sigT1(:));
            sig_T2_21(i_TRs, i_TRd, i_s1, i_s2, i_d1) = max(worst_sigT2(:));
        end
        end
    end
end
end

% Compute objective function, Psi(sig_T1, sig_T2)
Psi_21 = c * sig_T1_21 + sig_T2_21;

% Find Psi_21 minimum-variance indices
[TRs_21, TRd_21, s1_21, s2_21, d1_21] = ind2sub(size(Psi_21), ...
    find(Psi_21 == min(Psi_21(:)), 1));

%% Robustness Comparison

% Set tolerances
del_flip = 0.2;
del_TR = 0.2;

% Get corresponding indices
TRs_idx_11 = (TR_spgr >= (1-del_TR)*TR_spgr(TRs_11)) & ...
    (TR_spgr <= (1+del_TR)*TR_spgr(TRs_11));
TRs_idx_21 = (TR_spgr >= (1-del_TR)*TR_spgr(TRs_21)) & ...
    (TR_spgr <= (1+del_TR)*TR_spgr(TRs_21));

TRd_idx_11 = (TR_dess >= (1-del_TR)*TR_dess(TRd_11)) & ...
    (TR_dess <= (1+del_TR)*TR_dess(TRd_11));
TRd_idx_21 = (TR_dess >= (1-del_TR)*TR_dess(TRd_21)) & ...
    (TR_dess <= (1+del_TR)*TR_dess(TRd_21));
TRd_idx_02 = (TR_dess >= (1-del_TR)*TR_dess(TRd_02)) & ...
    (TR_dess <= (1+del_TR)*TR_dess(TRd_02));

s1_idx_11 = (spgr1 >= (1-del_flip)*spgr1(s1_11)) & ...
    (spgr1 <= (1+del_flip)*spgr1(s1_11));
s1_idx_21 = (spgr1 >= (1-del_flip)*spgr1(s1_21)) & ...
    (spgr1 <= (1+del_flip)*spgr1(s1_21));
s2_idx_21 = (spgr2 >= (1-del_flip)*spgr2(s2_21)) & ...
    (spgr2 <= (1+del_flip)*spgr2(s2_21));

d1_idx_11 = (dess1 >= (1-del_flip)*dess1(d1_11)) & ...
    (dess1 <= (1+del_flip)*dess1(d1_11));
d1_idx_21 = (dess1 >= (1-del_flip)*dess1(d1_21)) & ...
    (dess1 <= (1+del_flip)*dess1(d1_21));
d1_idx_02 = (dess1 >= (1-del_flip)*dess1(d1_02)) & ...
    (dess1 <= (1+del_flip)*dess1(d1_02));
d2_idx_02 = (dess2 >= (1-del_flip)*dess2(d2_02)) & ...
    (dess2 <= (1+del_flip)*dess2(d2_02));

% Find worse case values
sigT1_11_wst = max(col(sig_T1_11(TRs_idx_11, TRd_idx_11, s1_idx_11, d1_idx_11)));
sigT2_11_wst = max(col(sig_T2_11(TRs_idx_11, TRd_idx_11, s1_idx_11, d1_idx_11)));
Psi_11_wst = max(col(Psi_11(TRs_idx_11, TRd_idx_11, s1_idx_11, d1_idx_11)));

sigT1_02_wst = max(col(sig_T1_02(TRd_idx_02, d1_idx_02, d2_idx_02)));
sigT2_02_wst = max(col(sig_T2_02(TRd_idx_02, d1_idx_02, d2_idx_02)));
Psi_02_wst = max(col(Psi_02(TRd_idx_02, d1_idx_02, d2_idx_02)));

sigT1_21_wst = max(col(sig_T1_21(TRs_idx_21, TRd_idx_21, s1_idx_21, s2_idx_21, d1_idx_21)));
sigT2_21_wst = max(col(sig_T2_21(TRs_idx_21, TRd_idx_21, s1_idx_21, s2_idx_21, d1_idx_21)));
Psi_21_wst = max(col(Psi_21(TRs_idx_21, TRd_idx_21, s1_idx_21, s2_idx_21, d1_idx_21)));

% Prepare table for output
% TRspgr; TRdess; spgr1; spgr2; dess1; dess2; sigT1; sigT2; Psi; w_sigT1; w_sigT2; w_Psi
table = [...
    TR_spgr(TRs_11), TR_dess(TRd_11), spgr1(s1_11), NaN, dess1(d1_11), NaN,...
    sig_T1_11(TRs_11, TRd_11, s1_11, d1_11)...
    sig_T2_11(TRs_11, TRd_11, s1_11, d1_11)...
    Psi_11(TRs_11, TRd_11, s1_11, d1_11)...
    sigT1_11_wst, sigT2_11_wst, Psi_11_wst;...
    ...
    TR_spgr(TRs_21), TR_dess(TRd_21), spgr1(s1_21), spgr2(s2_21), dess1(d1_21), NaN,...
    sig_T1_21(TRs_21, TRd_21, s1_21, s2_21, d1_21)...
    sig_T2_21(TRs_21, TRd_21, s1_21, s2_21, d1_21)...
    Psi_21(TRs_21, TRd_21, s1_21, s2_21, d1_21)...
    sigT1_21_wst, sigT2_21_wst, Psi_21_wst;...
    ...
    NaN, TR_dess(TRd_02), NaN, NaN, dess1(d1_02), dess2(d2_02),...
    sig_T1_02(TRd_02, d1_02, d2_02)...
    sig_T2_02(TRd_02, d1_02, d2_02)...
    Psi_02(TRd_02, d1_02, d2_02)...
    sigT1_02_wst, sigT2_02_wst, Psi_02_wst];
table(:, 3:6) = (180/pi) * table(:, 3:6);
table(:, 7:end) = log10(table(:, 7:end));

%% Graphical Output: 2D Heatmaps for 1 SPGR, 1 DESS
% T1 sigma vs. (spgr1/dess1)
figure(1); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T1_11(TRs_11, TRd_11, :, :)))', [4 7]);...
    xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 1 SPGR, 1 DESS'); axis xy square; colormap('hot'); 
    if (pr) print -depsc sig_t1_vs_s1d1_11.eps, end;

% T2 sigma vs. (spgr1/dess1)
figure(2); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T2_11(TRs_11, TRd_11, :, :)))', [4 7]);...
    xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 1 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t2_vs_s1d1_11.eps, end;

% Psi vs. (spgr1/dess1)
figure(3); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(Psi_11(TRs_11, TRd_11, :, :)))', [4 7]);...
    xlabel('SPGR flip (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 1 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc psi_vs_s1d1_11.eps, end;

%% Graphical Output: 2D Heatmaps for 0 SPGR, 2 DESS
% T1 sigma vs. (dess1/dess2)
figure(4); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sig_T1_02(TRd_02, :, :)))', [4 7]);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t1_vs_d1d2_02.eps, end;

% T2 sigma vs. (dess1/dess2)
figure(5); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(sig_T2_02(TRd_02, :, :)))', [4 7]);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t2_vs_d1d2_02.eps, end;

% Psi vs. (dess1/dess2)
figure(6); imagesc(180*dess1/pi, 180*dess2/pi, ...
    squeeze(log10(Psi_02(TRd_02, :, :)))', [4 7]);...
    xlabel('DESS flip1 (deg)'); ylabel('DESS flip2 (deg)'); colorbar;...
    title('log(Psi): 0 SPGR, 2 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc psi_vs_d1d2_02.eps, end;
    
%% Graphical Output: 2D Heatmaps for 2 SPGR, 1 DESS
% T1 sigma vs. (spgr1/spgr2)
figure(7); imagesc(180*spgr1/pi, 180*spgr2/pi, ...
    squeeze(log10(sig_T1_21(TRs_21, TRd_21, :, :, d1_21)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t1_vs_s1s2_21.eps, end;

% T2 sigma vs. (spgr1/spgr2)
figure(8); imagesc(180*spgr1/pi, 180*spgr2/pi, ...
    squeeze(log10(sig_T2_21(TRs_21, TRd_21, :, :, d1_21)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t2_vs_s1s2_21.eps, end;

% Psi vs. (spgr1/spgr2)
figure(9); imagesc(180*spgr1/pi, 180*spgr2/pi, ...
    squeeze(log10(Psi_21(TRs_21, TRd_21, :, :, d1_21)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('SPGR flip2 (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc psi_vs_s1s2_21.eps, end;

% T1 sigma vs. (spgr1/dess1)
figure(10); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T1_21(TRs_21, TRd_21, :, s2_21, :)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t1_vs_s1d1_21.eps, end;

% T2 sigma vs. (spgr1/dess1)
figure(11); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T2_21(TRs_21, TRd_21, :, s2_21, :)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t2_vs_s1d1_21.eps, end;

% Psi vs. (spgr1/dess1)
figure(12); imagesc(180*spgr1/pi, 180*dess1/pi, ...
    squeeze(log10(Psi_21(TRs_21, TRd_21, :, s2_21, :)))', [4 7]);...
    xlabel('SPGR flip1 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc psi_vs_s1d1_21.eps, end;


% T1 sigma vs. (spgr2/dess1)
figure(13); imagesc(180*spgr2/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T1_21(TRs_21, TRd_21, s1_21, :, :)))', [4 7]);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T1 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t1_vs_s2d1_21.eps, end;

% T2 sigma vs. (spgr2/dess1)
figure(14); imagesc(180*spgr2/pi, 180*dess1/pi, ...
    squeeze(log10(sig_T2_21(TRs_21, TRd_21, s1_21, :, :)))', [4 7]);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(T2 sigma): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc sig_t2_vs_s2d1_21.eps, end;

% Psi vs. (spgr2/dess1)
figure(15); imagesc(180*spgr2/pi, 180*dess1/pi, ...
    squeeze(log10(Psi_21(TRs_21, TRd_21, s1_21, :, :)))', [4 7]);...
    xlabel('SPGR flip2 (deg)'); ylabel('DESS flip (deg)'); colorbar;...
    title('log(Psi): 2 SPGR, 1 DESS'); axis xy square; colormap('hot');
    if (pr) print -depsc psi_vs_s2d1_21.eps, end;
            