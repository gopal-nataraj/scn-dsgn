%% Header Script for Regularized M0s/T1/T2/kappa Recon 
% 11/07/2014, in vivo, brain, patient: J.F. Nielsen
% Additional code added for ISMRM 2014: empirical std. dev. study
%
% Written by: Gopal Nataraj
% Modified 2015-05-19: tighter kappa ranges

%% Raw data extraction
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
addpath('Calculations/DESS/');
addpath('Helper/');
addpath('JFNielsen/common');
addpath('JFNielsen/img');

% Global TE value must be set for all scans
TE_global = 4.67;                   % ms

% Imaging Parameters
nx = 240; ny = 240; nz = 6; nc = 32;    
opflip_all = [5:5:90]';             % degrees
flip_d_all = opflip_all * pi/180;   % radians
nfd_all = length(flip_d_all);
TRd_all = 17.28 * ones(nfd_all, 1);     % ms
TEd_all = TE_global * ones(nfd_all, 1); % ms
wf = 0;                             % Zero off-resonance
sl = 6;                             % Latter slice chosen due to wrap-around
pr = 1;                             % Toggle printing on/off

% Load DESS dat from .mat files 
% Choose just one medioaxial slice
% Use ir_mri_coil_combine.m to merge nc-channel coil data 
if (~exist('yp_im_all', 'var'))
    try
        addpath('Brain_11,07,14/');
        load(sprintf('ims_coil_comb_slice%u.mat', sl));
        crop_lr = [46:205]; nx = length(crop_lr);
        crop_ap = [15:224]; ny = length(crop_ap);
        yp_im_all = yp_im;
        ym_im_all = ym_im; 
        clear yp_im ym_im;
    catch
        dir = pwd; 
        cd('/Volumes/General Storage/Documents/Research/DESS/2014,11,07_brain_32ch');

        yp_im_all = NaN(nx, ny, nfd_all);
        ym_im_all = NaN(nx, ny, nfd_all);
        for l = 1:nfd_all
            load(sprintf('dess_jfn_07Nov2014_frame%u', l));
            yp_coil = squeeze(ims_echo1(:,:,sl,:));             % [nx ny nc]
            ym_coil = squeeze(ims_echo2(:,:,sl,:));             % [nx ny nc]
            yp_im_all(:,:,l) = ir_mri_coil_combine(yp_coil);    % [nx ny nfd]
            ym_im_all(:,:,l) = ir_mri_coil_combine(ym_coil);    % [nx ny nfd]
            clear ims_echo1 ims_echo2 yp_coil ym_coil;
        end
        cd(dir);

        % Reorder the flip angles
        % Before, they were [5, 90, 15, 80, 25, 70,...]
        yp_im_all(:,:,2:2:end) = flipdim(yp_im_all(:,:,2:2:end), 3);
        ym_im_all(:,:,2:2:end) = flipdim(ym_im_all(:,:,2:2:end), 3);

        % Take transpose to orient anterior --> posterior
        % Flip image to compensate second echo inverted in readout direction
        % This is now coil-combined complex data
        yp_im_all = permute(yp_im_all, [2 1 3]);
        ym_im_all = permute(flipdim(ym_im_all, 1), [2 1 3]);
        
        % Crop image to a tighter FOV for display 
        crop_lr = [46:205]; nx = length(crop_lr);
        crop_ap = [15:224]; ny = length(crop_ap);
        yp_im_all = yp_im_all(crop_lr, crop_ap, :);
        ym_im_all = ym_im_all(crop_lr, crop_ap, :);
        save(sprintf('ims_coil_comb_slice%u.mat', sl), 'yp_im_all', 'ym_im_all');
    end
end

% Even though data complex, turn is_mag on because most of the signal
% energy is placed on the real axis by ir_mri_coi_combine.m
is_mag = 1;

% OPTIONAL: Select only two flip angles 
indices = [2 7]';                           % 10,35 degrees
% indices = [3 8]';                           % 15,40 degrees
opflip_all = opflip_all(indices);           % degrees
flip_d_all = opflip_all  * pi/180;          % radians
nfd_all = length(flip_d_all);
yp_im_all = yp_im_all(:,:,indices);
ym_im_all = ym_im_all(:,:,indices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART ONE: "GROUND TRUTH" FROM ALL 18 FLIP ANGLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total number of datasets
M_all = 2*nfd_all;

%% Dictionary creation (could be precomputed)
T1 = logspace(1, 3.5, 100);
T2 = logspace(0, 2.5, 100);
kap = linspace(0.9, 1.1, 11);
% kap = 2 .^ linspace(-0.5, 0.5, 21);
% kap = 2 .^ linspace(-0.25, 0.25, 11);
D_all = NaN(M_all, length(T1), length(T2), length(kap));
K = length(T1) * length(T2) * length(kap);

for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        for s = 1:length(kap)
            % DESS dictionary component
            for a = 1:nfd_all
                [D_all(a, t1, t2, s), D_all(nfd_all+a, t1, t2, s)] = ...
                    dess_fun(1, T1(t1), T2(t2), kap(s), flip_d_all(a),...
                    TRd_all(a), TEd_all(a), wf, is_mag);
            end
        end
    end
end
D_all = reshape(D_all, [M_all K]);

%% Maximum-likelihood Estimation
% Dictionary-based estimation via variable-projection method
weights = ones(M_all, 1);
W = spdiags(weights, 0, M_all, M_all);      % Weighting Matrix
y = reshape(permute(cat(3, yp_im_all, ym_im_all), [3 1 2]), [M_all nx*ny]);

tic; 
maxProd = zeros(1, nx*ny);
idx = zeros(1, nx*ny);
for k = 1:K
    % Compute kth inner product
    hess = abs(D_all(:,k)' * W * D_all(:,k));
    ytild = D_all(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;

    % If the kth inner product is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    idx(update) = k;
end
time_ml = toc;

% Extract indices for maximum-likelihood maps
[t1_idx, t2_idx, s_idx] = ind2sub([length(T1) length(T2) length(kap)], idx);
T1_ml = reshape(T1(t1_idx), [nx ny]); 
T2_ml = reshape(T2(t2_idx), [nx ny]);
kap_ml = reshape(kap(s_idx), [nx ny]);

% For now, take true value wf=0
wf_ml = zeros(nx, ny);

%% Preprocessing and Masking
% Project Images to within range
T1max = 5000;       T1_ml = min(T1_ml, T1max);
T1min = 5;          T1_ml = max(T1_ml, T1min);
T2max = 2000;       T2_ml = min(T2_ml, T2max);
T2min = 5;          T2_ml = max(T2_ml, T2min);

% Create tight and loose masks
tight_mask = imfill(squeeze(abs(yp_im_all(:,:,ceil(nfd_all/2)+1)) >...
    0.05*max(col(abs(yp_im_all(:,:,ceil(nfd_all/2)+1))))), 'holes');
loose_mask = imdilate(tight_mask, strel('disk', 10)); 

% Set voxels inside loose mask but outside tight mask to mean
T1_ml(~tight_mask & loose_mask) = mean(col(T1_ml));   
T2_ml(~tight_mask & loose_mask) = mean(col(T2_ml));   
wf_ml(~tight_mask & loose_mask) = mean(col(wf_ml)); 
kap_ml(~tight_mask & loose_mask) = mean(col(kap_ml));

% Set voxels outside both tight and loose masks to zero
T1_ml(~loose_mask) = 0;  
T2_ml(~loose_mask) = 0;  
kap_ml(~loose_mask) = 0;
wf_ml(~loose_mask) = 0;  

% Median filtering
T1_med = medfilt2(T1_ml);  
T2_med = medfilt2(T2_ml); 
kap_med = medfilt2(kap_ml);

%% Regularized, Joint M0s,T1,T2,kap Estimation
% Define iteration parameters
n_outer = 20;    
niterM = 50;   
niter1 = 100;   
niter2 = 100;   
niterK = 50;

tolM = 10^-6;   
tol1 = 10^-7;  
tol2 = 10^-7;  
tolk = 10^-4; 
disp = 0; 

% Define regularizers, Rm, R1, and R2
betaM = nfd_all * 2^-8; 
beta1 = nfd_all * 2^-15;%2^-14; 
beta2 = nfd_all * 2^-15;%2^-15;
betaK = nfd_all * 2^-4;

deltaM = 10^-2; 
delta1 = 10^-1;
delta2 = 10^-2;

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);
Rk = Reg1(loose_mask, 'pot_arg', {'quad'}, 'beta', betaK);

% Fictitious SPGR placeholders
nfs = 0;
ys_im = zeros(nx, ny, nfs);
flip_s = zeros(nfs, 1);
TRs = zeros(nfs, 1);
TEs = zeros(nfs, 1);

% Reshape data inputs to work with 3D implementation
ys_3D = permute(ys_im, [1 2 4 3]);
yp_3D = permute(yp_im_all, [1 2 4 3]);
ym_3D = permute(ym_im_all, [1 2 4 3]);

% Regularized Reconstruction
tic;
[M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
    mri_spgr_dess_m0st1t2kap_map(T1_med, T2_med, kap_med, wf_ml, flip_s,...
    flip_d_all, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
    TRs, TRd_all, TEs, TEd_all, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
    niter2, niterK, tolM, tol1, tol2, tolk, is_mag, disp);
time_reg = toc; 

% % Sanity check: note that a continuous kap_init needs to be used...
% [M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
%     mri_spgr_dess_m0st1t2kap_map(T1_true, T2_true, kap_ml, wf_ml, flip_s,...
%     flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
%     TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
%     niter2, niterK, tolM, tol1, tol2, tolk, is_mag, disp);

% Postprocessing for display
% For real data, can use this mask
reg_mask = imfill(abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:))), 'holes');

% Remove pixels outside reg_mask
M0s_ml(~reg_mask) = 0; M0s_med(~reg_mask) = 0; M0s_reg(~reg_mask) = 0;
T1_ml(~reg_mask) = 0;  T1_med(~reg_mask) = 0;  T1_reg(~reg_mask) = 0;
T2_ml(~reg_mask) = 0;  T2_med(~reg_mask) = 0;  T2_reg(~reg_mask) = 0; 
kap_ml(~reg_mask) = 0; kap_med(~reg_mask) = 0; kap_reg(~reg_mask) = 0;

% % OPTIONAL: Remove CSF pixels 
% mask_roi = (T2_reg > 20) & (T2_reg < 80);     % Range of interest
% M0s_ml(~mask_roi) = 0; M0s_med(~mask_roi) = 0; M0s_reg(~mask_roi) = 0;
% T1_ml(~mask_roi) = 0;  T1_med(~mask_roi) = 0;  T1_reg(~mask_roi) = 0;
% T2_ml(~mask_roi) = 0;  T2_med(~mask_roi) = 0;  T2_reg(~mask_roi) = 0; 
% kap_ml(~mask_roi) = 0; kap_med(~mask_roi) = 0; kap_reg(~mask_roi) = 0;

%% For comparison, Method-of-Moments T2 Estimation (biased)
E2sqrd = NaN(nx, ny);
for xx = 1:nx
    for yy = 1:ny
        Y = squeeze(ym_im_all(xx, yy, :));
        A = squeeze(yp_im_all(xx, yy, :));
        C = [A Y];
        sig_small = min(svds(C'*C));
        E2sqrd(xx, yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
    end
end
T2_mom = -2*(mean(TRd_all)-TE_global) ./ log(E2sqrd);
T2_mom(~reg_mask) = 0;
% T2_mom(~mask_roi) = 0;

%% Summary Statistics 
% Define ROIs (WM, GM, CSF)
ctrX = [100 85  72];
ctrY = [51  37  86];
rad  = [5   2   4 ];
nROI = length(rad);

% Extract data means and standard deviations from ROI
yp_mean = NaN(nROI, nfd_all); yp_std = NaN(nROI, nfd_all);
ym_mean = NaN(nROI, nfd_all); ym_std = NaN(nROI, nfd_all);

M0s_ml_mean = NaN(nROI, 1);     M0s_ml_std = NaN(nROI, 1);
M0s_reg_mean = NaN(nROI, 1);    M0s_reg_std = NaN(nROI, 1);

T1_ml_mean = NaN(nROI, 1);      T1_ml_std = NaN(nROI, 1);
T1_reg_mean = NaN(nROI, 1);     T1_reg_std = NaN(nROI, 1);

T2_ml_mean = NaN(nROI, 1);      T2_ml_std = NaN(nROI, 1);
T2_reg_mean = NaN(nROI, 1);     T2_reg_std = NaN(nROI, 1);
T2_mom_mean = NaN(nROI, 1);     T2_mom_std = NaN(nROI, 1);

kap_ml_mean = NaN(nROI, 1);     kap_ml_std = NaN(nROI, 1);
kap_reg_mean = NaN(nROI, 1);    kap_reg_std = NaN(nROI, 1);

for r = 1:nROI
    for l = 1:nfd_all
        [yp_mean(r,l), yp_std(r,l)] ...
            = multiMeans(yp_im_all(:,:,l), [ctrX(r) ctrY(r)], rad(r));
        [ym_mean(r,l), ym_std(r,l)] ...
            = multiMeans(ym_im_all(:,:,l), [ctrX(r) ctrY(r)], rad(r));
    end
    
    [M0s_ml_mean(r), M0s_ml_std(r)] ...
        = multiMeans(abs(M0s_ml),   [ctrX(r) ctrY(r)], rad(r));
    [M0s_reg_mean(r), M0s_reg_std(r)] ...
        = multiMeans(abs(M0s_reg),  [ctrX(r) ctrY(r)], rad(r));
    
    [T1_ml_mean(r), T1_ml_std(r)] ...
        = multiMeans(T1_ml,     [ctrX(r) ctrY(r)], rad(r));
    [T1_reg_mean(r), T1_reg_std(r)] ...
        = multiMeans(T1_reg,    [ctrX(r) ctrY(r)], rad(r));
    
    [T2_ml_mean(r), T2_ml_std(r)] ...
        = multiMeans(T2_ml,     [ctrX(r) ctrY(r)], rad(r));
    [T2_reg_mean(r), T2_reg_std(r)] ...
        = multiMeans(T2_reg,    [ctrX(r) ctrY(r)], rad(r));
    [T2_mom_mean(r), T2_mom_std(r)] ...
        = multiMeans(T2_mom,    [ctrX(r) ctrY(r)], rad(r));
    
    [kap_ml_mean(r), kap_ml_std(r)] ...
        = multiMeans(kap_ml,    [ctrX(r) ctrY(r)], rad(r));
    [kap_reg_mean(r), kap_reg_std(r)] ...
        = multiMeans(kap_reg,   [ctrX(r) ctrY(r)], rad(r));
end

% Export summary statistics to file
if (pr)
    fid = fopen('Brain_11,07,14_summary_10,35', 'a');
    fprintf(fid, 'M0s Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, 'Vial%u:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', r,...
            M0s_ml_mean(r), char(177), M0s_ml_std(r),...
            M0s_reg_mean(r), char(177), M0s_reg_std(r));
    end

    fprintf(fid, '\n\nT1 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, 'Vial%u:\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\n', r,...
            T1_ml_mean(r), char(177), T1_ml_std(r),...
            T1_reg_mean(r), char(177), T1_reg_std(r));
    end

    fprintf(fid, '\n\nT2 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\t\tMoM\n');
    for r = 1:nROI
        fprintf(fid, 'Vial%u:\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\n', r,...
            T2_ml_mean(r), char(177), T2_ml_std(r),...
            T2_reg_mean(r), char(177), T2_reg_std(r),...
            T2_mom_mean(r), char(177), T2_mom_std(r));
    end

    fprintf(fid, '\n\nKappa Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, 'Vial%u:\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\n', r,...
            kap_ml_mean(r), char(177), kap_ml_std(r),...
            kap_reg_mean(r), char(177), kap_reg_std(r));
    end
    fclose(fid);
end

%% Images, RMSE, and Comparisons
% ROI rectangles (x,y,w,h)
posWM  = [ctrX(1)-rad(1) ctrY(1)-rad(1) 2*rad(1)+1 2*rad(1)+1];
posGM  = [ctrX(2)-rad(2) ctrY(2)-rad(2) 2*rad(2)+1 2*rad(2)+1];
posCSF = [ctrX(3)-rad(3) ctrY(3)-rad(3) 2*rad(1)+1 2*rad(3)+1];

% M0s Images 
figure; im('notick', abs(M0s_ml),  [0 4], 'cbar', ' ');
if (pr), print('-deps','M0s_ml_10,35.eps'); end
figure; im('notick', abs(M0s_med), [0 4], 'cbar', ' ');
if (pr), print('-deps','M0s_med_10,35.eps'); end
figure; im('notick', abs(M0s_reg), [0 4], 'cbar', ' ');
if (pr), print('-deps','M0s_reg_10,35.eps'); end

% T1 Images
t1min = 400; t1max = 1200;
figure; im('notick', T1_ml,  [t1min t1max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T1_ml_10,35.eps'); end
figure; im('notick', T1_med, [t1min t1max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T1_med_10,35.eps'); end
figure; im('notick', T1_reg, [t1min t1max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T1_reg_10,35.eps'); end

% T2 Images w/ ROIs
t2min = 20; t2max = 100;
figure; im('notick', T2_mom, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_mom_10,35.eps'); end
figure; im('notick', T2_ml,  [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_ml_10,35.eps'); end 
figure; im('notick', T2_med, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_med_10,35.eps'); end
figure; im('notick', T2_reg, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 1], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 0.9 0], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_reg_10,35.eps'); end

% kap Images
figure; im('notick', kap_ml,  [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_ml_10,35.eps'); end
figure; im('notick', kap_med, [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_med_10,35.eps'); end
figure; im('notick', kap_reg, [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_reg_10,35.eps'); end

% Cost vs. Iteration
figure; hold on;
scatter(1:4:4*n_outer, cost(2:4:end), 'bo');
scatter(2:4:4*n_outer, cost(3:4:end), 'ro');
scatter(3:4:4*n_outer, cost(4:4:end), 'go');
scatter(4:4:4*n_outer, cost(5:4:end), 'mo');
plot(0:4*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M0s update', 'T1 update', 'T2 update', 'kap update');
% print('-depsc','cost_vs_iter.eps');

%% Sanity check: does the data fit the model?
% Note that we can only do this because TRs are fixed across scans
% Use average values in ROIs to get corresponding ML parameters
% This prevents data noise from being amplified in estimating curves
y_mean = [yp_mean ym_mean]'; 
maxProd_mean = zeros(1, nROI);
idx_mean = zeros(1, nROI);
for k = 1:K
    % Compute kth inner product
    hess_mean = abs(D_all(:,k)' * W * D_all(:,k));
    ytild_mean = D_all(:,k)' * W * y_mean / sqrt(hess_mean);
    newProd_mean = abs(ytild_mean).^2;
    
    % If the kth inner product is largest, save k
    update_mean = newProd_mean > maxProd_mean;
    maxProd_mean(update_mean) = newProd_mean(update_mean);
    idx_mean(update_mean) = k;
end

% Extract indices for ML maps
[t1_idx_mean, t2_idx_mean, s_idx_mean] ...
    = ind2sub([length(T1) length(T2) length(kap)], idx_mean);
T1_sig_mean = T1(t1_idx_mean);
T2_sig_mean = T2(t2_idx_mean);
kap_sig_mean = kap(s_idx_mean);

% M0s initial guess (from dictionary)
M0s_sig_mean = NaN(nROI, 1);
for q = 1:length(idx_mean)
    M0s_sig_mean(q) = (D_all(:,idx_mean(q))' * y_mean(:,q)) ./...
        (D_all(:,idx_mean(q))' * D_all(:,idx_mean(q)));
end

% Signal model vs. Data
for r = 1:nROI
    flip_model = [0:90] * (pi/180) * kap_sig_mean(r);
    flip_data = flip_d_all * kap_sig_mean(r);
    [S1mean, S2mean] = dess_fun(M0s_sig_mean(r), T1_sig_mean(r),...
        T2_sig_mean(r), 1, flip_model, mean(TRd_all), TE_global, wf, 1);
    
    figure; hold on;
    errorbar(flip_data * (180/pi), abs(yp_mean(r,:)), yp_std(r,:), 'b*');
    plot(flip_model * (180/pi), abs(S1mean), 'c');

    errorbar(flip_data * (180/pi), abs(ym_mean(r,:)), ym_std(r,:), 'r*');
    plot(flip_model * (180/pi), abs(S2mean), 'm');
    hold off;

    legend('SSFP-FID Data', 'SSFP-FID ML Model',...
        'SSFP-Echo Data', 'SSFP-Echo ML Model');
    tit = sprintf('NIST ROI %d: DESS Model vs. Data, (T1,T2,kap) = (%0.1fms,%0.2fms,%0.2f)',...
        r, T1_reg_mean(r), T2_reg_mean(r), kap_reg_mean(r)); title(tit); 
    xlabel('(Compensated) flip angle (deg)'); 
    ylabel('Magnitude Signal (a.u.)');
    if (pr), print('-depsc', sprintf('sig_vs_10,35_ROI%u.eps', r)); end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART TWO: TESTS WITH DIFFERENT SCAN COMBINATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save "ground truth" variables 
M0s_ml_18flip = M0s_ml;
T1_ml_18flip  = T1_ml;
T2_ml_18flip  = T2_ml;
kap_ml_18flip = kap_ml;

M0s_reg_18flip = M0s_reg;
T1_reg_18flip  = T1_reg;
T2_reg_18flip  = T2_reg;
kap_reg_18flip = kap_reg;

% Create masks for comparisons
mask_roi = (T2_reg_18flip > 20) & (T2_reg_18flip < 100); % Range of interest
mask_wm = false(nx, ny);  mask_wm(ctrX(1)-rad(1):ctrX(1)+rad(1),  ctrY(1)-rad(1):ctrY(1)+rad(1)) = true;
mask_gm = false(nx, ny);  mask_gm(ctrX(2)-rad(2):ctrX(2)+rad(2),  ctrY(2)-rad(2):ctrY(2)+rad(2)) = true;
mask_csf = false(nx, ny); mask_csf(ctrX(3)-rad(3):ctrX(3)+rad(3), ctrY(3)-rad(3):ctrY(3)+rad(3)) = true;
mask_wmgm = mask_wm | mask_gm;

% Create various performance-criterion placeholders
% Save std dev in WM, GM, WM/GM, CSF for ML/REG
% Save RMSE in ROI, WM, GM, WM/GM, CSF
std_t1_wm_ml    = Inf(nfd_all, nfd_all);
std_t1_wm_reg   = Inf(nfd_all, nfd_all);
std_t1_gm_ml    = Inf(nfd_all, nfd_all);
std_t1_gm_reg   = Inf(nfd_all, nfd_all);
std_t1_csf_ml   = Inf(nfd_all, nfd_all);
std_t1_csf_reg  = Inf(nfd_all, nfd_all);

std_t2_wm_ml    = Inf(nfd_all, nfd_all);
std_t2_wm_reg   = Inf(nfd_all, nfd_all);
std_t2_gm_ml    = Inf(nfd_all, nfd_all);
std_t2_gm_reg   = Inf(nfd_all, nfd_all);
std_t2_csf_ml   = Inf(nfd_all, nfd_all);
std_t2_csf_reg  = Inf(nfd_all, nfd_all);

rmse_t1_wm_ml    = Inf(nfd_all, nfd_all);
rmse_t1_wm_reg   = Inf(nfd_all, nfd_all);
rmse_t1_gm_ml    = Inf(nfd_all, nfd_all);
rmse_t1_gm_reg   = Inf(nfd_all, nfd_all);
rmse_t1_wmgm_ml  = Inf(nfd_all, nfd_all);
rmse_t1_wmgm_reg = Inf(nfd_all, nfd_all);
rmse_t1_csf_ml   = Inf(nfd_all, nfd_all);
rmse_t1_csf_reg  = Inf(nfd_all, nfd_all);
rmse_t1_roi_ml   = Inf(nfd_all, nfd_all);
rmse_t1_roi_reg  = Inf(nfd_all, nfd_all);

rmse_t2_wm_ml    = Inf(nfd_all, nfd_all);
rmse_t2_wm_reg   = Inf(nfd_all, nfd_all);
rmse_t2_gm_ml    = Inf(nfd_all, nfd_all);
rmse_t2_gm_reg   = Inf(nfd_all, nfd_all);
rmse_t2_wmgm_ml  = Inf(nfd_all, nfd_all);
rmse_t2_wmgm_reg = Inf(nfd_all, nfd_all);
rmse_t2_csf_ml   = Inf(nfd_all, nfd_all);
rmse_t2_csf_reg  = Inf(nfd_all, nfd_all);
rmse_t2_roi_ml   = Inf(nfd_all, nfd_all);
rmse_t2_roi_reg  = Inf(nfd_all, nfd_all);

% RMSE helper function
rmse  = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ numel(tru));

% Run the helper script for all combinations of two flip angles
for i = 1:nfd_all
    for j = 1:nfd_all
        if i > j
            printf('\n\nCopied Flips (%u,%u)', opflip_all(i), opflip_all(j));
            
            std_t1_wm_ml(i,j)    = std_t1_wm_ml(j,i);
            std_t1_wm_reg(i,j)   = std_t1_wm_reg(j,i);
            std_t1_gm_ml(i,j)    = std_t1_gm_ml(j,i);
            std_t1_gm_reg(i,j)   = std_t1_gm_reg(j,i);
            std_t1_csf_ml(i,j)   = std_t1_csf_ml(j,i);
            std_t1_csf_reg(i,j)  = std_t1_csf_reg(j,i);

            std_t2_wm_ml(i,j)    = std_t2_wm_ml(j,i);
            std_t2_wm_reg(i,j)   = std_t2_wm_reg(j,i);
            std_t2_gm_ml(i,j)    = std_t2_gm_ml(j,i);
            std_t2_gm_reg(i,j)   = std_t2_gm_reg(j,i);
            std_t2_csf_ml(i,j)   = std_t2_csf_ml(j,i);
            std_t2_csf_reg(i,j)  = std_t2_csf_reg(j,i);

            rmse_t1_wm_ml(i,j)    = rmse_t1_wm_ml(j,i);
            rmse_t1_wm_reg(i,j)   = rmse_t1_wm_reg(j,i);
            rmse_t1_gm_ml(i,j)    = rmse_t1_gm_ml(j,i);
            rmse_t1_gm_reg(i,j)   = rmse_t1_gm_reg(j,i);
            rmse_t1_wmgm_ml(i,j)  = rmse_t1_wmgm_ml(j,i);
            rmse_t1_wmgm_reg(i,j) = rmse_t1_wmgm_reg(j,i);
            rmse_t1_csf_ml(i,j)   = rmse_t1_csf_ml(j,i);
            rmse_t1_csf_reg(i,j)  = rmse_t1_csf_reg(j,i);
            rmse_t1_roi_ml(i,j)   = rmse_t1_roi_ml(j,i);
            rmse_t1_roi_reg(i,j)  = rmse_t1_roi_reg(j,i);
            
            rmse_t2_wm_ml(i,j)    = rmse_t2_wm_ml(j,i);
            rmse_t2_wm_reg(i,j)   = rmse_t2_wm_reg(j,i);
            rmse_t2_gm_ml(i,j)    = rmse_t2_gm_ml(j,i);
            rmse_t2_gm_reg(i,j)   = rmse_t2_gm_reg(j,i);
            rmse_t2_wmgm_ml(i,j)  = rmse_t2_wmgm_ml(j,i);
            rmse_t2_wmgm_reg(i,j) = rmse_t2_wmgm_reg(j,i);
            rmse_t2_csf_ml(i,j)   = rmse_t2_csf_ml(j,i);
            rmse_t2_csf_reg(i,j)  = rmse_t2_csf_reg(j,i);
            rmse_t2_roi_ml(i,j)   = rmse_t2_roi_ml(j,i);
            rmse_t2_roi_reg(i,j)  = rmse_t2_roi_reg(j,i);
        elseif i < j
            printf('\n\nComputing Flips (%u,%u)', opflip_all(i), opflip_all(j));
            
            % Select the right flip angles and data
            flip1 = flip_d_all(i);
            flip2 = flip_d_all(j);
            flip_d = [flip1 flip2]';
            opflip = flip_d * 180/pi;
            nfd = length(flip_d);
            TRd = 17.28 * ones(nfd, 1);     % ms
            TEd = TE_global * ones(nfd, 1); % ms
            yp_im = cat(3, yp_im_all(:,:,i), yp_im_all(:,:,j));
            ym_im = cat(3, ym_im_all(:,:,i), ym_im_all(:,:,j));
            
            % Dictionary extraction
            rows = [i j i+nfd_all j+nfd_all]';
            D = D_all(rows, :);
            M = 2*nfd;
            
            % Run the helper script for ML/reg recon
            brain_2014_11_07_dess_helper();
            
            % Compare T1 and T2 maps
            std_t1_wm_ml(i,j)    = std(T1_ml(mask_wm));
            std_t1_wm_reg(i,j)   = std(T1_reg(mask_wm));
            std_t1_gm_ml(i,j)    = std(T1_ml(mask_gm));
            std_t1_gm_reg(i,j)   = std(T1_reg(mask_gm));
            std_t1_csf_ml(i,j)   = std(T1_ml(mask_csf));
            std_t1_csf_reg(i,j)  = std(T1_reg(mask_csf));

            std_t2_wm_ml(i,j)    = std(T2_ml(mask_wm));
            std_t2_wm_reg(i,j)   = std(T2_reg(mask_wm));
            std_t2_gm_ml(i,j)    = std(T2_ml(mask_gm));
            std_t2_gm_reg(i,j)   = std(T2_reg(mask_gm));
            std_t2_csf_ml(i,j)   = std(T2_ml(mask_csf));
            std_t2_csf_reg(i,j)  = std(T2_reg(mask_csf));

            rmse_t1_wm_ml(i,j)    = rmse(T1_ml(mask_wm), T1_ml_18flip(mask_wm));
            rmse_t1_wm_reg(i,j)   = rmse(T1_reg(mask_wm), T1_reg_18flip(mask_wm));
            rmse_t1_gm_ml(i,j)    = rmse(T1_ml(mask_gm), T1_ml_18flip(mask_gm));
            rmse_t1_gm_reg(i,j)   = rmse(T1_reg(mask_gm), T1_reg_18flip(mask_gm));
            rmse_t1_wmgm_ml(i,j)  = rmse(T1_ml(mask_wmgm), T1_ml_18flip(mask_wmgm));
            rmse_t1_wmgm_reg(i,j) = rmse(T1_reg(mask_wmgm), T1_reg_18flip(mask_wmgm));
            rmse_t1_csf_ml(i,j)   = rmse(T1_ml(mask_csf), T1_ml_18flip(mask_csf));
            rmse_t1_csf_reg(i,j)  = rmse(T1_reg(mask_csf), T1_reg_18flip(mask_csf));
            rmse_t1_roi_ml(i,j)   = rmse(T1_ml(mask_roi), T1_ml_18flip(mask_roi));
            rmse_t1_roi_reg(i,j)  = rmse(T1_reg(mask_roi), T1_reg_18flip(mask_roi));
            
            rmse_t2_wm_ml(i,j)    = rmse(T2_ml(mask_wm), T2_ml_18flip(mask_wm));
            rmse_t2_wm_reg(i,j)   = rmse(T2_reg(mask_wm), T2_reg_18flip(mask_wm));
            rmse_t2_gm_ml(i,j)    = rmse(T2_ml(mask_gm), T2_ml_18flip(mask_gm));
            rmse_t2_gm_reg(i,j)   = rmse(T2_reg(mask_gm), T2_reg_18flip(mask_gm));
            rmse_t2_wmgm_ml(i,j)  = rmse(T2_ml(mask_wmgm), T2_ml_18flip(mask_wmgm));
            rmse_t2_wmgm_reg(i,j) = rmse(T2_reg(mask_wmgm), T2_reg_18flip(mask_wmgm));
            rmse_t2_csf_ml(i,j)   = rmse(T2_ml(mask_csf), T2_ml_18flip(mask_csf));
            rmse_t2_csf_reg(i,j)  = rmse(T2_reg(mask_csf), T2_reg_18flip(mask_csf));
            rmse_t2_roi_ml(i,j)   = rmse(T2_ml(mask_roi), T2_ml_18flip(mask_roi));
            rmse_t2_roi_reg(i,j)  = rmse(T2_reg(mask_roi), T2_reg_18flip(mask_roi));
        else 
            printf('\n\nSkipping Flips (%u,%u)', opflip_all(i), opflip_all(j));
        end
    end
end

% Take max standard deviations over WM/GM
std_t1_wmgm_ml  = max(std_t1_wm_ml, std_t1_gm_ml);
std_t2_wmgm_ml  = max(std_t2_wm_ml, std_t2_gm_ml);
std_t1_wmgm_reg = max(std_t1_wm_reg, std_t1_gm_reg);
std_t2_wmgm_reg = max(std_t2_wm_reg, std_t2_gm_reg);

% Save the reconstruction results
save('empirical_stats', 'std_t*', 'rmse_t*');

% Display figures for selection
figure; imagesc(opflip_all, opflip_all, std_t2_wm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_wm_ml');
figure; imagesc(opflip_all, opflip_all, std_t2_gm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_gm_ml');
figure; imagesc(opflip_all, opflip_all, std_t2_wmgm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_wmgm_ml');
figure; imagesc(opflip_all, opflip_all, std_t2_csf_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_csf_ml');

figure; imagesc(opflip_all, opflip_all, std_t2_wm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_wm_reg');
figure; imagesc(opflip_all, opflip_all, std_t2_gm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_gm_reg');
figure; imagesc(opflip_all, opflip_all, std_t2_wmgm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_wmgm_reg');
figure; imagesc(opflip_all, opflip_all, std_t2_csf_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('std_t2_csf_reg');

figure; imagesc(opflip_all, opflip_all, rmse_t2_wm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_wm_ml');
figure; imagesc(opflip_all, opflip_all, rmse_t2_gm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_gm_ml');
figure; imagesc(opflip_all, opflip_all, rmse_t2_wmgm_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_wmgm_ml');
figure; imagesc(opflip_all, opflip_all, rmse_t2_csf_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_csf_ml');
figure; imagesc(opflip_all, opflip_all, rmse_t2_roi_ml, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_roi_ml');

figure; imagesc(opflip_all, opflip_all, rmse_t2_wm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_wm_reg');
figure; imagesc(opflip_all, opflip_all, rmse_t2_gm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_gm_reg');
figure; imagesc(opflip_all, opflip_all, rmse_t2_wmgm_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_wmgm_reg');
figure; imagesc(opflip_all, opflip_all, rmse_t2_csf_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_csf_reg');
figure; imagesc(opflip_all, opflip_all, rmse_t2_roi_reg, [0 20]);...
    colormap hot; colorbar; axis xy; title('rmse_t2_roi_reg'); 

% Fancy outputs for abstract
xmin = [15 45]; 
ymin = [45 15];
figure; hold on; imagesc(opflip_all, opflip_all, std_t2_wmgm_ml, [0 20]);...
    scatter(xmin, ymin, 150, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
    xlabel('DESS flip1 (deg)', 'FontSize', 14);...
    ylabel('DESS flip2 (deg)', 'FontSize', 14);...
    title('Empirical ML T2 Max Standard Deviation (ms)', 'FontSize', 14);...
    axis xy square; axis([5 90 5 90]); colormap('hot'); colorbar;...
    if (pr) print('-depsc', 'std_t2_wmgm_ml.eps'), end;

xmin = [20 50]; 
ymin = [50 20];
figure; hold on; imagesc(opflip_all, opflip_all, rmse_t2_wmgm_ml, [0 20]);...
    scatter(xmin, ymin, 150, 'p', 'MarkerEdgeColor', 'c', 'MarkerFaceColor', 'c');...
    xlabel('DESS flip1 (deg)', 'FontSize', 14);...
    ylabel('DESS flip2 (deg)', 'FontSize', 14);...
    title('Empirical ML T2 RMSD (ms)', 'FontSize', 14);...
    axis xy square; axis([5 90 5 90]); colormap('hot'); colorbar;...
    if (pr) print('-depsc', 'rmsd_t2_wmgm_ml.eps'), end;