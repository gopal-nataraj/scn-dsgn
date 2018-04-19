%% Header Script for Regularized M0s/T1/T2/kap Reconstruction from Synthetic Data
% For comparison of (2,1) optimized SPGR/DESS with (0,2) profile
% Dictionary initialization method and joint estimation method constant
%
% Written by: Gopal Nataraj
% Last modified 2015-05-09
% todo: update using 3parm version

%% Load true parameter maps and set scan protocol
% Load digital phantom
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;
addpath('Calculations/SPGR/');
addpath('Calculations/DESS/');
addpath('Helper/');

% Get true values
slice = 95; [nx, ny, nz] = size(labels);
M0_true = double(squeeze(labels(:,:,slice)));
T1_true = double(squeeze(labels(:,:,slice)));
T2_true = double(squeeze(labels(:,:,slice)));
T2s_true = double(squeeze(labels(:,:,slice)));
for idx = 0:10
    f = mri_brainweb_params(idx);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1; 
    T2_true(T2_true == idx) = f.t2;
    T2s_true(T2s_true == idx) = f.t2s;
end
wf_true = 0;

% Create a true excitation profile map
[xx, yy] = ndgrid(linspace(-0.5, 0.5, nx), linspace(-0.5, 0.5, ny));
kap_true = 1.5 - 2*(xx.^2 + yy.^2); 
kap_true(M0_true == 0) = 0; 

% Global TE value must be set for all scans
TE_global = 4;                          % ms

% SPGR imaging parameters
flip_s = []' * pi/180;                  % rad
nfs = length(flip_s);  
TRs = []';                              % ms
TEs = TE_global * ones(nfs, 1);         % ms
SNR_s = 60;                             % dB

% DESS imaging parameters
% flip_d = [45 10 90]' * pi/180;          % rad
flip_d = [90 45 10]' * pi/180;
nfd = length(flip_d);                   
% TRd = [20 8 22]';                       % ms
TRd = [8 10 12]';
TEd = TE_global * ones(nfd, 1);         % ms (symmetric echoes)
SNR_d = 60;                             % dB

% Total number of datasets 
M = length(flip_s) + 2*length(flip_d);

%% Dictionary creation (could be precomputed)
T1 = logspace(1, 3.5, 50);
T2 = logspace(0, 2.5, 50);
kap = 2 .^ linspace(-0.5, 1, 50);
D = NaN(M, length(T1), length(T2), length(kap));
K = length(T1) * length(T2) * length(kap);

for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        for s = 1:length(kap)
            % SPGR dictionary component
            for a = 1:nfs
                D(a, t1, t2, s) = spgr_fun(1, T1(t1), kap(s),...
                    flip_s(a), TRs(a), TEs(a), wf_true, 1);
            end
            
            % DESS dictionary component
            for a = 1:nfd
                [D(nfs+a, t1, t2, s), D(nfs+nfd+a, t1, t2, s)] = ...
                    dess_fun(1, T1(t1), T2(t2), kap(s), flip_d(a),...
                    TRd(a), TEd(a), wf_true, 1);
            end
        end
    end
end
D = reshape(D, [M, K]);

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = T2_true ~= T2s_true; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TE_global ./ T2s_true(T2s_msk));

%% Synthesize Phantom SPGR Data
% Forward model: make true data
ys_true = NaN(nx, ny, nfs);
for a = 1:nfs
    ys_true(:,:,a) = fft2(spgr_fun(M0s_true, T1_true, kap_true,...
        flip_s(a), TRs(a), TEs(a), wf_true, 1));
end

% Add complex white gaussian noise
var_s_im = 1.31e-7; sigma_s_k = sqrt(var_s_im * (nx*ny));
% sigma_s = exp(-SNR_s/20) * norm(ys_true(:)) / sqrt(2*numel(ys_true));
ys = ys_true + sigma_s_k * (randn(size(ys_true)) + 1i * randn(size(ys_true)));
printm('snr = %g dB', 20*log(norm(ys_true(:)) / norm(col(ys_true-ys))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
ys_im = abs(ifft2(ys));

%% Synthesize Phantom DESS Data
yp_true = NaN(nx, ny, nfd);
ym_true = NaN(nx, ny, nfd);
for a = 1:nfd
    [S1true, S2true] = dess_fun(M0s_true, T1_true, T2_true, kap_true,...
        flip_d(a), TRd(a), TEd(a), wf_true, 1);
    yp_true(:,:,a) = fft2(S1true);
    ym_true(:,:,a) = fft2(S2true);
end

% Add complex white gaussian noise
var_p_im = 1.31e-7; sigma_p_k = sqrt(var_p_im * (nx*ny));
var_m_im = 1.31e-7; sigma_m_k = sqrt(var_m_im * (nx*ny));
% sigma_p = exp(-SNR_d/20) * norm(yp_true(:)) / sqrt(2*numel(yp_true));
% sigma_m = exp(-SNR_d/20) * norm(ym_true(:)) / sqrt(2*numel(ym_true));

yp = yp_true + sigma_p_k * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
ym = ym_true + sigma_m_k * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
printm('snr_p = %g dB', 20*log(norm(yp_true(:)) / norm(col(yp_true-yp))));
printm('snr_m = %g dB', 20*log(norm(ym_true(:)) / norm(col(ym_true-ym))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
yp_im = abs(ifft2(yp));
ym_im = abs(ifft2(ym));

%% Maximum-Likelihood Estimation
% Dictionary-based estimation via variable projection method
weights = ones(M, 1);
W = spdiags(weights, 0, M, M);    % Weighting matrix
y = reshape(permute(cat(3, ys_im, yp_im, ym_im), [3 1 2]), [M nx*ny]);

maxProd = zeros(1, nx*ny);
idx = zeros(1, nx*ny);
for k = 1:K
    % Compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;
    
    % If the kth inner produt is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    idx(update) = k;
end

% Extract indices for maximum-likelihood maps
[t1_idx, t2_idx, s_idx] = ind2sub([length(T1) length(T2) length(kap)], idx);
T1_ml = reshape(T1(t1_idx), [nx ny]); 
T2_ml = reshape(T2(t2_idx), [nx ny]);
kap_ml = reshape(kap(s_idx), [nx ny]);

% % M0s initial guess
% M0s_ml = NaN(nx*ny, 1);
% for q = 1:length(idx)
%     M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
% end
% M0s_ml = reshape(M0s_ml, [nx ny]);

% For now, take true value wf = 0;
wf_ml = zeros(nx, ny);

%% Preprocessing and Masking
% Project Images to within range
T1max = 5000;       
T1min = 5;          
T2max = 500;    
T2min = 5;     

T1_ml = max(T1_ml, T1min);
T1_ml = min(T1_ml, T1max);
T2_ml = max(T2_ml, T2min);
T2_ml = min(T2_ml, T2max); 

% Preprocessing and Masking
tight_mask = T2s_msk;
loose_mask = imdilate(tight_mask, strel('disk', 5)); 

% Set voxels inside loose mask but outside tight mask to mean
T1_ml(~tight_mask & loose_mask)  = mean(col(T1_ml));   
T2_ml(~tight_mask & loose_mask)  = mean(col(T2_ml));   
kap_ml(~tight_mask & loose_mask) = mean(col(kap_ml));
wf_ml(~tight_mask & loose_mask)  = mean(col(wf_ml));   

% Set voxels outside both tight and loose masks to zero
T1_ml(~loose_mask)  = 0;  
T2_ml(~loose_mask)  = 0;  
kap_ml(~loose_mask) = 0;
wf_ml(~loose_mask)  = 0;  

% Median filtering (M0s_med might not be right...)
T1_med = medfilt2(T1_ml);  
T2_med = medfilt2(T2_ml);  
kap_med = medfilt2(kap_ml);

%% Regularized, Joint M0s,T1,T2,kap Estimation
% Define iteration parameters
n_outer = 10;    
niterM = 100;   
niter1 = 100;   
niter2 = 100;   
niterk = 50;

tolM = 10^-6;   
tol1 = 10^-7;  
tol2 = 10^-7;  
tolk = 10^-4; 

is_mag = 1;
disp = 1; 

% Define regularizers, Rm, R1, and R2
betaM = 2^-4; 
beta1 = 2^-7; 
beta2 = 2^-10;
betaK = 2^-3;

deltaM = 10^-2; 
delta1 = 10^-1;
delta2 = 10^-2;

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);
Rk = Reg1(loose_mask, 'pot_arg', {'quad'}, 'beta', betaK);

% Reshape data inputs to work with 3D implementation
ys_3D = permute(ys_im, [1 2 4 3]);
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);

% Regularized Reconstruction
tic; 
[M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
    mri_spgr_dess_m0st1t2kap_map(T1_ml, T2_ml, kap_ml, wf_ml, flip_s,...
    flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
    TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
    niter2, niterk, tolM, tol1, tol2, tolk, is_mag, disp);
recon_time = toc; 

% % Sanity check: note that a continuous kap_init needs to be used...
% [M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
%     mri_spgr_dess_m0st1t2kap_map(T1_true, T2_true, kap_ml, wf_ml, flip_s,...
%     flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
%     TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
%     niter2, niterK, tolM, tol1, tol2, tolk, is_mag, disp);

% Postprocessing for display
% For real data, can use this mask
reg_mask = abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:)));

% % Remove pixels outside reg_mask
% M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  kap_ml(~reg_mask) = 0;
% M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; kap_med(~reg_mask) = 0;
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0; kap_reg(~reg_mask) = 0;

%% Images, RMSE, and Comparisons
% Noise scaling of display
eps_scl = 0.2;

% M0s Images 
figure; im('notick', fliplr(abs(M0s_ml)), [0 1], 'cbar', ' ');
% print -deps M0s_ml_both.eps;
figure; im('notick', fliplr(abs(M0s_med)), [0 1], 'cbar', ' ');
% print -deps M0s_med_both.eps;
figure; im('notick', fliplr(abs(M0s_reg)), [0 1], 'cbar', ' ');
% print -deps M0s_reg_both.eps;
figure; im('notick', fliplr(abs(M0s_true)), [0 1], 'cbar', ' ');
% print -deps M0s_tru.eps;

figure; im('notick', fliplr(abs(M0s_ml-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_ml_err_both.eps;
figure; im('notick', fliplr(abs(M0s_med-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_med_err_both.eps;
figure; im('notick', fliplr(abs(M0s_reg-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_reg_err_both.eps;

% T1 Images
T1_med(~loose_mask) = 0;
figure; im('notick', fliplr(T1_ml),  [0 2000], 'cbar', ' ');
% print -deps T1_ml_both.eps;
figure; im('notick', fliplr(T1_med), [0 2000], 'cbar', ' ');
% print -deps T1_med_both.eps;
figure; im('notick', fliplr(T1_reg), [0 2000], 'cbar', ' ');
% print -deps T1_reg_both.eps;
figure; im('notick', fliplr(T1_true), [0 2000], 'cbar', ' ');
% print -deps T1_tru.eps;

figure; im('notick', fliplr(abs(T1_ml-T1_true)),  [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_ml_err_both.eps;
figure; im('notick', fliplr(abs(T1_med-T1_true)), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_med_err_both.eps;
figure; im('notick', fliplr(abs(T1_reg-T1_true)), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_reg_err_both.eps;

% T2 Images
figure; im('notick', fliplr(T2_ml),  [0 200], 'cbar', ' ');
% print -deps T2_ml_both.eps;
figure; im('notick', fliplr(T2_med), [0 200], 'cbar', ' ');
% print -deps T2_med_both.eps;
figure; im('notick', fliplr(T2_reg), [0 200], 'cbar', ' ');
% print -deps T2_reg_both.eps;
figure; im('notick', fliplr(T2_true), [0 200], 'cbar', ' ');
% print -deps T2_tru.eps;

figure; im('notick', fliplr(abs(T2_ml-T2_true)),  [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_ml_err_both.eps;
figure; im('notick', fliplr(abs(T2_med-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_med_err_both.eps;
figure; im('notick', fliplr(abs(T2_reg-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_reg_err_both.eps;

% kap Images
figure; im('notick', fliplr(kap_ml),  [0 2], 'cbar', ' ');
% print -deps kap_ml_both.eps;
figure; im('notick', fliplr(kap_med), [0 2], 'cbar', ' ');
% print -deps kap_med_both.eps;
figure; im('notick', fliplr(kap_reg), [0 2], 'cbar', ' ');
% print -deps kap_reg_both.eps;
figure; im('notick', fliplr(kap_true), [0 2], 'cbar', ' ');
% print -deps kap_tru_both.eps;

figure; im('notick', fliplr(abs(kap_ml-kap_true)),  [0 eps_scl*2], 'cbar', ' ');
% print -deps kap_ml_err_both.eps;
figure; im('notick', fliplr(abs(kap_med-kap_true)), [0 eps_scl*2], 'cbar', ' ');
% print -deps kap_med_err_both.eps;
figure; im('notick', fliplr(abs(kap_reg-kap_true)), [0 eps_scl*2], 'cbar', ' ');
% print -deps kap_reg_err_both.eps;

% Cost vs. Iteration
figure; hold on;
scatter(1:4:4*n_outer, cost(2:4:end), 'bo');
scatter(2:4:4*n_outer, cost(3:4:end), 'ro');
scatter(3:4:4*n_outer, cost(4:4:end), 'go');
scatter(4:4:4*n_outer, cost(5:4:end), 'mo');
plot(0:4*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M0s update', 'T1 update', 'T2 update', 'kap update');
% print -depsc cost_vs_iter.eps

% Compute NRMSE and RMSE
nrmse = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ sum(tru(:).^2));
rmse  = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ numel(tru));

% GM and WM masks for task-based performance assessment
mask_gm = labels(:,:,slice) == 2;
mask_wm = labels(:,:,slice) == 3;

% M0s NRMSE values
nrmse_m0s_ml_gm     = nrmse(M0s_ml(mask_gm),  M0s_true(mask_gm));
nrmse_m0s_med_gm    = nrmse(M0s_med(mask_gm), M0s_true(mask_gm));
nrmse_m0s_reg_gm    = nrmse(M0s_reg(mask_gm), M0s_true(mask_gm));
nrmse_m0s_tru_gm    = nrmse(M0s_true(mask_gm), M0s_true(mask_gm));

nrmse_m0s_ml_wm     = nrmse(M0s_ml(mask_wm),  M0s_true(mask_wm));
nrmse_m0s_med_wm    = nrmse(M0s_med(mask_wm), M0s_true(mask_wm));
nrmse_m0s_reg_wm    = nrmse(M0s_reg(mask_wm), M0s_true(mask_wm));
nrmse_m0s_tru_wm    = nrmse(M0s_true(mask_wm), M0s_true(mask_wm));

% T1 RMSE values
rmse_t1_ml_gm       = rmse(T1_ml(mask_gm),  T1_true(mask_gm));
rmse_t1_med_gm      = rmse(T1_med(mask_gm), T1_true(mask_gm));
rmse_t1_reg_gm      = rmse(T1_reg(mask_gm), T1_true(mask_gm));
rmse_t1_tru_gm      = rmse(T1_true(mask_gm), T1_true(mask_gm));

rmse_t1_ml_wm       = rmse(T1_ml(mask_wm),  T1_true(mask_wm));
rmse_t1_med_wm      = rmse(T1_med(mask_wm), T1_true(mask_wm));
rmse_t1_reg_wm      = rmse(T1_reg(mask_wm), T1_true(mask_wm));
rmse_t1_tru_wm      = rmse(T1_true(mask_wm), T1_true(mask_wm));

% T2 RMSE values
rmse_t2_ml_gm       = rmse(T2_ml(mask_gm),  T2_true(mask_gm));
rmse_t2_med_gm      = rmse(T2_med(mask_gm), T2_true(mask_gm));
rmse_t2_reg_gm      = rmse(T2_reg(mask_gm), T2_true(mask_gm));
rmse_t2_tru_gm      = rmse(T2_true(mask_gm), T2_true(mask_gm));

rmse_t2_ml_wm       = rmse(T2_ml(mask_wm),  T2_true(mask_wm));
rmse_t2_med_wm      = rmse(T2_med(mask_wm), T2_true(mask_wm));
rmse_t2_reg_wm      = rmse(T2_reg(mask_wm), T2_true(mask_wm));
rmse_t2_tru_wm      = rmse(T2_true(mask_wm), T2_true(mask_wm));

% kap NRMSE values
nrmse_kap_ml_gm     = nrmse(kap_ml(mask_gm),  kap_true(mask_gm));
nrmse_kap_med_gm    = nrmse(kap_med(mask_gm), kap_true(mask_gm));
nrmse_kap_reg_gm    = nrmse(kap_reg(mask_gm), kap_true(mask_gm));
nrmse_kap_tru_gm    = nrmse(kap_true(mask_gm), kap_true(mask_gm));

nrmse_kap_ml_wm     = nrmse(kap_ml(mask_wm),  kap_true(mask_wm));
nrmse_kap_med_wm    = nrmse(kap_med(mask_wm), kap_true(mask_wm));
nrmse_kap_reg_wm    = nrmse(kap_reg(mask_wm), kap_true(mask_wm));
nrmse_kap_tru_wm    = nrmse(kap_true(mask_wm), kap_true(mask_wm));

% M0s Means and Standard Deviations
m0s_ml_gm           = [mean(M0s_ml(mask_gm)),   std(M0s_ml(mask_gm))];
m0s_med_gm          = [mean(M0s_med(mask_gm)),  std(M0s_med(mask_gm))];
m0s_reg_gm          = [mean(M0s_reg(mask_gm)),  std(M0s_reg(mask_gm))];
m0s_tru_gm          = [mean(M0s_true(mask_gm)), std(M0s_true(mask_gm))];

m0s_ml_wm           = [mean(M0s_ml(mask_wm)),   std(M0s_ml(mask_wm))];
m0s_med_wm          = [mean(M0s_med(mask_wm)),  std(M0s_med(mask_wm))];
m0s_reg_wm          = [mean(M0s_reg(mask_wm)),  std(M0s_reg(mask_wm))];
m0s_tru_wm          = [mean(M0s_true(mask_wm)), std(M0s_true(mask_wm))];

% T1 Means and Standard Deviations
t1_ml_gm            = [mean(T1_ml(mask_gm)),   std(T1_ml(mask_gm))];
t1_med_gm           = [mean(T1_med(mask_gm)),  std(T1_med(mask_gm))];
t1_reg_gm           = [mean(T1_reg(mask_gm)),  std(T1_reg(mask_gm))];
t1_tru_gm           = [mean(T1_true(mask_gm)), std(T1_true(mask_gm))];

t1_ml_wm            = [mean(T1_ml(mask_wm)),   std(T1_ml(mask_wm))];
t1_med_wm           = [mean(T1_med(mask_wm)),  std(T1_med(mask_wm))];
t1_reg_wm           = [mean(T1_reg(mask_wm)),  std(T1_reg(mask_wm))];
t1_tru_wm           = [mean(T1_true(mask_wm)), std(T1_true(mask_wm))];

% T2 Means and Standard Deviations
t2_ml_gm            = [mean(T2_ml(mask_gm)),   std(T2_ml(mask_gm))];
t2_med_gm           = [mean(T2_med(mask_gm)),  std(T2_med(mask_gm))];
t2_reg_gm           = [mean(T2_reg(mask_gm)),  std(T2_reg(mask_gm))];
t2_tru_gm           = [mean(T2_true(mask_gm)), std(T2_true(mask_gm))];

t2_ml_wm            = [mean(T2_ml(mask_wm)),   std(T2_ml(mask_wm))];
t2_med_wm           = [mean(T2_med(mask_wm)),  std(T2_med(mask_wm))];
t2_reg_wm           = [mean(T2_reg(mask_wm)),  std(T2_reg(mask_wm))];
t2_tru_wm           = [mean(T2_true(mask_wm)), std(T2_true(mask_wm))];

% Kappa means and Standard Deviations
kap_ml_gm           = [mean(kap_ml(mask_gm)),  std(kap_ml(mask_gm))];
kap_med_gm          = [mean(kap_med(mask_gm)), std(kap_med(mask_gm))];
kap_reg_gm          = [mean(kap_reg(mask_gm)), std(kap_reg(mask_gm))];
kap_tru_gm          = [mean(kap_true(mask_gm)), std(kap_true(mask_gm))];

kap_ml_wm           = [mean(kap_ml(mask_wm)),  std(kap_ml(mask_wm))];
kap_med_wm          = [mean(kap_med(mask_wm)), std(kap_med(mask_wm))];
kap_reg_wm          = [mean(kap_reg(mask_wm)), std(kap_reg(mask_wm))];
kap_tru_wm          = [mean(kap_true(mask_wm)), std(kap_true(mask_wm))];

% Print Statistics
fprintf('\n\nM0s Statistics\n');
fprintf('\t\tGM Mean\t\tGM SD\t\tGM NRMSE\tWM Mean\t\tWM SD\t\tWM NRMSE\n');
fprintf('ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    m0s_ml_gm(1), m0s_ml_gm(2), nrmse_m0s_ml_gm,...
    m0s_ml_wm(1), m0s_ml_wm(2), nrmse_m0s_ml_wm);
fprintf('Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    m0s_med_gm(1), m0s_med_gm(2), nrmse_m0s_med_gm,...
    m0s_med_wm(1), m0s_med_wm(2), nrmse_m0s_med_wm);
fprintf('Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    m0s_reg_gm(1), m0s_reg_gm(2), nrmse_m0s_reg_gm,...
    m0s_reg_wm(1), m0s_reg_wm(2), nrmse_m0s_reg_wm);
fprintf('True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
    m0s_tru_gm(1), m0s_tru_gm(2), nrmse_m0s_tru_gm,...
    m0s_tru_wm(1), m0s_tru_wm(2), nrmse_m0s_tru_wm);

fprintf('T1 Statistics\n');
fprintf('\t\tGM Mean\t\tGM SD\t\tGM NRMSE\tWM Mean\t\tWM SD\t\tWM NRMSE\n');
fprintf('ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t1_ml_gm(1), t1_ml_gm(2), rmse_t1_ml_gm,...
    t1_ml_wm(1), t1_ml_wm(2), rmse_t1_ml_wm);
fprintf('Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t1_med_gm(1), t1_med_gm(2), rmse_t1_med_gm,...
    t1_med_wm(1), t1_med_wm(2), rmse_t1_med_wm);
fprintf('Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t1_reg_gm(1), t1_reg_gm(2), rmse_t1_reg_gm,...
    t1_reg_wm(1), t1_reg_wm(2), rmse_t1_reg_wm);
fprintf('True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
    t1_tru_gm(1), t1_tru_gm(2), rmse_t1_tru_gm,...
    t1_tru_wm(1), t1_tru_wm(2), rmse_t1_tru_wm);

fprintf('T2 Statistics\n');
fprintf('\t\tGM Mean\t\tGM SD\t\tGM NRMSE\tWM Mean\t\tWM SD\t\tWM NRMSE\n');
fprintf('ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t2_ml_gm(1), t2_ml_gm(2), rmse_t2_ml_gm,...
    t2_ml_wm(1), t2_ml_wm(2), rmse_t2_ml_wm);
fprintf('Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t2_med_gm(1), t2_med_gm(2), rmse_t2_med_gm,...
    t2_med_wm(1), t2_med_wm(2), rmse_t2_med_wm);
fprintf('Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    t2_reg_gm(1), t2_reg_gm(2), rmse_t2_reg_gm,...
    t2_reg_wm(1), t2_reg_wm(2), rmse_t2_reg_wm);
fprintf('True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
    t2_tru_gm(1), t2_tru_gm(2), rmse_t2_tru_gm,...
    t2_tru_wm(1), t2_tru_wm(2), rmse_t2_tru_wm);

% Print Statistics
fprintf('kappa Statistics\n');
fprintf('\t\tGM Mean\t\tGM SD\t\tGM NRMSE\tWM Mean\t\tWM SD\t\tWM NRMSE\n');
fprintf('ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    kap_ml_gm(1), kap_ml_gm(2), nrmse_kap_ml_gm,...
    kap_ml_wm(1), kap_ml_wm(2), nrmse_kap_ml_wm);
fprintf('Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    kap_med_gm(1), kap_med_gm(2), nrmse_kap_med_gm,...
    kap_med_wm(1), kap_med_wm(2), nrmse_kap_med_wm);
fprintf('Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    kap_reg_gm(1), kap_reg_gm(2), nrmse_kap_reg_gm,...
    kap_reg_wm(1), kap_reg_wm(2), nrmse_kap_reg_wm);
fprintf('True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
    kap_tru_gm(1), kap_tru_gm(2), nrmse_kap_tru_gm,...
    kap_tru_wm(1), kap_tru_wm(2), nrmse_kap_tru_wm);