%% Header Script for Regularized M0s/T1/T2 Recon from Synthetic SPGR/DESS Data
% For comparison of optimized (2,1), (0,2), (1,1) SPGR/DESS profiles
% Dictionary initialization method and joint estimation method constant
%
% Written by: Gopal Nataraj
% Originally created 2014-06-17
% Last modified 2015-11-11
%
% Mild flip-angle variation included
% todo: coil combination initial step (body coil data assumed)

%% Load true parameter maps and set scan protocol
% Load digital phantom
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../matlab;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;
addpath('Calculations/SPGR/');
addpath('Calculations/DESS/');
addpath('Helper/');
addpath('CRLB_Analysis/');

% Output desired? 
pr = 1;

% Get true values
b0 = 3.0;                               % Tesla
slice = 81; [nx, ny, nz] = size(labels);
M0_true = double(squeeze(labels(:,:,slice)));
T1_true = double(squeeze(labels(:,:,slice)));
T2_true = double(squeeze(labels(:,:,slice)));
T2s_true = double(squeeze(labels(:,:,slice)));
for idx = 0:10
    f = mri_brainweb_params(idx, 'b0', b0);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1; 
    T2_true(T2_true == idx) = f.t2;
    T2s_true(T2s_true == idx) = f.t2s;
end
wf_true = 0; 

% Common imaging parameters (fixed for all scans)
TE_global = 4.67;                       % ms

% SPGR imaging parameters
% TRs = [12.2 12.2]';                     % ms (2,1) Profile
% TRs = [13.9]';                          % ms (1,1) Profile
TRs = []';                              % ms (0,2) Profile

% flip_s = [15 5]' * pi/180;              % radians (2,1) Profile
% flip_s = [15]' * pi/180;                % radians (1,1) Profile
flip_s = []' * pi/180;                  % radians (0,2) Profile

nfs = length(flip_s);   
TEs = TE_global * ones(nfs, 1);         % ms

% DESS imaging parameters
% TRd = [17.5]';                          % ms (2,1) Profile
% TRd = [28.0]';                          % ms (1,1) Profile
TRd = [24.4 17.5]';                     % ms (0,2) Profile

% flip_d = [30]' * pi/180;                % radians (2,1) Profile
% flip_d = [10]' * pi/180;                % radians (1,1) Profile
flip_d = [35 10]' * pi/180;             % radians (0,2) Profile

nfd = length(flip_d);          
TEd = TE_global * ones(nfd, 1);         % ms

% Total number of datasets 
M = length(flip_s) + 2*length(flip_d);

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = M0_true ~= 0; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TE_global ./ T2s_true(T2s_msk));

% Create a true excitation profile map
% Cannot allow too much variation in flip because then ML estimate biased
perc_var = 0.0;                         % allow +/-0% variation in flip
[xx, yy] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
kap_true = (1+perc_var/2) - perc_var*(xx.^2 + yy.^2); 
kap_true(~T2s_msk) = 0;                 % No excitation where M0 = 0;


%% Dictionary creation (could be precomputed)
% Estimate a constant flip-angle scaling factor to create dictionary
kap_scale = mean(kap_true(T2s_msk));    % Needs to be approximately uniform

T1 = logspace(2, 3.5, 1000);
T2 = logspace(1, 2.5, 1000);
D = NaN(M, length(T1), length(T2));
K = length(T1) * length(T2);

for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        % SPGR dictionary component
        for a = 1:nfs
            D(a, t1, t2) = spgr_fun(1, T1(t1), kap_scale,...
                flip_s(a), TRs(a), TEs(a), wf_true, 1);
        end

        % DESS dictionary component
        for a = 1:nfd
            [D(nfs+a, t1, t2), D(nfs+nfd+a, t1, t2)] = ...
                dess_fun(1, T1(t1), T2(t2), kap_scale,...
                flip_d(a), TRd(a), TEd(a), wf_true, 1);
        end
    end
end
D = reshape(D, [M, K]);


%% Synthesize Phantom SPGR Data
% Forward model: make true data
ys_k_true = NaN(nx, ny, nfs);
ys_im_true = NaN(nx, ny, nfs);
for a = 1:nfs
    ys_im_true(:,:,a) = spgr_fun(M0s_true, T1_true, kap_true,...
        flip_s(a), TRs(a), TEs(a), wf_true, 1);
    ys_k_true(:,:,a) = fft2(ys_im_true(:,:,a));
end

% Rescale measured noise variance to current voxel size
oldVoxelVol = 1 * 1 * 5;                        % mm^3
newVoxelVol = (240/256) * (240/256) * (30/6);   % mm^3
oldvar_im = 1.31e-7;                            
newvar_im = oldvar_im * (oldVoxelVol / newVoxelVol);

% Add complex white gaussian noise
var_s_im = newvar_im; sigma_s_k = sqrt(var_s_im * (nx*ny));
ys = ys_k_true + sigma_s_k * (randn(size(ys_k_true)) + 1i * randn(size(ys_k_true)));
printm('snr_s = %g dB', 20*log10(norm(ys_k_true(:)) / norm(col(ys_k_true-ys))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
ys_im = abs(ifft2(ys));


%% Synthesize Phantom DESS Data
yp_k_true = NaN(nx, ny, nfd);
ym_k_true = NaN(nx, ny, nfd);
yp_im_true = NaN(nx, ny, nfd);
ym_im_true = NaN(nx, ny, nfd);
for a = 1:nfd
    [yp_im_true(:,:,a), ym_im_true(:,:,a)] = dess_fun(M0s_true, T1_true, T2_true, kap_true,...
        flip_d(a), TRd(a), TEd(a), wf_true, 1);
    yp_k_true(:,:,a) = fft2(yp_im_true(:,:,a));
    ym_k_true(:,:,a) = fft2(ym_im_true(:,:,a));
end

% Add complex white gaussian noise
var_p_im = newvar_im; sigma_p_k = sqrt(var_p_im * (nx*ny));
var_m_im = newvar_im; sigma_m_k = sqrt(var_m_im * (nx*ny));
yp = yp_k_true + sigma_p_k * (randn(size(yp_k_true)) + 1i * randn(size(yp_k_true)));
ym = ym_k_true + sigma_m_k * (randn(size(ym_k_true)) + 1i * randn(size(ym_k_true)));
printm('snr_p = %g dB', 20*log10(norm(yp_k_true(:)) / norm(col(yp_k_true-yp))));
printm('snr_m = %g dB', 20*log10(norm(ym_k_true(:)) / norm(col(ym_k_true-ym))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
yp_im = abs(ifft2(yp));
ym_im = abs(ifft2(ym));

% For paper, compute SNR for full (image-domain) dataset
Y_im = cat(3, ys_im, yp_im, ym_im);
Y_im_true = cat(3, ys_im_true, yp_im_true, ym_im_true);
[PSNR, SNR] = psnr(Y_im, Y_im_true);
printm('Image domain PSNR = %0.3f dB and SNR = %0.3f dB', PSNR, SNR);


%% Maximum-Likelihood Estimation
% Dictionary-based estimation via variable projection method
tic;
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
[t1_idx, t2_idx] = ind2sub([length(T1) length(T2)], idx);
T1_ml = reshape(T1(t1_idx), [nx ny]); 
T2_ml = reshape(T2(t2_idx), [nx ny]);
time_ml = toc;

% % M0s initial guess
% M0s_ml = NaN(nx*ny, 1);
% for q = 1:length(idx)
%     M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
% end
% M0s_ml = reshape(M0s_ml, [nx ny]);

% For now, take true value wf=0
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
wf_ml(~tight_mask & loose_mask)  = mean(col(wf_ml));  

% Set voxels outside both tight and loose masks to zero
T1_ml(~loose_mask) = 0;  
T2_ml(~loose_mask) = 0;  
wf_ml(~loose_mask) = 0;  

% Median filtering (M0s_med might not be right...)
T1_med = medfilt2(T1_ml);  
T2_med = medfilt2(T2_ml);  


%% Regularized, Joint M0s,T1,T2 Estimation 
% REGULARIZED ESTIMATES NOT INCLUDED IN PAPER
% Define iteration parameters
n_outer = 10;    
niterM = 100;   
niter1 = 100;   
niter2 = 100;  
niterk = 0;         % kap not updated

tolM = 10^-6;   
tol1 = 10^-7;  
tol2 = 10^-7;
tolk = 10^-4;       % irrelevant

is_mag = 1;
disp = 0; 

% Define regularizers, Rm, R1, and R2
betaM = 2^-16; 
beta1 = 2^-11;      % 2^-7
beta2 = 2^-14;      % 2^-10
betaK = 2^-3;       % irrelevant

deltaM = 10^-4;
delta1 = 10^-1;     
delta2 = 10^-2;     

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);
Rk = Reg1(loose_mask, 'pot_arg', {'quad'}, 'beta', betaK);

% Reshape inputs to work with 3-D implementation
ys_3D = permute(ys_im, [1 2 4 3]);
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);

% % Method One: Sequential Estimation (old code, needs edit)
% % Need a M0s initial guess
% M0s_ml = NaN(nx*ny, 1);
% for q = 1:length(idx)
%     M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
% end
% M0s_ml = reshape(M0s_ml, [nx ny]);
% M0s_med = M0s_ml;
% 
% % M0s, T1 regularized reconstruction from SPGR
% [M0s_reg, T1_reg, ~] = mri_despot1_m0t1map(M0s_ml, E1s_med,...
%     ys_3D, flip_s_3D, loose_mask, Rm, R1, T1max, TR_s, n_outer,...
%     niterM, niter1, tolM, tol1, disp);
% 
% % T2 regularized reconstruction from DESS
% E1d_reg = exp(-TR_d ./ T1_reg);
% [M0s_reg, T1_reg, T2_reg, wf_reg, ~] = ...
%     mri_dess_t2map_nest83_M0s(M0s_reg, E1d_reg, E2d_med, wf_ml, ...
%     flip_d_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
%     TR_d, TE, R2, n_outer, niter2, tol2, disp);

% Method Two: Joint Estimation
% Assume we know the ground truth flip angle scaling
tic; 
[M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
    mri_spgr_dess_m0st1t2kap_map(T1_ml, T2_ml, kap_true, wf_ml, flip_s,...
    flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
    TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
    niter2, niterk, tolM, tol1, tol2, tolk, is_mag, disp);
time_reg = toc; 

% Postprocessing for display
% For real data, can use this mask
reg_mask = abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:)));

% Remove pixels outside reg_mask
M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;
M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0;
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0;

%% Images, RMSE, and Comparisons
% Noise scaling of display
eps_scl = 0.2;
suffix = sprintf('%uSPGR%uDESS', nfs, nfd);

% M0s Images 
m0smin = 0; m0smax = 1; m0s_rng = [m0smin m0smax];
figure; im('notick', fliplr(abs(M0s_ml)), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_ml_', suffix)), end;
figure; im('notick', fliplr(abs(M0s_med)), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_med_', suffix)), end;
figure; im('notick', fliplr(abs(M0s_reg)), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_reg_', suffix)), end;
figure; im('notick', fliplr(abs(M0s_true)), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_tru_', suffix)), end;

figure; im('notick', fliplr(abs(M0s_ml-M0s_true)), eps_scl*m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_ml_err_', suffix)), end;
figure; im('notick', fliplr(abs(M0s_med-M0s_true)), eps_scl*m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_med_err_', suffix)), end;
figure; im('notick', fliplr(abs(M0s_reg-M0s_true)), eps_scl*m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_reg_err_', suffix)), end;

% T1 Images
T1_med(~loose_mask) = 0;
t1min = 0; t1max = 2000; t1_rng = [t1min t1max];
figure; im('notick', fliplr(T1_ml), t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_ml_', suffix)), end;
figure; im('notick', fliplr(T1_med), t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_med_', suffix)), end;
figure; im('notick', fliplr(T1_reg), t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_reg_', suffix)), end;
figure; im('notick', fliplr(T1_true), t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_tru_', suffix)), end;

figure; im('notick', fliplr(abs(T1_ml-T1_true)), eps_scl*t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_ml_err_', suffix)), end;
figure; im('notick', fliplr(abs(T1_med-T1_true)), eps_scl*t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_med_err_', suffix)), end;
figure; im('notick', fliplr(abs(T1_reg-T1_true)), eps_scl*t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_reg_err_', suffix)), end;

% T2 Images
t2min = 0; t2max = 200; t2_rng = [t2min t2max];
figure; im('notick', fliplr(T2_ml), t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_ml_', suffix)), end;
figure; im('notick', fliplr(T2_med), t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_med_', suffix)), end;
figure; im('notick', fliplr(T2_reg), t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_reg_', suffix)), end;
figure; im('notick', fliplr(T2_true), t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_tru_', suffix)), end;

figure; im('notick', fliplr(abs(T2_ml-T2_true)), eps_scl*t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_ml_err_', suffix)), end;
figure; im('notick', fliplr(abs(T2_med-T2_true)), eps_scl*t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_med_err_', suffix)), end;
figure; im('notick', fliplr(abs(T2_reg-T2_true)), eps_scl*t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_reg_err_', suffix)), end;

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


%% Task-based performance assessment: bias, variance, (N)RMSE
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


% Print Statistics (previous contents overwritten!)
if (pr)
    fid = fopen(strcat('DigiPhant_3parm_summary_', suffix), 'w');
    fprintf(fid, strcat('Estimation Statistics for (', suffix, ') Profile'));
    fprintf(fid, sprintf('\n\tML run time: %0.2f seconds', time_ml));
    fprintf(fid, sprintf('\n\tReg run time: %0.2f seconds', time_reg));
    fprintf(fid, '\n\nM0s Statistics\n');
    fprintf(fid, '\t\tGM Mean\t\tGM SD\t\tGM NRMSE\tWM Mean\t\tWM SD\t\tWM NRMSE\n');
    fprintf(fid, 'ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        m0s_ml_gm(1), m0s_ml_gm(2), nrmse_m0s_ml_gm,...
        m0s_ml_wm(1), m0s_ml_wm(2), nrmse_m0s_ml_wm);
    % fprintf(fid, 'Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    %     m0s_med_gm(1), m0s_med_gm(2), nrmse_m0s_med_gm,...
    %     m0s_med_wm(1), m0s_med_wm(2), nrmse_m0s_med_wm);
    fprintf(fid, 'Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        m0s_reg_gm(1), m0s_reg_gm(2), nrmse_m0s_reg_gm,...
        m0s_reg_wm(1), m0s_reg_wm(2), nrmse_m0s_reg_wm);
    fprintf(fid, 'True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
        m0s_tru_gm(1), m0s_tru_gm(2), nrmse_m0s_tru_gm,...
        m0s_tru_wm(1), m0s_tru_wm(2), nrmse_m0s_tru_wm);

    fprintf(fid, 'T1 Statistics\n');
    fprintf(fid, '\t\tGM Mean\t\tGM SD\t\tGM RMSE\t\tWM Mean\t\tWM SD\t\tWM RMSE\n');
    fprintf(fid, 'ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        t1_ml_gm(1), t1_ml_gm(2), rmse_t1_ml_gm,...
        t1_ml_wm(1), t1_ml_wm(2), rmse_t1_ml_wm);
    % fprintf(fid, 'Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    %     t1_med_gm(1), t1_med_gm(2), rmse_t1_med_gm,...
    %     t1_med_wm(1), t1_med_wm(2), rmse_t1_med_wm);
    fprintf(fid, 'Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        t1_reg_gm(1), t1_reg_gm(2), rmse_t1_reg_gm,...
        t1_reg_wm(1), t1_reg_wm(2), rmse_t1_reg_wm);
    fprintf(fid, 'True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
        t1_tru_gm(1), t1_tru_gm(2), rmse_t1_tru_gm,...
        t1_tru_wm(1), t1_tru_wm(2), rmse_t1_tru_wm);

    fprintf(fid, 'T2 Statistics\n');
    fprintf(fid, '\t\tGM Mean\t\tGM SD\t\tGM RMSE\t\tWM Mean\t\tWM SD\t\tWM RMSE\n');
    fprintf(fid, 'ML:\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        t2_ml_gm(1), t2_ml_gm(2), rmse_t2_ml_gm,...
        t2_ml_wm(1), t2_ml_wm(2), rmse_t2_ml_wm);
    % fprintf(fid, 'Med:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
    %     t2_med_gm(1), t2_med_gm(2), rmse_t2_med_gm,...
    %     t2_med_wm(1), t2_med_wm(2), rmse_t2_med_wm);
    fprintf(fid, 'Reg:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n',...
        t2_reg_gm(1), t2_reg_gm(2), rmse_t2_reg_gm,...
        t2_reg_wm(1), t2_reg_wm(2), rmse_t2_reg_wm);
    fprintf(fid, 'True:\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\t\t%0.2f\n\n',...
        t2_tru_gm(1), t2_tru_gm(2), rmse_t2_tru_gm,...
        t2_tru_wm(1), t2_tru_wm(2), rmse_t2_tru_wm);
    fclose(fid);
end


%% Comparison of CRLB theory vs. empirical results
% Note: norm_crlb_dess_3parm() yields standard deviations for abs(M0s) = 1. 
% We must scale this CRLB to match the signal level of the data.

% Worst-case three parameter time-uncompensated CRLB
T1 = linspace(800, 1400, 5);        % WM/GM voxels of interest
T2 = linspace(50, 120, 5);          % WM/GM voxels of interest
kap = linspace(0.9, 1.1, 5);        % a.u. (10% variation)

% Noise covariance matrix 
Sigma_inv = diag(1./[var_s_im*ones(nfs,1); 
                     var_p_im*ones(nfd,1); 
                     var_m_im*ones(nfd,1)]);
% No time compensation
time_comp = 0;

worst_sigT1  = NaN(length(T1), length(T2), length(kap));
worst_sigT2  = NaN(length(T1), length(T2), length(kap));
for i_t1 = 1:length(T1)
for i_t2 = 1:length(T2)
for i_kap = 1:length(kap)
    [~, ~, worst_sigT1(i_t1,i_t2,i_kap), worst_sigT2(i_t1,i_t2,i_kap)]...
        = norm_crlb_dess_3parm(T1(i_t1), T2(i_t2), mean(wf_true(:)),...
        kap(i_kap)*flip_s, kap(i_kap)*flip_d, TRs, TRd, TE_global,...
        Sigma_inv, time_comp);
end
end
end
sigwt1_norm = max(worst_sigT1(:));
sigwt2_norm = max(worst_sigT2(:));

% WM/GM three-parameter time-uncompensated CRLB
[~, ~, sigt1_gm_norm, sigt2_gm_norm] = ...
    norm_crlb_dess_3parm(t1_tru_gm(1), t2_tru_gm(1), mean(wf_true(:)),...
    1*flip_s, 1*flip_d, TRs, TRd, TE_global, Sigma_inv, time_comp);
[~, ~, sigt1_wm_norm, sigt2_wm_norm] = ...
    norm_crlb_dess_3parm(t1_tru_wm(1), t2_tru_wm(1), mean(wf_true(:)),...
    1*flip_s, 1*flip_d, TRs, TRd, TE_global, Sigma_inv, time_comp);

% Scale the normalized CRLB to match signal level of the data
% Note M0s_theory_sd is the std. dev. of M0s*sqrt(SoS(sensitivities)).
% M0s_theory_sd / sqrt(SoS(sens)) would be LATENT std. dev. of M0s
% M0s_theory_norm = M0s_theory_sd / M0s_ml_mean is std. dev. of M0s = 1.
sigwt1_gm = sigwt1_norm / m0s_tru_gm(1);
sigwt1_wm = sigwt1_norm / m0s_tru_wm(1);
sigwt2_gm = sigwt2_norm / m0s_tru_gm(1);
sigwt2_wm = sigwt2_norm / m0s_tru_wm(1);

sigt1_gm  = sigt1_gm_norm / m0s_tru_gm(1);
sigt1_wm  = sigt1_wm_norm / m0s_tru_wm(1);
sigt2_gm  = sigt2_gm_norm / m0s_tru_gm(1);
sigt2_wm  = sigt2_wm_norm / m0s_tru_wm(1);

% Code for MATLAB 2014b
% Set T1 and T2 bin widths
bin_wdth_t1 = 10;
bin_wdth_t2 = 1;

% Histograms of T1 GM/WM ML estimates (will be discretized)
min_t1_gm = floor(min(T1_ml(mask_gm)) / bin_wdth_t1) * bin_wdth_t1 - bin_wdth_t1/2;
max_t1_gm = ceil(max(T1_ml(mask_gm))  / bin_wdth_t1) * bin_wdth_t1 - bin_wdth_t1/2;

X1_ml_gm = linspace(min_t1_gm, max_t1_gm, 1000);
Y1_ml_gm_crb = normpdf(X1_ml_gm, t1_tru_gm(1), sigwt1_gm);
Y1_ml_gm_tru = normpdf(X1_ml_gm, t1_tru_gm(1), sigt1_gm);

figure; hold on; cmp = get(gca, 'ColorOrder');...
    histogram(T1_ml(mask_gm), 'BinLimits', [min_t1_gm max_t1_gm], 'BinWidth',...
        bin_wdth_t1, 'Normalization', 'pdf', 'FaceAlpha', 0.4, 'FaceColor', cmp(3,:));...
    plot(X1_ml_gm, Y1_ml_gm_crb, 'LineWidth', 2, 'LineStyle', '-',  'Color', cmp(5,:));...
    plot(X1_ml_gm, Y1_ml_gm_tru, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmp(4,:));...
    hold off; axis tight;
    set(gca, 'XLim', [min_t1_gm max_t1_gm]);
    set(gca, 'XTick', [min_t1_gm + bin_wdth_t1/2 : 2*bin_wdth_t1 : max_t1_gm - bin_wdth_t1/2]);
    xlabel('', 'FontSize', 16); xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
	ylabel('', 'FontSize', 16); ytick([]);
    l1_gm = legend('Empirical', 'Worst-case SD', 'Latent SD in GM'); 
    set(l1_gm, 'FontSize', 16);
    if (pr), print('-depsc', strcat('T1_hist_ml_gm_3parm_', suffix, '.eps')), end;

    
min_t1_wm = floor(min(T1_ml(mask_wm)) / bin_wdth_t1) * bin_wdth_t1 - bin_wdth_t1/2;
max_t1_wm = ceil(max(T1_ml(mask_wm))  / bin_wdth_t1) * bin_wdth_t1 + bin_wdth_t1/2;    
    
X1_ml_wm = linspace(min_t1_wm, max_t1_wm, 1000);
Y1_ml_wm_crb = normpdf(X1_ml_wm, t1_tru_wm(1), sigwt1_wm);
Y1_ml_wm_tru = normpdf(X1_ml_wm, t1_tru_wm(1), sigt1_wm);

figure; hold on; cmp = get(gca, 'ColorOrder');...
    histogram(T1_ml(mask_wm), 'BinLimits', [min_t1_wm max_t1_wm], 'BinWidth',...
        bin_wdth_t1, 'Normalization', 'pdf', 'FaceAlpha', 0.4, 'FaceColor', cmp(3,:));...
    plot(X1_ml_wm, Y1_ml_wm_crb, 'LineWidth', 2, 'LineStyle', '-',  'Color', cmp(5,:));...
    plot(X1_ml_wm, Y1_ml_wm_tru, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmp(4,:));...
    hold off; axis tight;
    set(gca, 'XLim', [min_t1_wm max_t1_wm]);
    set(gca, 'XTick', [min_t1_wm + bin_wdth_t1/2 : 2*bin_wdth_t1 : max_t1_wm - bin_wdth_t1/2]);
    xlabel('', 'FontSize', 16); xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
	ylabel('', 'FontSize', 16); ytick([]);
    l1_wm = legend('Empirical', 'Worst-case SD', 'Latent SD in WM');
    set(l1_wm, 'FontSize', 16);
    if (pr), print('-depsc', strcat('T1_hist_ml_wm_3parm_', suffix, '.eps')), end;
    
    
% Histograms of T2 GM/WM ML estimates (will be discretized)
min_t2_gm = floor(min(T2_ml(mask_gm)) / bin_wdth_t2) * bin_wdth_t2 - bin_wdth_t2/2;
max_t2_gm = ceil(max(T2_ml(mask_gm))  / bin_wdth_t2) * bin_wdth_t2 + bin_wdth_t2/2;

X2_ml_gm = linspace(min_t2_gm, max_t2_gm, 1000);
Y2_ml_gm_crb = normpdf(X2_ml_gm, t2_tru_gm(1), sigwt2_gm);
Y2_ml_gm_tru = normpdf(X2_ml_gm, t2_tru_gm(1), sigt2_gm);

figure; hold on; cmp = get(gca, 'ColorOrder');...
    histogram(T2_ml(mask_gm), 'BinLimits', [min_t2_gm max_t2_gm], 'BinWidth',...
        bin_wdth_t2, 'Normalization', 'pdf', 'FaceAlpha', 0.4, 'FaceColor', cmp(3,:));...
    plot(X2_ml_gm, Y2_ml_gm_crb, 'LineWidth', 2, 'LineStyle', '-',  'Color', cmp(5,:));...
    plot(X2_ml_gm, Y2_ml_gm_tru, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmp(4,:));...
    hold off; axis tight;
    set(gca, 'XLim', [min_t2_gm max_t2_gm]);
    set(gca, 'XTick', [min_t2_gm + bin_wdth_t2/2 : 2*bin_wdth_t2 : max_t2_gm - bin_wdth_t2/2]);
    xlabel('', 'FontSize', 16); xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
	ylabel('', 'FontSize', 16); ytick([]);
    l1_gm = legend('Empirical', 'Worst-case SD', 'Latent SD in GM'); 
    set(l1_gm, 'FontSize', 16);
    if (pr), print('-depsc', strcat('T2_hist_ml_gm_3parm_', suffix, '.eps')), end;

    
min_t2_wm = floor(min(T2_ml(mask_wm)) / bin_wdth_t2) * bin_wdth_t2 + bin_wdth_t2/2;
max_t2_wm = ceil(max(T2_ml(mask_wm))  / bin_wdth_t2) * bin_wdth_t2 + bin_wdth_t2/2;    
    
X2_ml_wm = linspace(min_t2_wm, max_t2_wm, 1000);
Y2_ml_wm_crb = normpdf(X2_ml_wm, t2_tru_wm(1), sigwt2_wm);
Y2_ml_wm_tru = normpdf(X2_ml_wm, t2_tru_wm(1), sigt2_wm);

figure; hold on; cmp = get(gca, 'ColorOrder');...
    histogram(T2_ml(mask_wm), 'BinLimits', [min_t2_wm max_t2_wm], 'BinWidth',...
        bin_wdth_t2, 'Normalization', 'pdf', 'FaceAlpha', 0.4, 'FaceColor', cmp(3,:));...
    plot(X2_ml_wm, Y2_ml_wm_crb, 'LineWidth', 2, 'LineStyle', '-',  'Color', cmp(5,:));...
    plot(X2_ml_wm, Y2_ml_wm_tru, 'LineWidth', 2, 'LineStyle', '--', 'Color', cmp(4,:));...
    hold off; axis tight;
    set(gca, 'XLim', [min_t2_wm max_t2_wm]);
    set(gca, 'XTick', [min_t2_wm + bin_wdth_t2/2 : 2*bin_wdth_t2 : max_t2_wm - bin_wdth_t2/2]);
    xlabel('', 'FontSize', 16); xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
	ylabel('', 'FontSize', 16); ytick([]);
    l1_wm = legend('Empirical', 'Worst-case SD', 'Latent SD in WM');
    set(l1_wm, 'FontSize', 16);
    if (pr), print('-depsc', strcat('T2_hist_ml_wm_3parm_', suffix, '.eps')), end;
    
% % Code for MATLAB 2013a (old)
% % Histograms of T1/T2 GM ML estimates (will be discretized)
% X1_ml_gm = linspace(min(T1_ml(mask_gm)), max(T1_ml(mask_gm)), 1000);
% Y1_ml_gm_crb = normpdf(X1_ml_gm, t1_tru_gm(1), sigwt1_gm);
% Y1_ml_gm_tru = normpdf(X1_ml_gm, t1_tru_gm(1), sigt1_gm);
% [N1_ml_gm, B1_ml_gm] = hist(T1_ml(mask_gm), floor(sum(mask_gm(:))/500));...
% figure; bar(B1_ml_gm, N1_ml_gm/max(N1_ml_gm)*max(Y1_ml_gm_tru)); hold on;...
%     plot(X1_ml_gm, Y1_ml_gm_crb, 'r', 'LineWidth', 2);...
%     plot(X1_ml_gm, Y1_ml_gm_tru, 'g', 'LineWidth', 2); hold off;...
%     legend('Empirical', 'CRLB', 'Latent');
%     if (~pr), title('T1 ML Estimates in GM');
%     else print('-depsc', strcat('T1_hist_ml_gm_3parm_', suffix, '.eps')), end;
% 
% X2_ml_gm = linspace(min(T2_ml(mask_gm)), max(T2_ml(mask_gm)), 1000);...
% Y2_ml_gm_crb = normpdf(X2_ml_gm, t2_tru_gm(1), sigwt2_gm);
% Y2_ml_gm_tru = normpdf(X2_ml_gm, t2_tru_gm(1), sigt2_gm);
% [N2_ml_gm, B2_ml_gm] = hist(T2_ml(mask_gm), floor(sum(mask_gm(:))/500));
% figure; bar(B2_ml_gm, N2_ml_gm/max(N2_ml_gm)*max(Y2_ml_gm_tru)); hold on;...
%     plot(X2_ml_gm, Y2_ml_gm_crb, 'r', 'LineWidth', 2);...
%     plot(X2_ml_gm, Y2_ml_gm_tru, 'g', 'LineWidth', 2); hold off;...
%     legend('Empirical', 'CRLB', 'Latent');
%     if (~pr), title('T2 ML Estimates in GM');
%     else print('-depsc', strcat('T2_hist_ml_gm_3parm_', suffix, '.eps')), end;
% 
% % Histograms of T1/T2 WM ML estimates (will be discretized)
% X1_ml_wm = linspace(min(T1_ml(mask_wm)), max(T1_ml(mask_wm)), 1000);
% Y1_ml_wm_crb = normpdf(X1_ml_wm, t1_tru_wm(1), sigwt1_wm);
% Y1_ml_wm_tru = normpdf(X1_ml_wm, t1_tru_wm(1), sigt1_wm);
% [N1_ml_wm, B1_ml_wm] = hist(T1_ml(mask_wm), floor(sum(mask_wm(:))/500));...
% figure; bar(B1_ml_wm, N1_ml_wm/max(N1_ml_wm)*max(Y1_ml_wm_tru)); hold on;...
%     plot(X1_ml_wm, Y1_ml_wm_crb, 'r', 'LineWidth', 2);...
%     plot(X1_ml_wm, Y1_ml_wm_tru, 'g', 'LineWidth', 2); hold off;...
%     legend('Empirical', 'CRLB', 'Latent');
%     if (~pr), title('T1 ML Estimates in WM');
%     else print('-depsc', strcat('T1_hist_ml_wm_3parm_', suffix, '.eps')), end;
% 
% X2_ml_wm = linspace(min(T2_ml(mask_wm)), max(T2_ml(mask_wm)), 1000);...
% Y2_ml_wm_crb = normpdf(X2_ml_wm, t2_tru_wm(1), sigwt2_wm);
% Y2_ml_wm_tru = normpdf(X2_ml_wm, t2_tru_wm(1), sigt2_wm);
% [N2_ml_wm, B2_ml_wm] = hist(T2_ml(mask_wm), floor(sum(mask_wm(:))/500));
% figure; bar(B2_ml_wm, N2_ml_wm/max(N2_ml_wm)*max(Y2_ml_wm_tru)); hold on;...
%     plot(X2_ml_wm, Y2_ml_wm_crb, 'r', 'LineWidth', 2);...
%     plot(X2_ml_wm, Y2_ml_wm_tru, 'g', 'LineWidth', 2); hold off;...
%     legend('Empirical', 'CRLB', 'Latent');
%     if (~pr), title('T2 ML Estimates in WM');
%     else print('-depsc', strcat('T2_hist_ml_wm_3parm_', suffix, '.eps')), end;