%% Header Script for Regularized M0/T1/T2 Recon from Synthetic Data
% This script resolves the model mismatch: the SPGR output is effectively
% M0star, and the expected input into DESS is also M0star.
% 
% Written By: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

%% Synthesis of Phantom DESS Data
% Load digital phantom
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;

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

% % Make a T2p_true map, where 1/T2s = 1/T2p + 1/T2
% T2p_msk = T2_true ~= T2s_true; 
% T2p_true = zeros(size(T2_true));
% T2p_true(T2p_msk) = (T2_true(T2p_msk) .* T2s_true(T2p_msk)) ./ ...
%     (T2_true(T2p_msk) - T2s_true(T2p_msk));

%% DESPOT1 Scan Simulation
% Imaging Parameters
flip = [5 22]' * pi/180;                % radians 
% flip = [5 30]' * pi/180;                % radians 
nf = length(flip);
TR = 12;                                % ms
% TR = 20;                                % ms
TE = 5;                                 % ms 
E1_true = exp(-TR./T1_true);
E2_true = exp(-TR./T2_true);

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = T2_true ~= T2s_true; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TE ./ T2s_true(T2s_msk));

% Forward model: make true data
y_true = NaN(nx, ny, nf);
for a = 1:nf
%     y_true(:,:,a) = fft2(M0_true .* exp(-TE./T2s_true) .* (1-E1_true)...
%         .* sin(flip(a)) ./ (1 - E1_true .* cos(flip(a))));
    y_true(:,:,a) = fft2(M0s_true .* (1-E1_true)...
        .* sin(flip(a)) ./ (1 - E1_true .* cos(flip(a))));
end

% Add complex white gaussian noise
% SNR = 40;                               % dB
% sigma = exp(-SNR/20) * norm(y_true(:)) / sqrt(2*numel(y_true));
sigma = 1;

y = y_true + sigma * (randn(size(y_true)) + 1i * randn(size(y_true)));
printm('snr = %g dB', 20*log(norm(y_true(:)) / norm(col(y_true-y))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
y_im = ifft2(y);

% Method-of-Moments M0,T1 Estimation
E1_mom = NaN(nx, ny);
intercept_map_mom = NaN(nx, ny);
flipmap = permute(repmat(flip, [1 nx ny]), [2 3 1]);

for xx = 1:nx
    for yy = 1:ny
        Y = squeeze(y_im(xx,yy,:) ./ sin(flipmap(xx,yy,:)));
        X = squeeze(y_im(xx,yy,:) ./ tan(flipmap(xx,yy,:)));
        W = diag(squeeze((1-cos(flipmap(xx,yy,:)))./sin(flipmap(xx,yy,:)))); 
        A = [X ones(length(flip), 1)];
        regression = (A' * W * A) \ (A' * W * Y);
        E1_mom(xx,yy) = regression(1);
        intercept_map_mom(xx,yy) = regression(2);
    end
end

M0s_mom = intercept_map_mom ./ (1-E1_mom);
T1_mom = -TR ./ log(abs(E1_mom));     % Should be pure real

% Project Images to within range 
T1_max = 5000; 
E1_max = exp(-TR ./ T1_max);
M0s_max = 2;

T1_mom = max(T1_mom, 0);
T1_mom = min(T1_mom, T1_max);
M0s_mom = max(M0s_mom, 0);
M0s_mom = min(M0s_mom, M0s_max);

% Preprocessing and Masking
tight_mask = (M0s_mom > 0.1) & (T1_mom < T1_max) & ...
    (T1_mom > 0) & ~isnan(T1_mom); 
loose_mask = imdilate(tight_mask, strel('disk', 5));

% Set voxels inside loose mask but outside tight mask to mean
M0s_mom(~tight_mask & loose_mask) = mean(col(M0s_mom));
T1_mom(~tight_mask & loose_mask) = mean(col(T1_mom));

% Set voxels outside both masks to zero
M0s_mom(~loose_mask) = 0;
T1_mom(~loose_mask) = 0;

% Median filtering (M0s_med might not be right...)
M0s_med = medfilt2(real(M0s_mom)) + 1i*medfilt2(imag(M0s_mom)); 
T1_med = medfilt2(T1_mom);

% Update the E1 maps
E1_mom = exp(-TR ./ T1_mom);
E1_med = exp(-TR ./ T1_med);

% Define iteration parameters
n_outer_iter = 30;
niter1 = 200;
niter2 = 200; 
tol_1 = 10^-5;
tol_2 = 10^-5;
disp = 0;

% Define regularizers, R1 and R2
beta1 = 2^-1;   %2^-2;
beta2 = 2^9;    %2^7;
delta1 = 10^-6; 
delta2 = 10^-5;
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1, 'type_penal', 'mat');
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Regularized reconstruction
% Note that M0s_reg is actually M0* := M0exp(-TE/T2s)
y_tmp = permute(y_im, [1 2 4 3]);
flipmap_tmp = permute(flipmap, [1 2 4 3]);
mask_tmp = permute(permute(loose_mask, [3 2 1]), [3 2 1]);

% [M0s_reg, T1_reg, cost] = mri_despot1_m0t1map(M0s_mom, E1_med,...
%     y_tmp, flipmap_tmp, loose_mask, R1, R2, T1_max, TR, n_outer_iter, ...
%     niter1, niter2, tol_1, tol_2, disp);
tic;
[M0s_reg, T1_reg, cost] = mri_despot1_m0t1map(M0s_med, E1_med,...
    y_tmp, flipmap_tmp, loose_mask, R1, R2, T1_max, TR, n_outer_iter, ...
    niter1, niter2, tol_1, tol_2, disp);
time_despot1 = toc;

% Postprocessing for display and analysis
% T1 values in voxels with M0 = 0 should also be 0
thresh = 0.1;
T1_reg(M0s_reg < thresh * max(col(M0s_reg))) = 0;

% % Median-filter the outputs
% M0s_reg = medfilt2(real(M0s_reg)) + 1i*medfilt2(imag(M0s_reg)); 
% T1_reg = medfilt2(T1_reg);

% M0s and T1 values outside tight_mask should also be 0
M0s_reg(~tight_mask) = 0;
T1_reg(~tight_mask) = 0;

% Only examine relevant voxels; mask out rest
post_mask = labels(:,:,slice) > 0; 
M0s_mom(~post_mask) = 0; M0s_med(~post_mask) = 0; M0s_reg(~post_mask) = 0; 
T1_mom(~post_mask) = 0; T1_med(~post_mask) = 0; T1_reg(~post_mask) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DESS Scan Simulation
% Imaging Parameters 
flip = [47]' * pi/180;                  % radians
nf = length(flip);
TR = 16;                                % ms 
% TR = 20;                                % ms 
TE = 5;                                 % ms
wf_true = 0;                            % Off-resonance
SNR = 60;                               % dB

% Forward model: make true data, using M0s_true
yp_true = NaN(nx, ny, nf);
ym_true = NaN(nx, ny, nf);
for a = 1:nf
    v1 = (1 - E1_true * cos(flip(a))) ./ (E1_true - cos(flip(a)));
    [S1true, S2true] = dess_fun_M0star(TR, TE, M0s_true, wf_true, flip(a), v1, T2_true);
    yp_true(:,:,a) = fft2(S1true);
    ym_true(:,:,a) = fft2(S2true);
end

% Add complex white gaussian noise
% sigma_p = exp(-SNR/20) * norm(yp_true(:)) / sqrt(2*numel(yp_true));
% sigma_m = exp(-SNR/20) * norm(ym_true(:)) / sqrt(2*numel(ym_true));
sigma_p = 1;
sigma_m = 1;

yp = yp_true + sigma_p * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
ym = ym_true + sigma_m * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
printm('snr_p = %g dB', 20*log(norm(yp_true(:)) / norm(col(yp_true-yp))));
printm('snr_m = %g dB', 20*log(norm(ym_true(:)) / norm(col(ym_true-ym))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
yp_im = ifft2(yp);
ym_im = ifft2(ym);

% Method-of-Moments T2 Estimation
if (nf == 1)
    T2_mom = -2*(TR-TE) ./ log(abs(ym_im ./ yp_im));
else
    E2squared = NaN(nx, ny);
    for xx = 1:nx
        for yy = 1:ny 
            Y = squeeze(ym_im(xx,yy,:));
            A = squeeze(yp_im(xx,yy,:));
            C = [A Y];
            sig_small = min(svds(C'*C));
            E2squared(xx,yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
        end
    end
    T2_mom = -2*(TR-TE) ./ log(E2squared);
end

% For now, take true value wf=0
wf_mom = zeros(nx, ny);

% Project Images to within range
T2max = 500; E2max = exp(-TR ./ T2max);
T2min = 5;  E2min = exp(-TR ./ T2min);
T2_mom = max(T2_mom, T2min);
T2_mom = min(T2_mom, T2max); 

% Median filtering of initial estimate
T2_med = medfilt2(T2_mom);

% Preprocessing and Masking
tight_mask = M0s_reg > 0.1;
loose_mask = imdilate(tight_mask, strel('disk', 5));

% Set voxels inside loose mask but outside tight mask to mean
% M0s_reg(~tight_mask & loose_mask) = mean(col(M0s_reg));
% T1_reg(~tight_mask & loose_mask) = mean(col(T1_reg));
T2_mom(~tight_mask & loose_mask) = mean(col(T2_mom));
wf_mom(~tight_mask & loose_mask) = mean(col(wf_mom));

% Set voxels outside both tight and loose masks to zero
% M0s_reg(~loose_mask) = 0;
% T1_reg(~loose_mask) = 0;
T2_mom(~loose_mask) = 0;
wf_mom(~loose_mask) = 0;

% Update the E1 and E2 maps
E1_reg = exp(-TR ./ T1_reg);
E2_mom = exp(-TR ./ T2_mom);
E2_med = exp(-TR ./ T2_med);

% Define iteration parameters
n_outer = 30;
n_inner2 = 200;
tol2 = 10^-5;
disp = 0;

% Define regularizers
beta2 = 2^7;
delta2 = 10^-5;
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Reshape inputs to work with 3-D implementation
flip_im = permute(repmat(flip, [1 nx ny]), [2 3 1]);
flip_3D = permute(flip_im, [1 2 4 3]);
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);
%mask_tmp = permute(permute(loose_mask, [3 2 1]), [3 2 1]);

% Regularized reconstruction (Nesterov momentum only)
% Now the input and output are BOTH expected to be M0s := M0exp(-TE/T2s)
% Note that for this to work, TE must remain the same for all scans.
tic;
[M0s_reg_n83, T1_reg_n83, T2_reg_n83, wf_reg_n83, cost_n83] = ...
    mri_dess_t2map_nest83_M0s(M0s_reg, E1_reg, E2_med, wf_mom, ...
    flip_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
    TR, TE, R2, n_outer, n_inner2, tol2, disp);
time_dess = toc;

% % Sanity check: perturbing the solution shouldn't change the answer
% SNR_e2 = 50;
% sigma_e2 = exp(-SNR_e2/20) * norm(E2_true(:)) / sqrt(numel(E2_true));
% E2_med = E2_true + sigma_e2 * randn(size(E2_true));
% printm('snr = %g dB', 20*log(norm(E2_true(:)) / norm(col(E2_true-E2_med))));

% % Sanity check: regularized reconstruction with true M0, T1
% [M0s_reg_n83, T1_reg_n83, T2_reg_n83, wf_reg_n83, cost_n83] = ...
%     mri_dess_t2map_nest83_M0s(M0s_true, E1_true, E2_med, wf_mom, ...
%     flip_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
%     TR, TE, R2, n_outer, n_inner2, tol2, disp);

% Postprocessing for display and analysis
% T2 values in voxels with M0 = 0 should also be 0
T2_reg_n83(M0s_reg_n83 <= 0.1*max(M0s_reg_n83(:))) = 0; 

% Remove pixels outside tight_mask
T2_mom(~tight_mask) = 0;
T2_med(~tight_mask) = 0;
T2_reg_n93(~tight_mask) = 0;

%% Images, RMSE, and Comparisons

% Noise scaling of display
eps_scl = 0.2;

% M0s Images 
figure; im('notick', fliplr(abs(M0s_mom)), [0 1], 'cbar', ' ');
% print -deps M0s_mom_seq_short.eps;
figure; im('notick', fliplr(abs(M0s_med)), [0 1], 'cbar', ' ');
% print -deps M0s_med_seq_short.eps;
figure; im('notick', fliplr(medfilt2(abs(M0s_reg))), [0 1], 'cbar', ' ');
% print -deps M0s_reg_seq_short.eps;
figure; im('notick', fliplr(abs(M0s_true)), [0 1], 'cbar', ' ');
% print -deps M0s_tru_seq_short.eps;

figure; im('notick', fliplr(abs(M0s_mom-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_mom_err_seq_short.eps;
figure; im('notick', fliplr(abs(M0s_med-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_med_err_seq_short.eps;
figure; im('notick', fliplr(abs(M0s_reg-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_reg_err_seq_short.eps;

% T1 Images
T1_med(~loose_mask) = 0;
figure; im('notick', fliplr(T1_mom), [0 2000], 'cbar', ' ');
% print -deps T1_mom_seq_short.eps;
figure; im('notick', fliplr(T1_med), [0 2000], 'cbar', ' ');
% print -deps T1_med_seq_short.eps;
figure; im('notick', fliplr(T1_reg), [0 2000], 'cbar', ' ');
% print -deps T1_reg_seq_short.eps;
figure; im('notick', fliplr(T1_true), [0 2000], 'cbar', ' ');
% print -deps T1_tru_seq_short.eps;

figure; im('notick', fliplr(T1_mom-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_mom_err_seq_short.eps;
figure; im('notick', fliplr(T1_med-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_med_err_seq_short.eps;
figure; im('notick', fliplr(T1_reg-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_reg_err_seq_short.eps;

% T2 Images
figure; im('notick', fliplr(T2_mom), [0 200], 'cbar', ' ');
% print -deps T2_mom_seq_short.eps;
figure; im('notick', fliplr(T2_med), [0 200], 'cbar', ' ');
% print -deps T2_med_seq_short.eps;
figure; im('notick', fliplr(T2_reg_n83), [0 200], 'cbar', ' ');
% print -deps T2_reg_seq_short.eps;
figure; im('notick', fliplr(T2_true), [0 200], 'cbar', ' ');
% print -deps T2_tru_seq_short.eps;

figure; im('notick', fliplr(abs(T2_mom-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_mom_err_seq_short.eps;
figure; im('notick', fliplr(abs(T2_med-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_med_err_seq_short.eps;
figure; im('notick', fliplr(abs(T2_reg_n83-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_reg_err_seq_short.eps;

% Compute NRMSE and RMSE
nrmse = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ sum(tru(:).^2));
rmse  = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ numel(tru));

% GM and WM masks for task-based performance assessment
mask_gm = labels(:,:,slice) == 2;
mask_wm = labels(:,:,slice) == 3;

% M0s NRMSE values
nrmse_m0s_mom_gm = nrmse(M0s_mom(mask_gm), M0s_true(mask_gm));
nrmse_m0s_med_gm = nrmse(M0s_med(mask_gm), M0s_true(mask_gm));
nrmse_m0s_reg_gm = nrmse(M0s_reg(mask_gm), M0s_true(mask_gm));

nrmse_m0s_mom_wm = nrmse(M0s_mom(mask_wm), M0s_true(mask_wm));
nrmse_m0s_med_wm = nrmse(M0s_med(mask_wm), M0s_true(mask_wm));
nrmse_m0s_reg_wm = nrmse(M0s_reg(mask_wm), M0s_true(mask_wm));

% T1 RMSE values
rmse_t1_mom_gm = rmse(T1_mom(mask_gm), T1_true(mask_gm));
rmse_t1_med_gm = rmse(T1_med(mask_gm), T1_true(mask_gm));
rmse_t1_reg_gm = rmse(T1_reg(mask_gm), T1_true(mask_gm));

rmse_t1_mom_wm = rmse(T1_mom(mask_wm), T1_true(mask_wm));
rmse_t1_med_wm = rmse(T1_med(mask_wm), T1_true(mask_wm));
rmse_t1_reg_wm = rmse(T1_reg(mask_wm), T1_true(mask_wm));

% T2 RMSE values
rmse_t2_mom_gm = rmse(T2_mom(mask_gm), T2_true(mask_gm));
rmse_t2_med_gm = rmse(T2_med(mask_gm), T2_true(mask_gm));
rmse_t2_reg_gm = rmse(T2_reg_n83(mask_gm), T2_true(mask_gm));

rmse_t2_mom_wm = rmse(T2_mom(mask_wm), T2_true(mask_wm));
rmse_t2_med_wm = rmse(T2_med(mask_wm), T2_true(mask_wm));
rmse_t2_reg_wm = rmse(T2_reg_n83(mask_wm), T2_true(mask_wm));

% M0s Means and Standard Deviations
m0s_mom_gm = [mean(M0s_mom(mask_gm)), std(M0s_mom(mask_gm))];
m0s_med_gm = [mean(M0s_med(mask_gm)), std(M0s_med(mask_gm))];
m0s_reg_gm = [mean(M0s_reg(mask_gm)), std(M0s_reg(mask_gm))];
m0s_tru_gm = [mean(M0s_true(mask_gm)), std(M0s_true(mask_gm))];

m0s_mom_wm = [mean(M0s_mom(mask_wm)), std(M0s_mom(mask_wm))];
m0s_med_wm = [mean(M0s_med(mask_wm)), std(M0s_med(mask_wm))];
m0s_reg_wm = [mean(M0s_reg(mask_wm)), std(M0s_reg(mask_wm))];
m0s_tru_wm = [mean(M0s_true(mask_wm)), std(M0s_true(mask_wm))];

% T1 Means and Standard Deviations
t1_mom_gm = [mean(T1_mom(mask_gm)), std(T1_mom(mask_gm))];
t1_med_gm = [mean(T1_med(mask_gm)), std(T1_med(mask_gm))];
t1_reg_gm = [mean(T1_reg(mask_gm)), std(T1_reg(mask_gm))];
t1_tru_gm = [mean(T1_true(mask_gm)), std(T1_true(mask_gm))];

t1_mom_wm = [mean(T1_mom(mask_wm)), std(T1_mom(mask_wm))];
t1_med_wm = [mean(T1_med(mask_wm)), std(T1_med(mask_wm))];
t1_reg_wm = [mean(T1_reg(mask_wm)), std(T1_reg(mask_wm))];
t1_tru_wm = [mean(T1_true(mask_wm)), std(T1_true(mask_wm))];

% T2 Means and Standard Deviations
t2_mom_gm = [mean(T2_mom(mask_gm)), std(T2_mom(mask_gm))];
t2_med_gm = [mean(T2_med(mask_gm)), std(T2_med(mask_gm))];
t2_reg_gm = [mean(T2_reg_n83(mask_gm)), std(T2_reg_n83(mask_gm))];
t2_tru_gm = [mean(T2_true(mask_gm)), std(T2_true(mask_gm))];

t2_mom_wm = [mean(T2_mom(mask_wm)), std(T2_mom(mask_wm))];
t2_med_wm = [mean(T2_med(mask_wm)), std(T2_med(mask_wm))];
t2_reg_wm = [mean(T2_reg_n83(mask_wm)), std(T2_reg_n83(mask_wm))];
t2_tru_wm = [mean(T2_true(mask_wm)), std(T2_true(mask_wm))];

