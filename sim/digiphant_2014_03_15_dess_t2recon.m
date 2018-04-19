%% Header Script for Regularized T2 Recon from Synthetic Data

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

% Make a T2p_true map, where 1/T2s = 1/T2p + 1/T2
T2p_msk = T2_true ~= T2s_true; 
T2p_true = zeros(size(T2_true));
T2p_true(T2p_msk) = (T2_true(T2p_msk) .* T2s_true(T2p_msk)) ./ ...
    (T2_true(T2p_msk) - T2s_true(T2p_msk));

% Imaging Parameters
%flip = [15 30 45 60 75 90]' * pi/180;   % radians 
flip = [45]' * pi/180;                  % radians
bst_flp_idx = 1;                        % choose a flip for initial guess
nf = length(flip);
TR = 20;                                % ms
TE = 5;                                 % ms
wf_true = 0;                            % Off-resonance
E1_true = exp(-TR./T1_true);
E2_true = exp(-TR./T2_true);

% Make a M0p_true map, where M0p = M0 * exp(-TE/T2p)
% This is to be used for generating the forward model 
M0p_true = M0_true;
M0p_true(T2p_msk) = M0_true(T2p_msk) .* exp(-TE ./ T2p_true(T2p_msk));

% Forward model: make true data
% SHOULD be using M0p_true; instead using M0_true --> better RMSE...
yp_true = NaN(nx, ny, nf);
ym_true = NaN(nx, ny, nf);
for a = 1:nf
    v1 = (1 - E1_true * cos(flip(a))) ./ (E1_true - cos(flip(a)));
    [S1true, S2true] = dess_fun(TR, TE, M0_true, wf_true, flip(a), v1, T2_true);
    yp_true(:,:,a) = fft2(S1true);
    ym_true(:,:,a) = fft2(S2true);
end

% Add complex white gaussian noise
SNR = 40;                               % dB
sigma_p = exp(-SNR/20) * norm(yp_true(:)) / sqrt(2*numel(yp_true));
sigma_m = exp(-SNR/20) * norm(ym_true(:)) / sqrt(2*numel(ym_true));
yp = yp_true + sigma_p * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
ym = ym_true + sigma_m * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
printm('snr_p = %g dB', 20*log(norm(yp_true(:)) / norm(col(yp_true-yp))));
printm('snr_m = %g dB', 20*log(norm(ym_true(:)) / norm(col(ym_true-ym))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
yp_im = ifft2(yp);
ym_im = ifft2(ym);

%% Method-of-Moments T2 Estimation
T2_mom = -2*(TR-TE) ./ log(abs(ym_im(:,:,bst_flp_idx) ./ yp_im(:,:,bst_flp_idx)));

% %% Method-of-Moments T1 Estimation
% T1_mom = NaN(nx, ny, nf);
% M0_mom = ones(nx, ny);                  % Rough guess
% for a = 1:nf
%     E2E_mom = exp(-TE ./ T2_mom(:,:,a));
%     beta = (1- (squeeze(ym_im(:,:,a)).*E2E_mom./M0_mom./tan(flip(a)/2))) ./ ...
%         (1 - (squeeze(yp_im(:,:,a))./E2E_mom./M0_mom./tan(flip(a)/2)));
%     T1_mom(:,:,a) = -TR ./ log((1+beta*cos(flip(a))) ./ (beta+cos(flip(a))));
% end
% figure; im(T1_mom, [0 3000], 'cbar');
% figure; im(T1_true, [0 3000], 'cbar');

% For now, take true values of M0, T1 and wf=0
M0_mom = M0_true;
T1_mom = T1_true;
wf_mom = zeros(nx, ny);

% Project Images to within range
T2max = 1000; E2max = exp(-TR ./ T2max);
T2min = 5;  E2min = exp(-TR ./ T2min);
T2_mom = max(T2_mom, T2min);
T2_mom = min(T2_mom, T2max); 

%% Preprocessing

% Median filtering of initial estimate
T2_med = medfilt2(T2_mom);

% Masking
tight_mask = M0_mom > 0.1;
loose_mask = imdilate(tight_mask, strel('disk', 5));

% Set voxels inside loose mask but outside tight mask to mean
M0_mom(~tight_mask & loose_mask) = mean(col(M0_mom));
T1_mom(~tight_mask & loose_mask) = mean(col(T1_mom));
T2_mom(~tight_mask & loose_mask) = mean(col(T2_mom));
wf_mom(~tight_mask & loose_mask) = mean(col(wf_mom));

% Set voxels outside both tight and loose masks to zero
M0_mom(~loose_mask) = 0;
T1_mom(~loose_mask) = 0;
T2_mom(~loose_mask) = 0;
wf_mom(~loose_mask) = 0;

% Update the E1 and E2 maps
E1_mom = exp(-TR ./ T1_mom);
E2_mom = exp(-TR ./ T2_mom);
E2_med = exp(-TR ./ T2_med);

%% Iterative Reconstruction
% Define iteration parameters
n_outer = 20;
n_inner2 = 200;
tol2 = 10^-6;
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

% Regularized reconstruction comparing methods
% [M0_reg_pgd, T1_reg_pgd, T2_reg_pgd, wf_reg_pgd, cost_pgd] = ...
%     mri_dess_t2map_pgd(M0_mom, E1_mom, E2_med, wf_mom, ...
%     flip_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
%     TR, TE, R2, n_outer, n_inner2, tol2, disp);
[M0_reg_n83, T1_reg_n83, T2_reg_n83, wf_reg_n83, cost_n83] = ...
    mri_dess_t2map_nest83(M0_mom, E1_mom, E2_med, wf_mom, ...
    flip_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
    TR, TE, R2, n_outer, n_inner2, tol2, disp);

% Output images
% figure; im(abs(M0_mom), [0.4 1], 'cbar');
% figure; im(abs(M0_reg), [0.4 1], 'cbar');
% figure; im(abs(M0_true), [0.4 1], 'cbar');
% figure; im(T1_mom, [0 3000], 'cbar');
% figure; im(T1_reg, [0 3000], 'cbar');
% figure; im(T1_true, [0 3000], 'cbar');

T2_med(~loose_mask) = 0;
figure; im('notick', T2_mom, [0 200], 'cbar', ' ');
% print -deps T2_mom_phant.eps;
figure; im('notick', T2_med, [0 200], 'cbar', ' ');
% print -deps T2_med_phant.eps;
figure; im('notick', T2_reg_n83, [0 200], 'cbar', ' ');
% print -deps T2_reg_phant.eps;
figure; im('notick', T2_true, [0 200], 'cbar', ' ');
% print -deps T2_tru_phant.eps;

figure; im('notick', abs(T2_mom-T2_true), [0 200], 'cbar', ' ');
% print -deps T2_mom_err_phant.eps;
figure; im('notick', abs(T2_med-T2_true), [0 200], 'cbar', ' ');
% print -deps T2_med_err_phant.eps;
figure; im('notick', abs(T2_reg_n83-T2_true), [0 200], 'cbar', ' ');
% print -deps T2_reg_err_phant.eps;
figure; im('notick', abs(T2_true-T2_true), [0 200], 'cbar', ' ');
% print -deps T2_tru_err_phant.eps;

% figure; im(wf_mom, [-pi/2 pi/2], 'cbar');
% figure; im(wf_reg, [-pi/2 pi/2], 'cbar');

% % Plotting cost vs. iteration
% figure; hold on;
% scatter([0:n_inner2:n_inner2*n_outer], cost_pgd, 'bo');
% scatter([0:n_inner2:n_inner2*n_outer], cost_n83, 'ro');
% plot([0:n_inner2:n_inner2*n_outer], cost_pgd, 'k');
% plot([0:n_inner2:n_inner2*n_outer], cost_n83, 'k');
% hold off; title('Cost vs. iteration');
% legend('PGD','Nesterov-83');
% % print -depsc cost_phant.eps;

% Compute NRMSE
nrmse_mom = sqrt(sum((abs(T2_mom(:) - T2_true(:))).^2) ./ sum(T2_true(:).^2))
nrmse_med = sqrt(sum((abs(T2_med(:) - T2_true(:))).^2) ./ sum(T2_true(:).^2))
nrmse_reg = sqrt(sum((abs(T2_reg_n83(:) - T2_true(:))).^2) ./ sum(T2_true(:).^2))

rmse_mom = sqrt(sum((abs(T2_mom(:) - T2_true(:))).^2) ./ numel(T2_true))
rmse_med = sqrt(sum((abs(T2_med(:) - T2_true(:))).^2) ./ numel(T2_true))
rmse_reg = sqrt(sum((abs(T2_reg_n83(:) - T2_true(:))).^2) ./ numel(T2_true))

% % Check to make sure that the signal model matches the raw data
% Sp = fp(M0_mom, E1_mom, E2_true, pi/4, TR, TE, wf_mom);
% Sm = fm(M0_mom, E1_mom, E2_true, pi/4, TR, TE, wf_mom);
% figure; im(abs(Sp), [0 0.2], 'cbar'); title('Sp');
% figure; im(abs(Sm), [0 0.2], 'cbar'); title('Sm'); %nan
% figure; im(abs(ifft2(yp_true(:,:,3))), [0 0.2], 'cbar'); title('yp');
% figure; im(abs(ifft2(ym_true(:,:,3))), [0 0.2], 'cbar'); title('ym');
% figure; im(abs(Sp - ifft2(yp_true(:,:,3))), [0 0.2], 'cbar');
% figure; im(abs(Sm - ifft2(ym_true(:,:,3))), [0 0.2], 'cbar');

