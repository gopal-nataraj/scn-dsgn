%% Test Script for Regularized M0/T1 Recon from Synthetic BrainWeb Data
% Written by: Gopal Nataraj and Jeffrey A. Fessler, Copyright 2013

%% Synthesis of Phantom SPGR Data
% Load digital phantom; will need to modify for your directory...
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;

% Extract true M0, T1, T2* values; these are used only to simulate SPGR data
slice = 95; [nx, ny, nz] = size(labels);
M0_true = double(squeeze(labels(:,:,slice)));
T1_true = double(squeeze(labels(:,:,slice)));
T2s_true = double(squeeze(labels(:,:,slice)));
for idx = 0:10
    f = mri_brainweb_params(idx);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1;
    T2s_true(T2s_true == idx) = f.t2s;
end

% Imaging Parameters (feel free to vary flip angles, TR, TE!)
flip = [5 30]' * pi/180;                % radians 
nf = length(flip);
TR = 22.1;                              % ms
TE = 5.1;                               % ms
E1_true = exp(-TR./T1_true);

% Forward model: make true data
y_true = NaN(nx, ny, nf);
for a = 1:nf
    y_true(:,:,a) = fft2(M0_true .* exp(-TE./T2s_true) .* (1-E1_true)...
        .* sin(flip(a)) ./ (1 - E1_true .* cos(flip(a))));
end

% Add complex white gaussian noise (feel free to vary SNR!)
SNR = 40;                               % dB
sigma = exp(-SNR/20) * norm(y_true(:)) / sqrt(2*numel(y_true));
y = y_true + sigma * (randn(size(y_true)) + 1i * randn(size(y_true)));
printm('snr = %g dB', 20*log(norm(y_true(:)) / norm(col(y_true-y))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(sum-of-squares) is just the magnitude
y_im = ifft2(y);

%% (Voxel-by-voxel) Method-of-Moments Estimation
% Used as initialization for regularized algorithm; also for comparison
% Note that loops are slow in MatLab, but in practice are massively parallelizable.
E1_map_mom = NaN(nx, ny);
intercept_map_mom = NaN(nx, ny);
flipmap = permute(repmat(flip, [1 nx ny]), [2 3 1]);

for xx = 1:nx
    for yy = 1:ny
        Y = squeeze(y_im(xx,yy,:) ./ sin(flipmap(xx,yy,:)));
        X = squeeze(y_im(xx,yy,:) ./ tan(flipmap(xx,yy,:)));
        W = diag(squeeze((1-cos(flipmap(xx,yy,:)))./sin(flipmap(xx,yy,:)))); 
        A = [X ones(length(flip), 1)];
        regression = (A' * W * A) \ (A' * W * Y);
        E1_map_mom(xx,yy) = regression(1);
        intercept_map_mom(xx,yy) = regression(2);
    end
end

M0_map_mom = intercept_map_mom ./ (1-E1_map_mom);
T1_map_mom = -TR ./ log(abs(E1_map_mom));     % Should be pure real

% Project Images to within some loose, predefined range 
T1_max = 5000; 
E1_max = exp(-TR ./ T1_max);
M0_max = 2;

T1_map_mom = max(T1_map_mom, 0);
T1_map_mom = min(T1_map_mom, T1_max);
M0_map_mom = max(M0_map_mom, 0);
M0_map_mom = min(M0_map_mom, M0_max);

%% Preprocessing and Masking
% Masking
tight_mask = (M0_map_mom > 0.1) & (T1_map_mom < T1_max) & ...
    (T1_map_mom > 0) & ~isnan(T1_map_mom); 
loose_mask = imdilate(tight_mask, strel('disk', 5));

% Set voxels inside loose mask but outside tight mask to mean
M0_map_mom(~tight_mask & loose_mask) = mean(col(M0_map_mom));
T1_map_mom(~tight_mask & loose_mask) = mean(col(T1_map_mom));

% Set voxels outside both masks to zero
M0_map_mom(~loose_mask) = 0;
T1_map_mom(~loose_mask) = 0;

% Update the E1 map
E1_map_mom = exp(-TR ./ T1_map_mom);

%% Iterative Reconstruction
% Define iteration parameters (speed vs. completeness of convergence)
n_outer_iter = 30;
niter1 = 100;
niter2 = 100; 
tol_1 = 10^-6;
tol_2 = 10^-6;
disp = 1;

% Define regularizer objects R1 and R2 for M0 and T1, respectively.
% beta values control how much emphasis is placed on smoothness vs. data fidelity
% delta values control the "roundedness" of hyperbolic penalty function
beta1 = 2^-2;  
beta2 = 2^9;   
delta1 = 10^-6;
delta2 = 10^-5;
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1, 'type_penal', 'mat');
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Reshape inputs to work with 3D implementation
y_3D = permute(y_im, [1 2 4 3]);
flipmap_3D = permute(flipmap, [1 2 4 3]);

% Regularized reconstruction
[M0_map, T1_map, cost] = mri_despot1_m0t1map(M0_map_mom, E1_map_mom,...
    y_3D, flipmap_3D, loose_mask, R1, R2, T1_max, TR, n_outer_iter, ...
    niter1, niter2, tol_1, tol_2, disp);

%% Postprocessing for display and analysis

% T1 values in voxels with M0 = 0 should also be 0
thresh = 0.1;
T1_map_mom(M0_map_mom < thresh * max(col(M0_map_mom))) = 0;
T1_map(M0_map < thresh * max(col(M0_map))) = 0;

% Only examine relevant voxels; mask out rest
post_mask = labels(:,:,slice) > 0; 
M0_map_mom(~post_mask) = 0; M0_map(~post_mask) = 0; 
T1_map_mom(~post_mask) = 0; T1_map(~post_mask) = 0;

% figure; im(['notick'], fliplr(M0_map_mom), [0.4 1], ['cbar'], ' ');
% print -deps m0_mom_phant.eps;
% figure; im(['notick'], fliplr(M0_map), [0.4 1], ['cbar'], ' ');
% print -deps m0_reg_phant.eps;
% figure; im(['notick'], fliplr(M0_true), [0.4 1], ['cbar'], ' ');
% print -deps m0_tru_phant.eps;
% figure; im(['notick'], fliplr(T1_map_mom), [0 3000], ['cbar'], ' ');
% print -deps t1_mom_phant.eps;
% figure; im(['notick'], fliplr(T1_map), [0 3000], ['cbar'], ' ');
% print -deps t1_reg_phant.eps;
% figure; im(['notick'], fliplr(T1_true), [0 3000], ['cbar'], ' ');
% print -deps t1_tru_phant.eps

% Images 
figure; im(abs(M0_map_mom), [0.4 1], 'cbar');
figure; im(abs(M0_map), [0.4 1], 'cbar');
figure; im(abs(M0_true), [0.4 1], 'cbar');
figure; im(T1_map_mom, [0 4000], 'cbar');
figure; im(T1_map, [0 4000], 'cbar');
figure; im(T1_true, [0 4000], 'cbar');

% Difference Images
figure; im(abs(M0_map_mom - M0_true), [0 1], 'cbar');
figure; im(abs(M0_map - M0_true), [0 1], 'cbar');
figure; im(abs(T1_map_mom - T1_true), [0 4000], 'cbar');
figure; im(abs(T1_map - T1_true), [0 4000], 'cbar');

% Plotting cost vs. iteration
figure; hold on;
scatter([1:2:2*n_outer_iter], cost(2:2:end), 'bo');
scatter([2:2:2*n_outer_iter], cost(3:2:end), 'ro');
plot([0:2*n_outer_iter], cost, 'k'); hold off;
title('Cost vs. iteration');
legend('After x1 (M0) update','After x2 (E1) update');

% %% Analysis
% 
% % UNREALISTICALLY compensate for T2s 
% comp_mask = T2s_true > 0;
% M0_map_comp = zeros(size(M0_map));
% M0_map_comp(comp_mask) = M0_map(comp_mask) ./ exp(-TE ./ T2s_true(comp_mask));
% figure; im(abs(M0_map_comp), [0.4 1], 'cbar');
% figure; im(abs(M0_map_comp - M0_true), [0 1], 'cbar');
% 
% % NRMSE Comparison
% nrmse = @(tst, tru) sqrt(sum((abs(tst(:))-tru(:)).^2) ./ sum(abs(tru(:)).^2));
% err_m0_mom = nrmse(M0_map_mom, M0_true);
% err_m0_reg = nrmse(M0_map, M0_true);
% err_m0_corr = nrmse(M0_map_comp, M0_true);
% 
% err_t1_mom = nrmse(T1_map_mom, T1_true);
% err_t1_reg = nrmse(T1_map, T1_true);
