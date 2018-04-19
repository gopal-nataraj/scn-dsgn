% script digiphant_2016_02_17_se_model_test.m
% simulations to find whether model mismatch due to finite TR and T1
% effects can significantly influence T2 estimates
%
% Written by: Gopal Nataraj
% Last modified: 2016-02-17

%% load true parameter maps and set scan protocol
% load digital phantom
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../matlab;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;
addpath('Calculations/SE');
addpath('Helper/');

% get true values
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

% SE imaging parameters
TE = [10 150]';                 % ms
nTE = length(TE);
TR = 3000 * ones(nTE, 1);       % ms
flip_ex = pi/2;                 % rad   (nominal)
flip_ref = pi;                  % rad   (nominal)
dw = 0;                         % kHz   (no off-resonance)
is_mag = 1;

% create a true excitation profile map
% cannot allow too much variation in flip because then ML estimate biased
mask = M0_true == 0;
perc_var = 0.30;                % allow +/-0% variation in flip
[xx, yy] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
kap_true = (1+perc_var/2) - perc_var*(xx.^2 + yy.^2); 
kap_true(mask) = 0;             % No excitation where M0 = 0;
figure, im('notick', fliplr(kap_true), [0.9 1.1], 'cbar', ' ');

%% dictionary creation (could be precomputed)
% to test model mismatch, make dictionary with 'old' model
T2 = logspace(1, 3, 10000);
D = NaN(nTE, length(T2));
K = length(T2);

for t2 = 1:length(T2)
    for i = 1:nTE
        D(i, t2) = SE_fun(1, T2(t2), 1, flip_ex, TE(i), is_mag);
    end
end

%% sythesize SE data
% to test model mismatch, synthesize with 'new' model
y_true = NaN(nx, ny, nTE);
for i = 1:nTE
%     y_true(:,:,i) = fft2(SE_fun_v2(M0_true, T1_true, T2_true,...
%         TE(i), TR(i), is_mag));
%     y_true(:,:,i) = fft2(SE_fun_v3(M0_true, T1_true, T2_true, kap_true,...
%         TE(i), TR(i), flip_ex, flip_ref, dw, is_mag));
    y_im_true = SE_fun_v3(M0_true, T1_true, T2_true, kap_true,...
        TE(i), TR(i), flip_ex, flip_ref, dw, is_mag);
    y_true(:,:,i) = fft2(y_im_true);
end

% add complex white gaussian noise
var_im = 0;
% var_im = 1.31e-7; 
sigma_k = sqrt(var_im * (nx*ny));
y = y_true + sigma_k * (randn(size(y_true)) + 1i*randn(size(y_true)));
printm('snr = %g dB', 20*log10(norm(y_true(:)) / norm(col(y_true-y))));

% take noisy data back to image domain
% since simulating only body coil, sqrt(SoS) is just the magnitude
y_im = abs(ifft2(y));

%% maximum-likelihood estimation
% dictionary-based estimation via variable projection method
weights = ones(nTE, 1);
W = spdiags(weights, 0, nTE, nTE);
y = reshape(permute(y_im, [3 1 2]), [nTE nx*ny]);

tic;
maxProd = zeros(1, nx*ny);
t2_idx = ones(1, nx*ny);
for k = 1:K
    % compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;
    
    % if the kth inner produt is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    t2_idx(update) = k;
end
time_ml = toc;

% extract indices for maximum-likelihood maps
T2_ml = reshape(T2(t2_idx), [nx ny]);

%% images, rmse, and comparisons
% noise scaling of display
eps_scl = 0.1;

% postprocessing for display
T2_ml(mask) = 0;

% t2 images
figure; im('notick', fliplr(T2_ml), [0 120], 'cbar', ' ');
figure; im('notick', fliplr(T2_true), [0 120], 'cbar', ' ');

figure; im('notick', fliplr(abs(T2_ml-T2_true)), [0 eps_scl*120], 'cbar', ' ');

% compute nrmse and rmse
nrmse = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ sum(tru(:).^2));
rmse  = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ numel(tru));

% gm and wm masks for task-based performance assessment
mask_gm = labels(:,:,slice) == 2;
mask_wm = labels(:,:,slice) == 3;

% t2 bias values
bias_t2_ml_gm       = mean(T2_ml(mask_gm) - T2_true(mask_gm))
bias_t2_tru_gm      = mean(T2_true(mask_gm) - T2_true(mask_gm));

bias_t2_ml_wm       = mean(T2_ml(mask_wm) - T2_true(mask_wm))
bias_t2_tru_wm      = mean(T2_true(mask_wm) - T2_true(mask_wm));

% t2 rmse values
rmse_t2_ml_gm       = rmse(T2_ml(mask_gm),  T2_true(mask_gm));
rmse_t2_tru_gm      = rmse(T2_true(mask_gm), T2_true(mask_gm));

rmse_t2_ml_wm       = rmse(T2_ml(mask_wm),  T2_true(mask_wm));
rmse_t2_tru_wm      = rmse(T2_true(mask_wm), T2_true(mask_wm));

% t2 means and standard deviations
t2_ml_gm            = [mean(T2_ml(mask_gm)),   std(T2_ml(mask_gm))];
t2_tru_gm           = [mean(T2_true(mask_gm)), std(T2_true(mask_gm))];

t2_ml_wm            = [mean(T2_ml(mask_wm)),   std(T2_ml(mask_wm))];
t2_tru_wm           = [mean(T2_true(mask_wm)), std(T2_true(mask_wm))];