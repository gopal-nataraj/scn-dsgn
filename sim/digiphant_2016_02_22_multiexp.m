% script digiphant_2016_02_22_multiexp.m
% simulations to find whether model mismatch due to multiexponential T2
%
% Written by: Gopal Nataraj
% Copyright 2016, University of Michigan
%
% Version control
%   v1.1      2016-02-22      spin echo modeling only
%   v1.2      2016-07-16      minor tweaks for paper

%% load true parameter maps 
% load digital phantom
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/DigitalPhantom');
addpath('../../model/se');
addpath('../../etc');

% load phantom parameters
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
labels = fld_read(f.filename);

% make true multicompartmental t2 maps
% tissue weights (myelin-rich, myelin-starved, csf)
wghts_wm = [0.20 0.75 0.05]';
wghts_gm = [0.00 0.95 0.05]';
wghts_sf = [0.00 0.00 1.00]';
wghts_df = [0.00 0.00 0.00]';           % set everything else to zero

% compartment t2 (at 1.5T)
t2_myl = exp(mean([log(10) log(20)]));
t2_dmy = exp(mean([log(70) log(100)]));
t2_csf = 2000;
t2_cmp = [t2_myl t2_dmy t2_csf]';

% make true parameter maps
b0 = 3.0;                               % tesla                        
slice = 81; 
[nx, ny, nz] = size(labels);
M0_true     = double(squeeze(labels(:,:,slice)));
T1_true     = double(squeeze(labels(:,:,slice)));
T2s_true    = double(squeeze(labels(:,:,slice)));

T2eff_true  = double(squeeze(labels(:,:,slice)));
myl_true    = double(squeeze(labels(:,:,slice)));
dmy_true    = double(squeeze(labels(:,:,slice)));
csf_true    = double(squeeze(labels(:,:,slice)));

for idx = 0:10
  f = mri_brainweb_params(idx, 'b0', b0);
  M0_true(M0_true == idx) = f.pd;
  T1_true(T1_true == idx) = f.t1; 
  T2s_true(T2s_true == idx) = f.t2s;
    
  % t2eff uses the volumetric mean of the rates in each voxel
  switch idx
    case 1  % csf
      T2eff_true(T2eff_true == idx) = div0(1, wghts_sf' * div0(1, t2_cmp));
      myl_true(myl_true == idx) = wghts_sf(1);
      dmy_true(dmy_true == idx) = wghts_sf(2);
      csf_true(csf_true == idx) = wghts_sf(3);
            
    case 2  % gm
      T2eff_true(T2eff_true == idx) = div0(1, wghts_gm' * div0(1, t2_cmp));
      myl_true(myl_true == idx) = wghts_gm(1);
      dmy_true(dmy_true == idx) = wghts_gm(2);
      csf_true(csf_true == idx) = wghts_gm(3);

    case 3  % wm
      T2eff_true(T2eff_true == idx) = div0(1, wghts_wm' * div0(1, t2_cmp));
      myl_true(myl_true == idx) = wghts_wm(1);
      dmy_true(dmy_true == idx) = wghts_wm(2);
      csf_true(csf_true == idx) = wghts_wm(3);

    otherwise
      T2eff_true(T2eff_true == idx) = div0(1, wghts_df' * div0(1, t2_cmp));
      myl_true(myl_true == idx) = wghts_df(1);
      dmy_true(dmy_true == idx) = wghts_df(2);
      csf_true(csf_true == idx) = wghts_df(3);
  end
end
wf_true = 0;

% iterate for different long TE times
TE_short = 10;                    % ms
TE_long = [30 60 150]';           % ms
T2_ml = NaN(nx, ny, length(TE_long));

for idx_TE = 1:length(TE_long)
  %% dictionary creation (could be precomputed)
  % SE imaging parameters
  TE = [TE_short TE_long(idx_TE)]';       % ms
  nTE = length(TE);
  TR = 3000 * ones(nTE, 1);       % ms
  flip_ex = pi/2;                 % rad   (nominal)
  flip_ref = pi;                  % rad   (nominal)
  dw = 0;                         % kHz   (no off-resonance)
  is_mag = 1;

  % to test model mismatch, make dictionary with 'old' model
  T2 = logspace(1, 3, 10000);
  D = NaN(nTE, length(T2));
  K = length(T2);

  for t2 = 1:length(T2)
    for i = 1:nTE
      D(i, t2) = SE_fun_v1(1, T2(t2), 1, flip_ex, TE(i), is_mag);
    end
  end

  % Create a true excitation profile map
  % Cannot allow too much variation in flip because then ML estimate biased
  mask = M0_true == 0;
  perc_var = 0.00;                % allow +/-0% variation in flip
  [xx, yy] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
  kap_true = (1+perc_var/2) - perc_var*(xx.^2 + yy.^2); 
  kap_true(mask) = 0;             % No excitation where M0 = 0;
  % figure, im('notick', fliplr(kap_true), [0.9 1.1], 'cbar', ' ');


  %% synthesize SE data
  % to test model mismatch, synthesize with 'new' model
  y_true = NaN(nx, ny, nTE);
  for i = 1:nTE
    y_im_true = M0_true .* ...
     (myl_true .* exp(-TE(i) / t2_myl) + ...
      dmy_true .* exp(-TE(i) / t2_dmy) + ...
      csf_true .* exp(-TE(i) / t2_csf));
    y_true(:,:,i) = fft2(y_im_true);
  end

  % add complex white gaussian noise
  var_im = 0;
  % var_im = 1.31e-7; 
  sigma_k = sqrt(var_im * (nx*ny));
  y = y_true + sigma_k * (randn(size(y_true)) + 1i*randn(size(y_true)));
  printm('snr = %g dB', 20*log10(norm(y_true(:)) / norm(col(y_true-y))));

  % noisy data back to image domain
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

  % extract indices for maximum-likelihood single-component t2 map
  T2_ml(:,:,idx_TE) = reshape(T2(t2_idx), [nx ny]);
end

%% images, rmse, and comparisons
% noise scaling of display
eps_scl = 0.2;

% postprocessing for display
T2_ml(repmat(mask, [1 1 length(TE_long)])) = 0;

% t2 images
figure; im('notick', 'row', 1, T2_ml, [20 120], 'cbar', ' ');
figure; im('notick', 'row', 1, T2eff_true, [20 120], 'cbar', ' ');

figure; im('notick', 'row', 1, abs(T2_ml-repmat(T2eff_true, [1 1 length(TE_long)])),...
  [0 eps_scl*120], 'cbar', ' ');

% compute nrmse and rmse
nrmse = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ sum(tru(:).^2));
rmse  = @(tst,tru) sqrt(sum((abs(tst(:) - tru(:))).^2) ./ numel(tru));

% gm and wm masks for task-based performance assessment
mask_gm = labels(:,:,slice) == 2;
mask_wm = labels(:,:,slice) == 3;

bias_t2_ml_gm   = NaN(length(TE_long), 1);
bias_t2_ml_wm   = NaN(length(TE_long), 1);
rmse_t2_ml_gm   = NaN(length(TE_long), 1);
rmse_t2_ml_wm   = NaN(length(TE_long), 1);
t2_ml_gm        = NaN(length(TE_long), 2);
t2_ml_wm        = NaN(length(TE_long), 2);

for idx_TE = 1:length(TE_long)
    T2_ml_curr = squeeze(T2_ml(:,:,idx_TE));
    
    % t2 bias values
    bias_t2_ml_gm(idx_TE)   = mean(T2_ml_curr(mask_gm) - T2eff_true(mask_gm));
    bias_t2_ml_wm(idx_TE)   = mean(T2_ml_curr(mask_wm) - T2eff_true(mask_wm));

    % t2 rmse values
    rmse_t2_ml_gm(idx_TE)   = rmse(T2_ml_curr(mask_gm),  T2eff_true(mask_gm));
    rmse_t2_ml_wm(idx_TE)   = rmse(T2_ml_curr(mask_wm),  T2eff_true(mask_wm));

    % t2 means and standard deviations
    t2_ml_gm(idx_TE, :)     = [mean(T2_ml_curr(mask_gm)),   std(T2_ml_curr(mask_gm))];
    t2_ml_wm(idx_TE, :)     = [mean(T2_ml_curr(mask_wm)),   std(T2_ml_curr(mask_wm))];
end

% plot bias vs. long TE
figure; plot(TE_long, bias_t2_ml_gm, TE_long, bias_t2_ml_wm, 'linewidth', 2);
xlabel('Long TE time (ms)');
ylabel('Bias, w.r.t. t2 ground truth (ms) from volumetric voxel r2 mean');
legend('GM', 'WM');

