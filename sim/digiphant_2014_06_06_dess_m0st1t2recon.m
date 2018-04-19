%% Header Script for Regularized, Joint M0s, T1, T2 Recon from Synthetic DESS Data

%% Synthesis of Phantom DESS Data
% Load digital phantom
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;
addpath('Calculations/DESS/');

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

% Imaging Parameters
flip = [11 37]' * pi/180;               % radians
%flip = [17 43]' * pi/180;               % radians
nf = length(flip);
TR = 20;                                % ms
TE = 5;                                 % ms
wf_true = 0;                            % Off-resonance
SNR = 40;                               % dB
E1_true = exp(-TR./T1_true);
E2_true = exp(-TR./T2_true);

% Dictionary creation (could be precomputed)
T1 = logspace(1, 4, 50);
T2 = logspace(0, 3, 50);
D = NaN(2*length(flip), length(T1), length(T2));
for t1 = 1:length(T1)
    E1 = exp(-TR / T1(t1));
    for t2 = 1:length(T2)
        for a = 1:nf
            v1 = (1 - E1*cos(flip(a))) ./ (E1 - cos(flip(a)));
            [D(a, t1, t2), D(nf+a, t1, t2)] = ...
                 dess_fun_M0star(TR, TE, 1, wf_true, flip(a), v1, T2(t2));
        end
    end
end
D = reshape(D, [2*length(flip) length(T1)*length(T2)]);

% Make a M0s_true map, where M0s = M0 * exp(-TE/T2s)
% This is to be used for generating the forward model 
T2s_msk = T2_true ~= T2s_true; 
M0s_true = M0_true;
M0s_true(T2s_msk) = M0_true(T2s_msk) .* exp(-TE ./ T2s_true(T2s_msk));

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
% figure; im('notick', fliplr3D(abs(yp_im)), [0 0.2], 'cbar', ' '); print -deps yp_im.eps
% figure; im('notick', fliplr3D(abs(ym_im)), [0 0.2], 'cbar', ' '); print -deps ym_im.eps

%% Method-of-Moments Estimation
% Dictionary-based estimation via variable-projection method
W = diag([1 1 1 1]);            % Weighting matrix
y = reshape(permute(cat(3, yp_im, ym_im), [3 1 2]), [2*nf nx*ny]);

inProdTbl = NaN(length(T1)*length(T2), nx*ny);
for k = 1:(length(T1)*length(T2))
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    inProdTbl(k,:) = sum(abs(ytild.^2), 1);
end

[~, idx] = max(inProdTbl, [], 1);
[t1_idx, t2_idx] = ind2sub([length(T1) length(T2)], idx);
T1_mom = reshape(T1(t1_idx), [nx ny]); T1_mom(~T2s_msk) = 0;
T2_mom = reshape(T2(t2_idx), [nx ny]); T2_mom(~T2s_msk) = 0;

% % Method-of-Moments T2 Estimation
% if (nf == 1)
%     T2_mom = -2*(TR-TE) ./ log(abs(ym_im ./ yp_im));
% else
%     E2squared = NaN(nx, ny);
%     for xx = 1:nx
%         for yy = 1:ny 
%             Y = squeeze(ym_im(xx,yy,:));
%             A = squeeze(yp_im(xx,yy,:));
%             C = [A Y];
%             sig_small = min(svds(C'*C));
%             E2squared(xx,yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
%         end
%     end
%     T2_mom = -2*(TR-TE) ./ log(E2squared);
% end
% 
% % Method-of-Moments T1 Estimation (Method One)
% E2E_mom = exp(-TE ./ T2_mom);
% yp_bst1 = abs(squeeze(yp_im(:,:,1))) ./ tan(flip(1)/2);
% ym_bst1 = abs(squeeze(ym_im(:,:,1))) ./ tan(flip(1)/2) .* E2E_mom.^2;
% yp_bst2 = abs(squeeze(yp_im(:,:,2))) ./ tan(flip(2)/2);
% ym_bst2 = abs(squeeze(ym_im(:,:,2))) ./ tan(flip(2)/2) .* E2E_mom.^2;
% 
% flip_avg = mean([flip(1) flip(2)]);
% v1est = (ym_bst1 - ym_bst2) ./ (yp_bst1 - yp_bst2);
% % T1_mom = abs(-TR ./ log((1+v1est*cos(flip_avg)) ./ (v1est+cos(flip_avg))));
% T1_mom = abs(-TR ./ log((1+v1est*cos(flip(1))) ./ (v1est+cos(flip(1)))));
% T1_mom(~T2s_msk) = 0; 

% Method-of-Moment T1 Estimation (Method Two)
% bst_flp_idx = 1;
% E2E_mom = exp(-TE ./ T2_mom);
% yp_bst = abs(squeeze(yp_im(:,:,bst_flp_idx)));
% ym_bst = abs(squeeze(ym_im(:,:,bst_flp_idx)));
% M0s_t1est = M0s_mom .* tan(flip(bst_flp_idx)/2);
% v1est = (1 - ym_bst .* E2E_mom.^2 ./ M0s_t1est) ./ (1 - yp_bst ./ M0s_t1est);
% T1_mom = -TR ./ log((1+v1est*cos(flip(bst_flp_idx))) ./ (v1est+cos(flip(bst_flp_idx))));
% T1_mom(isnan(T1_mom)) = 0;

% % Method-of-Moments M0s Estimation
% M0s_mom = zeros(nx, ny);                 
% M0s_mom(T2s_msk) = 1;

% For now, take true value wf=0
wf_mom = zeros(nx, ny);

%% Preprocessing and Masking
% Project Images to within range
T1max = 5000;       E1max = exp(-TR ./ T1max);
T1min = 5;          E1min = exp(-TR ./ T1min);
T2max = 500;        E2max = exp(-TR ./ T2max);
T2min = 5;          E2min = exp(-TR ./ T2min);

T1_mom = max(T1_mom, T1min);
T1_mom = min(T1_mom, T1max);
T2_mom = max(T2_mom, T2min);
T2_mom = min(T2_mom, T2max); 

% Preprocessing and Masking
tight_mask = T2s_msk;
loose_mask = imdilate(tight_mask, strel('disk', 5)); 

% Set voxels inside loose mask but outside tight mask to mean
T1_mom(~tight_mask & loose_mask) = mean(col(T1_mom));   
T2_mom(~tight_mask & loose_mask) = mean(col(T2_mom));   
wf_mom(~tight_mask & loose_mask) = mean(col(wf_mom));   

% Set voxels outside both tight and loose masks to zero
T1_mom(~loose_mask) = 0;  
T2_mom(~loose_mask) = 0;  
wf_mom(~loose_mask) = 0;  

% Median filtering (M0s_med might not be right...)
T1_med = medfilt2(T1_mom);  
T2_med = medfilt2(T2_mom);  

% Update the E1 and E2 maps
E1_mom = exp(-TR ./ T1_mom);  
E1_med = exp(-TR ./ T1_med);  
E2_mom = exp(-TR ./ T2_mom);  
E2_med = exp(-TR ./ T2_med);  

%% Regularized, Joint M0s,T1,T2 Estimation
% Define iteration parameters
n_outer = 30;    
niterM = 200;   
niter1 = 200;   
niter2 = 200;   
tolM = 10^-5;   
tol1 = 10^-5;  
tol2 = 10^-5;  
T2meth = 'BBGM'; 
disp = 1; 

% Define regularizers, Rm, R1, and R2
betaM = 2^-6; 
beta1 = 2^-2; %2^-3; 
beta2 = 2^-6; %2^-7;

deltaM = 10^-2;
delta1 = 10^0;
delta2 = 10^-1; 

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Reshape inputs to work with 3-D implementation
flip_im = permute(repmat(flip, [1 nx ny]), [2 3 1]);
flip_3D = permute(flip_im, [1 2 4 3]);
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);

% % Perturb flip angle maps to see if recon can handle it
% perturb = 0.05;
% flip_3D = flip_3D + perturb * randn(size(flip_3D));

% Regularized reconstruction
tic;
[M0s_mom, M0s_med, M0s_reg, T1_reg, T2_reg, wf_reg, cost] = ...
    mri_dess_m0st1t2map(E1_med, E2_med, wf_mom, flip_3D, yp_3D, ym_3D,...
    loose_mask, T1max, T1min, T2max, T2min, TR, TE, Rm, R1, R2,...
    n_outer, niterM, niter1, niter2, tolM, tol1, tol2, T2meth, disp);
time_joint = toc;

% Postprocessing for display
% For real data, can use this mask
reg_mask = abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:)));

% Remove pixels outside reg_mask
M0s_mom(~reg_mask) = 0; T1_mom(~reg_mask) = 0; T2_mom(~reg_mask) = 0;
M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0;
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0;

%% Images, RMSE, and Comparisons

% Noise scaling of display
eps_scl = 0.2;

% M0s Images 
% figure; im('notick', fliplr(abs(M0s_mom)), [0 1], 'cbar', ' ');
% print -deps M0s_mom_dess.eps;
% figure; im('notick', fliplr(abs(M0s_med)), [0 1], 'cbar', ' ');
% print -deps M0s_med_dess.eps;
figure; im('notick', fliplr(abs(M0s_reg)), [0 1], 'cbar', ' ');
% print -deps M0s_reg_dess.eps;
% figure; im('notick', fliplr(abs(M0s_true)), [0 1], 'cbar', ' ');
% print -deps M0s_tru.eps;

% figure; im('notick', fliplr(abs(M0s_mom-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_mom_err_dess.eps;
% figure; im('notick', fliplr(abs(M0s_med-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_med_err_dess.eps;
figure; im('notick', fliplr(abs(M0s_reg-M0s_true)), [0 eps_scl*1], 'cbar', ' ');
% print -deps M0s_reg_err_dess.eps;

% T1 Images
T1_med(~loose_mask) = 0;
% figure; im('notick', fliplr(T1_mom), [0 2000], 'cbar', ' ');
% print -deps T1_mom_dess.eps;
% figure; im('notick', fliplr(T1_med), [0 2000], 'cbar', ' ');
% print -deps T1_med_dess.eps;
figure; im('notick', fliplr(T1_reg), [0 2000], 'cbar', ' ');
% print -deps T1_reg_dess.eps;
% figure; im('notick', fliplr(T1_true), [0 2000], 'cbar', ' ');
% print -deps T1_tru.eps;

% figure; im('notick', fliplr(T1_mom-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_mom_err_dess.eps;
% figure; im('notick', fliplr(T1_med-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_med_err_dess.eps;
figure; im('notick', fliplr(T1_reg-T1_true), [0 eps_scl*2000], 'cbar', ' ');
% print -deps T1_reg_err_dess.eps;

% T2 Images
% figure; im('notick', fliplr(T2_mom), [0 200], 'cbar', ' ');
% print -deps T2_mom_dess.eps;
% figure; im('notick', fliplr(T2_med), [0 200], 'cbar', ' ');
% print -deps T2_med_dess.eps;
figure; im('notick', fliplr(T2_reg), [0 200], 'cbar', ' ');
% print -deps T2_reg_dess.eps;
% figure; im('notick', fliplr(T2_true), [0 200], 'cbar', ' ');
% print -deps T2_tru.eps;

% figure; im('notick', fliplr(abs(T2_mom-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_mom_err_dess.eps;
% figure; im('notick', fliplr(abs(T2_med-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_med_err_dess.eps;
figure; im('notick', fliplr(abs(T2_reg-T2_true)), [0 eps_scl*200], 'cbar', ' ');
% print -deps T2_reg_err_dess.eps;

% Cost vs. Iteration
figure; hold on;
scatter([1:3:3*n_outer], cost(2:3:end), 'bo');
scatter([2:3:3*n_outer], cost(3:3:end), 'ro');
scatter([3:3:3*n_outer], cost(4:3:end), 'go');
plot([0:3*n_outer], cost, 'k'); hold off;
title('Cost vs. iteration');
legend('After M0s update','After E1 update','After E2 update');
% print -depsc cost_vs_iter.eps

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
rmse_t2_reg_gm = rmse(T2_reg(mask_gm), T2_true(mask_gm));

rmse_t2_mom_wm = rmse(T2_mom(mask_wm), T2_true(mask_wm));
rmse_t2_med_wm = rmse(T2_med(mask_wm), T2_true(mask_wm));
rmse_t2_reg_wm = rmse(T2_reg(mask_wm), T2_true(mask_wm));

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
t2_reg_gm = [mean(T2_reg(mask_gm)), std(T2_reg(mask_gm))];
t2_tru_gm = [mean(T2_true(mask_gm)), std(T2_true(mask_gm))];

t2_mom_wm = [mean(T2_mom(mask_wm)), std(T2_mom(mask_wm))];
t2_med_wm = [mean(T2_med(mask_wm)), std(T2_med(mask_wm))];
t2_reg_wm = [mean(T2_reg(mask_wm)), std(T2_reg(mask_wm))];
t2_tru_wm = [mean(T2_true(mask_wm)), std(T2_true(mask_wm))];

