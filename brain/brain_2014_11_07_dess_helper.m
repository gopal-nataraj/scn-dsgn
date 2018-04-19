%% Subroutine: T1/T2 ML and Reg estimates at different flip angle combos

%% Method-of-Moments Estimation
% Dictionary-based estimation via variable-projection method
weights = ones(M, 1);
W = spdiags(weights, 0, M, M);      % Weighting Matrix
y = reshape(permute(cat(3, yp_im, ym_im), [3 1 2]), [M nx*ny]);

maxProd = zeros(1, nx*ny);
idx = zeros(1, nx*ny);
for k = 1:K
    % Compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;

    % If the kth inner product is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    idx(update) = k;
end

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

% % Create tight and loose masks
% tight_mask = imfill(squeeze(abs(yp_im(:,:,ceil(nfd/2)+1)) >...
%     0.1*max(col(abs(yp_im(:,:,ceil(nfd/2)+1))))), 'holes');
% loose_mask = imdilate(tight_mask, strel('disk', 10)); 

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

% Median filtering (M0s_med might not be right...)
T1_med = medfilt2(T1_ml);  
T2_med = medfilt2(T2_ml); 
kap_med = medfilt2(kap_ml);

%% Regularized, Joint M0s,T1,T2,kap Estimation
% Define iteration parameters
% n_outer = 10;    
% niterM = 50;   
% niter1 = 100;   
% niter2 = 100;   
% niterK = 50;
% 
% tolM = 10^-6;   
% tol1 = 10^-7;  
% tol2 = 10^-7;  
% tolk = 10^-4; 
disp = 1; 

% Define regularizers, Rm, R1, and R2
betaM = nfd * 2^-8; 
beta1 = nfd * 2^-15;%2^-14; 
beta2 = nfd * 2^-15;%2^-15;
betaK = nfd * 2^-4;

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
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);

% Regularized Reconstruction
tic; 
[M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
    mri_spgr_dess_m0st1t2kap_map(T1_med, T2_med, kap_med, wf_ml, flip_s,...
    flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
    TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
    niter2, niterK, tolM, tol1, tol2, tolk, is_mag, disp);
recon_time = toc; 

% Postprocessing for display
% For real data, can use this mask
% reg_mask = imfill(abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:))), 'holes');

% % Remove pixels outside reg_mask
M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  kap_ml(~reg_mask) = 0;
M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; kap_med(~reg_mask) = 0;
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0; kap_reg(~reg_mask) = 0;
