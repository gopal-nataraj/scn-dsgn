%% Header Script for Regularized M0s/T1/T2/kappa Recon
%       03.23.2015 -- Brain DESS data; patient: J. A. Fessler
%       03.24.2015 -- Data extracted and analyzed

%% Raw data extraction
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    cd ../IRT; setup(); cd ../Scripts;
end
addpath('Calculations/DESS/');
addpath('Helper/');
addpath('JFNielsen/common');
addpath('JFNielsen/img');

% Imaging Parameters
nx = 240; ny = 240; nz = 6; nc = 8;     
opflip = [5:5:90]';                 % degrees
flip_d = opflip * pi/180;           % radians
nfd = length(flip_d);
TRd = 20.8 * ones(nfd, 1);          % ms
TEd = 6.0 * ones(nfd, 1);           % ms
wf = 0;                             % Zero off-resonance assumed
ncycles = 2;                        % Spoiling cycles across voxel
sl = 4;                             % Central slice chosen
ROI = 'cerebrum';                   % Choose cerebrum/cerebellum
pr = 0;                             % Toggle printing on/off
all_flips = 1;                      % Use all 18 flip angles 

% Get the coil-combined image data
if (~exist('yp_im', 'var'))
    try
        % First try to load DESS data from .mat files
        addpath('Brain_03,23,15/');
        load(sprintf('ims_coil_comb_%s_slice%u.mat', ROI, sl));
    catch
        dir = pwd; 
        cd('/Volumes/a2/gnataraj/Data');
        cd(strcat('2015,03,23_brain_8ch/', ROI));
        
        % Choose just one mediocoronal slice
        yp_coil = NaN(nx, ny, nfd, nc);
        ym_coil = NaN(nx, ny, nfd, nc);
        for l = 1:nfd
            load(strcat('dess_23Mar2015_', ROI, sprintf('_frame%u', l)));
            yp_coil(:,:,l,:) = ims_echo1(:,:,sl,:);         % [nx ny nc]
            ym_coil(:,:,l,:) = ims_echo2(:,:,sl,:);         % [nx ny nc]
        end
        cd(dir);
        
        % Reorder the flip angles
        % Before, they were [5, 90, 15, 80, 25, 70,...]
        yp_coil(:,:,2:2:end,:) = flipdim(yp_coil(:,:,2:2:end,:), 3);
        ym_coil(:,:,2:2:end,:) = flipdim(ym_coil(:,:,2:2:end,:), 3);

        % Take transpose to orient anterior --> posterior
        % Flip image to compensate second echo inverted in readout direction
        % This is now coil-combined complex data
        yp_coil = permute(yp_coil, [2 1 3 4]);
        ym_coil = permute(flipdim(ym_coil, 1), [2 1 3 4]);
        
        % Use mri_multidata_coil_combine.m to merge nc-channel coil data
        y_coil = cat(3, yp_coil, ym_coil);                  % [nx ny M nc]
        [y_im, smap] = mri_multidata_coil_combine(y_coil);  % [nx ny M]
        yp_im = y_im(:,:,1:nfd);                            % [nx ny nfd]
        ym_im = y_im(:,:,nfd+1:end);                        % [nx ny nfd]
        
        % Save coil-combined data for future use
        save(sprintf('ims_coil_comb_%s_slice%u.mat', ROI, sl),...
            'yp_im', 'ym_im', 'smap');
    end
end

% Even though the data is complex, take the magnitude and turn is_mag on
%   1) Large linear phase observed in main field direction; as this phase
%      is not the same for both echoes, would cause phase errors. 
%   2) Susceptibility effects could cause large off-resonance, distorting phase
yp_im = abs(yp_im);
ym_im = abs(ym_im);
is_mag = 1;

% OPTIONAL: Select only a few flip angles
if (~all_flips)
    opflip = [10 40]';                  % degrees
    flip_d = opflip  * pi/180;          % radians
    nfd = length(flip_d);
    yp_im = cat(3, yp_im(:,:,2), yp_im(:,:,8));
    ym_im = cat(3, ym_im(:,:,2), ym_im(:,:,8));
end

% Total number of datasets
M = 2*nfd;

%% Dictionary creation (could be precomputed)
T1 = logspace(1, 3.5, 100);
T2 = logspace(0, 3.5, 100);
kap = 2 .^ linspace(-1, 1, 21);
D = NaN(M, length(T1), length(T2), length(kap));
K = length(T1) * length(T2) * length(kap);

for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        for s = 1:length(kap)
            % DESS dictionary component
            for a = 1:nfd
                [D(a, t1, t2, s), D(nfd+a, t1, t2, s)] = ...
                    dess_fun(1, T1(t1), T2(t2), kap(s), flip_d(a),...
                    TRd(a), TEd(a), wf, is_mag);
            end
        end
    end
end
D = reshape(D, [M K]);

%% Maximum-likelihood Estimation
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

% Create tight and loose masks
tight_mask = imfill(squeeze(abs(yp_im(:,:,ceil(nfd/2))) >...
    median(col(abs(yp_im(:,:,ceil(nfd/2)))))), 'holes');
loose_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 15)); 

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
n_outer = 10;
niterM = 50;   
niter1 = 100;   
niter2 = 100;   
niterK = 50;

tolM = 10^-6;   
tol1 = 10^-7;  
tol2 = 10^-7;  
tolk = 10^-4; 
disp = 1; 

% Define regularizers, Rm, R1, R2, and Rk
betaM = nfd * 2^-9;  %2^-5; 
beta1 = nfd * 2^-13; %2^-9; 
beta2 = nfd * 2^-14; %2^-10;
betaK = nfd * 2^-4;  %2^-3;

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
M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  kap_ml(~reg_mask) = 0;
M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; kap_med(~reg_mask) = 0;
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0; kap_reg(~reg_mask) = 0;

%% For comparison, Method-of-Moments T2 Estimation (biased)
E2sqrd = NaN(nx, ny);
for xx = 1:nx
    for yy = 1:ny
        Y = squeeze(ym_im(xx, yy, :));
        A = squeeze(yp_im(xx, yy, :));
        C = [A Y];
        sig_small = min(svds(C'*C));
        E2sqrd(xx, yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
    end
end
T2_mom = -2*(mean(TRd - TEd)) ./ log(E2sqrd);
T2_mom(~reg_mask) = 0;

%% Summary Statistics 
% Define ROIs (WM, GM, CSF)
if strcmp(ROI, 'cerebrum')
    ctrX = [145 82  102];
    ctrY = [35  71  56 ];
    rad  = [5   2   4  ];
else fprintf('Unknown ROI');
end
nROI = length(rad);

% Extract data means and standard deviations from ROI
yp_mean = NaN(nROI, nfd); yp_std = NaN(nROI, nfd);
ym_mean = NaN(nROI, nfd); ym_std = NaN(nROI, nfd);

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
    for l = 1:nfd
        [yp_mean(r,l), yp_std(r,l)] ...
            = multiMeans(yp_im(:,:,l), [ctrX(r) ctrY(r)], rad(r));
        [ym_mean(r,l), ym_std(r,l)] ...
            = multiMeans(ym_im(:,:,l), [ctrX(r) ctrY(r)], rad(r));
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
    fid = fopen(sprintf('NIST_03,17,15_summary_%uflips', nfd), 'a');
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
if (pr), print('-deps','M0s_ml_15,40.eps'); end
figure; im('notick', abs(M0s_med), [0 4], 'cbar', ' ');
if (pr), print('-deps','M0s_med_15,40.eps'); end
figure; im('notick', abs(M0s_reg), [0 4], 'cbar', ' ');
if (pr), print('-deps','M0s_reg_15,40.eps'); end

% T1 Images
figure; im('notick', T1_ml,  [0 3000], 'cbar', ' ');
if (pr), print('-deps','T1_ml_15,40.eps'); end
figure; im('notick', T1_med, [0 3000], 'cbar', ' ');
if (pr), print('-deps','T1_med_15,40.eps'); end
figure; im('notick', T1_reg, [0 3000], 'cbar', ' ');
if (pr), print('-deps','T1_reg_15,40.eps'); end

% T2 Images w/ ROIs
t2min = 30; t2max = 80;
figure; im('notick', T2_mom, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 0], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 1 0], 'LineWidth', 2);...
    rectangle('Position', posCSF, 'EdgeColor', [0 0 1], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_mom_15,40.eps'); end
figure; im('notick', T2_ml,  [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 0], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 1 0], 'LineWidth', 2);...
    rectangle('Position', posCSF, 'EdgeColor', [0 0 1], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_ml_15,40.eps'); end 
figure; im('notick', T2_med, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 0], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 1 0], 'LineWidth', 2);...
    rectangle('Position', posCSF, 'EdgeColor', [0 0 1], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_med_15,40.eps'); end
figure; im('notick', T2_reg, [t2min t2max], 'cbar', ' '); hold on;...
    rectangle('Position', posWM,  'EdgeColor', [1 0 0], 'LineWidth', 2);...
    rectangle('Position', posGM,  'EdgeColor', [0 1 0], 'LineWidth', 2);...
    rectangle('Position', posCSF, 'EdgeColor', [0 0 1], 'LineWidth', 2); hold off;
if (pr), print('-depsc','T2_reg_15,40.eps'); end

% kap Images
figure; im('notick', kap_ml,  [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_ml_15,40.eps'); end
figure; im('notick', kap_med, [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_med_15,40.eps'); end
figure; im('notick', kap_reg, [0 2], 'cbar', ' ');
if (pr), print('-deps','kap_reg_15,40.eps'); end

% Cost vs. Iteration
figure; hold on;
scatter(1:4:4*n_outer, cost(2:4:end), 'bo');
scatter(2:4:4*n_outer, cost(3:4:end), 'ro');
scatter(3:4:4*n_outer, cost(4:4:end), 'go');
scatter(4:4:4*n_outer, cost(5:4:end), 'mo');
plot(0:4*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M0s update', 'T1 update', 'T2 update', 'kap update');
%print('-depsc','cost_vs_iter.eps');

%% Sanity check: does the data fit the model?
% Note that we can only do this because TRs are fixed across scans
% Use average values in ROIs to get corresponding ML parameters
% This prevents data noise from being amplified in estimating curves
y_mean = [yp_mean ym_mean]'; 
maxProd_mean = zeros(1, nROI);
idx_mean = zeros(1, nROI);
for k = 1:K
    % Compute kth inner product
    hess_mean = abs(D(:,k)' * W * D(:,k));
    ytild_mean = D(:,k)' * W * y_mean / sqrt(hess_mean);
    newProd_mean = abs(ytild_mean).^2;
    
    % If the kth inner product is largest, save k
    update_mean = newProd_mean > maxProd_mean;
    maxProd_mean(update_mean) = newProd_mean(update_mean);
    idx_mean(update_mean) = k;
end

% Extract indices for method-of-moments maps
[t1_idx_mean, t2_idx_mean, s_idx_mean] ...
    = ind2sub([length(T1) length(T2) length(kap)], idx_mean);
T1_sig_mean = T1(t1_idx_mean);
T2_sig_mean = T2(t2_idx_mean);
kap_sig_mean = kap(s_idx_mean);

% M0s initial guess (from dictionary)
M0s_sig_mean = NaN(nROI, 1);
for q = 1:length(idx_mean)
    M0s_sig_mean(q) = (D(:,idx_mean(q))' * y_mean(:,q)) ./...
        (D(:,idx_mean(q))' * D(:,idx_mean(q)));
end

% Signal model vs. Data
for r = 1:nROI
    flip_model = [0:60] * (pi/180) * kap_sig_mean(r);
    flip_data = flip_d * kap_sig_mean(r);
    [S1mean, S2mean] = dess_fun(M0s_sig_mean(r), T1_sig_mean(r),...
        T2_sig_mean(r), 1, flip_model, mean(TRd), mean(TEd), wf, 1);
    
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
    if (pr), print('-depsc', sprintf('sig_vs_%uflips_ROI%02u.eps', nfd, r)); end;
end