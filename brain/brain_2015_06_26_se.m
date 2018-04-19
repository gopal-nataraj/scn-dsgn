%% Header Script for Regularized M0/T2 Estimation from 2D Single SE Data at Variable TE
%       05.15.2015 -- Pulse sequences prepared
%       06.26.2015 -- Brain Data acquired (J. A. Fessler)
%       08.04.2015 -- Choose reg parameters carefully
%                     Matched dictionary T1/T2 ranges to NIST 06/11 data
%       05.13.2016 -- Added manual roi selection capability
%       05.16.2016 -- Manual registration of ROI to fast scans
%
% Collected at two optimized TE times
% Written by: Gopal Nataraj

%% Raw data extraction
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    curdir = pwd; cd ../../../irt; setup(); cd(curdir);
end
addpath('../../data/Brain_06,26,15/');
addpath('../../model/se/');
addpath('../../map/b1/');
addpath('../../map/sense/');
addpath('../../map/t1-t2/old/');

% Imaging parameters
nx = 256; ny = 256; nz = 1; nc = 32;        
FOVx = 240; FOVy = 240; FOVz = 5;           % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

TE = [10 150]';                             % ms
nTE = length(TE);                           
TR = 3000 * ones(nTE, 1);                   % ms
flip_ex = pi/2;                             % rad   Nominal excitation flip
wf = 0;                                     % Hz    Assume zero off-resonance
kap = 1;                                    %       Flip angle scaling
% rs = 0;                                     %       Toggle ROI selection prompts on/off
% pr = 0;                                     % Toggle printing on/off


%% Coil-combine the SE raw data
if (~exist('y_im', 'var'))
    try
        % First try to load SE coil-combined data
        addpath('../../data/Brain_06,26,15/SE/');
        load('ims_coil_comb_SE.mat');
    catch
        % Grab the SE coil data
        dir = pwd; cd('/Volumes/a2/gnataraj/Data/2015,06,26_brain_32ch/SE');
        load('im_se_coil.mat'); cd(dir);
        
        % Use mri_multidata_coil_combine.m to merge nc-channel coil data
        y_coil = permute(im_se, [2 1 4 3]);                 % [nx ny nTE nc]
        [y_im, smap] = mri_multidata_coil_combine(y_coil);  % [nx ny nTE]
        
        % Save coil-combined data for future use
        dir = pwd; cd('../../data/Brain_06,26,15/SE/');
        save('ims_coil_comb_SE.mat', 'y_im', 'smap'); cd(dir);
    end
end

% Even though the data is complex, take the magnitude and turn is_mag on
% This is to remain consistent with the SPGR and DESS recons
y_im = abs(y_im);
is_mag = 1;

% Trick: normalize data by median of non-background data values
% so that the effective regularization beta is scale-invariant
tmp = abs(y_im);
tmp = median(tmp(tmp > 0.1 * max(tmp(:))));
y_im = y_im / tmp;

% Crop images for display
x_crop  = [38:217];  nx = length(x_crop);
% y_crop  = [8:227];   ny = length(y_crop);
y_crop = [9:228];  ny = length(y_crop);
y_im   = y_im(x_crop, y_crop, :);

% Create tight and loose masks
tight_mask = imfill(squeeze(abs(y_im(:,:,ceil(nTE/2))) >...
    0.02*max(col(abs(y_im(:,:,ceil(nTE/2)))))), 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 10));


%% Dictionary creation (could be precomputed)
% T2 = logspace(0, 3.5, 1000);
T2 = logspace(1, 3, 300);
D = NaN(nTE, length(T2));
K = length(T2);

for t2 = 1:length(T2)
    for i = 1:nTE
        D(i, t2) = SE_fun_v1(1, T2(t2), kap, flip_ex, TE(i), is_mag);
    end
end


%% Maximum-likelihood Estimation
% Dictionary-based estimation via variable-projection method
weights = ones(nTE, 1);
W = spdiags(weights, 0, nTE, nTE);      % Weighting matrix
y = reshape(permute(y_im, [3 1 2]), [nTE nx*ny]);

tic;
maxProd = zeros(1, nx*ny);
t2_idx = ones(1, nx*ny);
for k = 1:K
    % Compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;

    % If the kth inner product is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    t2_idx(update) = k;
end
time_ml = toc;

% Extract indices for maximum-likelihood maps
T2_ml = reshape(T2(t2_idx), [nx ny]);

% Assume negligible off-resonance effects (short TE) and kap=1
wf_ml = zeros(nx, ny);
kap_ml = ones(nx, ny);

% % M0 initial guess (from dictionary)
% M0_ml = NaN(nx*ny, 1);
% for q = 1:length(t2_idx)
%     M0_ml(q) = (D(:,t2_idx(q))' * y(:,q)) ./ (D(:,t2_idx(q))' * D(:,t2_idx(q)));
% end
% M0_ml = reshape(M0_ml, [nx ny]);


%% Preprocessing and Masking
% Project Images to within range
T2max = 2000;       T2_ml = min(T2_ml, T2max);
T2min = 1;          T2_ml = max(T2_ml, T2min);

% Set voxels inside loose mask but outside tight mask to mean
T2_ml(~tight_mask & loose_mask) = mean(col(T2_ml));  

% Set voxels outside both tight and loose masks to zero
T2_ml(~loose_mask) = 0;  
wf_ml(~loose_mask) = 0;  

% Median filtering
T2_med = medfilt2(T2_ml);


%% Regularized, Joint M0, T2 Estimation
% Define iteration parameters
n_outer = 100;
niterM = 100;
niter2 = 100;

tolM = 10^-7;
tol2 = 10^-7;
disp = 0;

% Define regularizers, Rm and R2
betaM = nTE * 10^-9;                        % nTE * 2^-10;  
beta2 = nTE * 10^-8;                        % nTE * 2^-10;   

% Delta values on the order of one standard deviation
deltaM = 10^-0.5;                           % 10^4;
delta2 = 10^0.5;                            % 10^-2;

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Reshape data inputs to work with 3D implementation
y_3D = permute(y_im, [1 2 4 3]);

% Regularized reconstruction
tic; 
% [M0_ml, M0_med, M0_reg, T2_reg, cost] = ...
%     mri_se_m0t2_map_v1(T2_med, y_3D, TE, loose_mask, T2max, T2min,...
%     Rm, R2, n_outer, niterM, niter2, tolM, tol2, disp);
[M0_ml, M0_med, M0_reg, T2_reg, cost] = ...
    mri_se_m0t2_map_v2(T2_med, kap_ml, y_3D, TE, loose_mask, T2max, T2min,...
    Rm, R2, n_outer, niterM, niter2, tolM, tol2, is_mag, disp);
time_reg = toc;

% Postprocessing for display
% For real data, can use this mask
% reg_mask = imfill(abs(M0_reg) >= 0.1*abs(max(M0_reg(:))), 'holes');
reg_mask = tight_mask;

% % Remove pixels outside reg_mask
M0_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  
M0_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; 
M0_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0; 


%% Images, RMSE, and Comparisons
suffix = sprintf('%uSE_brain', nTE);

% Define ROI masks
if rs
    prompt = 'ROI selection: [d]efault, [m]anual, or [s]aved? ';
    str = input(prompt, 's');
else
    str = 's';
end
roi_labels = {'WM','GM','SF'};
nROI = size(roi_labels, 2);

switch str
    case 'd'
        % Default circular ROIs
        ctrX = [109 142 101];
        ctrY = [47  157 69];
        rad  = [4   2    2];

        % ROI-only mask
        roi_mask = false(nx, ny, nROI);
        [xx, yy] = ndgrid(1:nx, 1:ny);
        for r = 1:nROI
            tmp = false(nx, ny);
            tmp((xx-ctrX(r)).^2 + (yy-ctrY(r)).^2 <= rad(r)^2) = true;
            roi_mask(:,:,r) = tmp;
        end
    
    case 'm'
        roi_mask = NaN(nx, ny, nROI);
        for r = 1:nROI
            % ROI manual selection
            fprintf('Select %s polygonal ROI.\n', roi_labels{r});
            roi_mask(:,:,r) = roipoly(T2_ml/(T2max/10));
        end
        roi_mask = roi_mask > 0;
        
        % Save for future use
        dir = pwd; 
        cd('Brain_06,26,15');
        if exist('roi_masks.mat', 'file')
            prompt = 'Overwrite previously definely roi_mask [y/n]? ';
            str2 = input(prompt, 's');
            switch str2
                case 'y'
                    save('roi_masks.mat', 'roi_mask'); 
                    fprintf('roi_masks.mat overwritten.\n');
                case 'n'
                    fprintf('New ROI mask not saved.\n');
                otherwise
                    error('Unrecognized input.');
            end
        else
            save('roi_masks.mat', 'roi_mask');
        end
        cd(dir);
        
    case 's'
        try
            addpath('Brain_06,26,15');
            load('roi_masks.mat');
        catch
            error('ROI masks not found.');
        end
        
    otherwise
        error('Unrecognized input.');
end

figure;
im(roi_mask);

% ROI boundaries
roi_bound.wm = bwboundaries(roi_mask(:,:,1), 'noholes');
roi_bound.gm = bwboundaries(roi_mask(:,:,2), 'noholes');

% M0 Images 
m0min = 0; m0max = 16; m0_rng = [m0min m0max];
figure; im('notick', abs(M0_ml), m0_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0_ml_', suffix, '.eps')), end;
figure; im('notick', abs(M0_med), m0_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0_med_', suffix, '.eps')), end;
figure; im('notick', abs(M0_reg), m0_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0_reg_', suffix, '.eps')), end;

% T2 Images
t2min = 20; t2max = 120; t2_rng = [t2min t2max]; 
figure; im('notick', T2_ml, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_ml_', suffix, '.eps')), end;
figure; im('notick', T2_med, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_med_', suffix, '.eps')), end;
figure; im('notick', T2_reg, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_reg_', suffix, '.eps')), end;

% T2 Images w/ ROIs
figure; im('notick', T2_ml, t2_rng, 'cbar', ' '); hold on;
    for i = 1:length(roi_bound.wm)
        bound_wm = roi_bound.wm{i};
        plot(bound_wm(:,1), bound_wm(:,2), 'm', 'LineWidth', 1);
    end
    for i = 1:length(roi_bound.gm)
        bound_gm = roi_bound.gm{i};
        plot(bound_gm(:,1), bound_gm(:,2), 'c', 'LineWidth', 1);
    end
    hold off;
if (pr), print('-depsc', strcat('T2_ml_ROI_', suffix, '.eps')), end;

% Cost vs. Iteration
figure; hold on;
scatter(1:2:2*n_outer, cost(2:2:end), 'bo');
scatter(2:2:2*n_outer, cost(3:2:end), 'ro');
plot(0:2*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M0 update', 'T2 update');
%print('-depsc','cost_vs_iter.eps');


%% Sanity check: does the data fit the model?
% Note: plots only make sense if TRs are fixed across scans

% Extract data, M0, T2 means and standard deviations from ROIs
y_mean = NaN(nROI, nTE); y_std = NaN(nROI, nTE);

M0_re_ml_mean = NaN(nROI, 1);   M0_re_ml_std = NaN(nROI, 1);
M0_im_ml_mean = NaN(nROI, 1);   M0_im_ml_std = NaN(nROI, 1);
M0_re_reg_mean = NaN(nROI, 1);  M0_re_reg_std = NaN(nROI, 1);
M0_im_reg_mean = NaN(nROI, 1);  M0_im_reg_std = NaN(nROI, 1);

T2_ml_mean = NaN(nROI, 1);      T2_ml_std = NaN(nROI, 1);
T2_reg_mean = NaN(nROI, 1);     T2_reg_std = NaN(nROI, 1);

for r = 1:nROI
    for i = 1:nTE
        y_mean(r,i)     = mean(masker(y_im(:,:,i), roi_mask(:,:,r)));
        y_std(r,i)      = std(masker(y_im(:,:,i), roi_mask(:,:,r)));
    end
    
    M0_re_ml_mean(r)   = mean(masker(real(M0_ml), roi_mask(:,:,r)));
    M0_im_ml_mean(r)   = mean(masker(imag(M0_ml), roi_mask(:,:,r)));
    M0_re_ml_std(r)    = std(masker(real(M0_ml), roi_mask(:,:,r)));
    M0_im_ml_std(r)    = std(masker(imag(M0_ml), roi_mask(:,:,r)));
    
    M0_re_reg_mean(r)  = mean(masker(real(M0_reg), roi_mask(:,:,r)));
    M0_im_reg_mean(r)  = mean(masker(imag(M0_reg), roi_mask(:,:,r)));
    M0_re_reg_std(r)   = std(masker(real(M0_reg), roi_mask(:,:,r)));
    M0_im_reg_std(r)   = std(masker(imag(M0_reg), roi_mask(:,:,r)));
    
    T2_ml_mean(r)       = mean(masker(T2_ml, roi_mask(:,:,r)));
    T2_ml_std(r)        = std(masker(T2_ml, roi_mask(:,:,r)));
    
    T2_reg_mean(r)      = mean(masker(T2_reg, roi_mask(:,:,r)));
    T2_reg_std(r)       = std(masker(T2_reg, roi_mask(:,:,r)));
    
%     for i = 1:nTE
%         [y_mean(r, i), y_std(r, i)]...
%             = multiMeans(y_im(:,:,i), [ctrX(r) ctrY(r)], rad(r));
%     end
%     
%     [M0_re_ml_mean(r), M0_re_ml_std(r)] ...
%         = multiMeans(real(M0_ml), [ctrX(r) ctrY(r)], rad(r));
%     [M0_im_ml_mean(r), M0_im_ml_std(r)] ...
%         = multiMeans(imag(M0_ml), [ctrX(r) ctrY(r)], rad(r));
%     [M0_re_reg_mean(r), M0_re_reg_std(r)] ...
%         = multiMeans(real(M0_reg), [ctrX(r) ctrY(r)], rad(r));
%     [M0_im_reg_mean(r), M0_im_reg_std(r)] ...
%         = multiMeans(imag(M0_reg), [ctrX(r) ctrY(r)], rad(r));
%     
%     [T2_ml_mean(r), T2_ml_std(r)] ...
%         = multiMeans(T2_ml, [ctrX(r) ctrY(r)], rad(r));
%     [T2_reg_mean(r), T2_reg_std(r)] ...
%         = multiMeans(T2_reg, [ctrX(r) ctrY(r)], rad(r));
end

% Export summary statistics to file
if (pr)
    fid = fopen(strcat('Summary_06,26,15_', suffix), 'w');
    fprintf(fid, sprintf('Estimation Statistics for 2 SE Scans'));
    fprintf(fid, sprintf('\n\tML run time: %0.2f seconds', time_ml));
    fprintf(fid, sprintf('\n\tReg run time: %0.2f seconds', time_reg));
    
    fprintf(fid, '\n\nRe[M0] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M0_re_ml_mean(r),   char(177),  M0_re_ml_std(r),...
            M0_re_reg_mean(r),  char(177),  M0_re_reg_std(r));
    end
    
    fprintf(fid, '\n\nIm[M0] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M0_im_ml_mean(r),   char(177),  M0_im_ml_std(r),...
            M0_im_reg_mean(r),  char(177),  M0_im_reg_std(r));
    end

    fprintf(fid, '\n\nT2 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            T2_ml_mean(r), char(177), T2_ml_std(r),...
            T2_reg_mean(r), char(177), T2_reg_std(r));
    end
    fclose(fid);
end

% Use average values in ROIs to get corresponding ML parameters
% This prevents data noise from being amplified in estimating curves
maxProd_mean = zeros(1, nROI);
t2_idx_mean = zeros(1, nROI);
for k = 1:K
    % Compute kth inner product
    hess_mean = abs(D(:,k)' * W * D(:,k));
    ytild_mean = D(:,k)' * W * y_mean' / sqrt(hess_mean);
    newProd_mean = abs(ytild_mean).^2;
    
    % If the kth inner product is largest, save k
    update_mean = newProd_mean > maxProd_mean;
    maxProd_mean(update_mean) = newProd_mean(update_mean);
    t2_idx_mean(update_mean) = k;
end

% Extract indices for method-of-moments maps
T2_sig_mean = T2(t2_idx_mean);

% M0 initial guess (from dictionary)
M0_sig_mean = NaN(nROI, 1);
for q = 1:length(t2_idx_mean)
    M0_sig_mean(q) = (D(:,t2_idx_mean(q))' * y_mean(q,:)') ./...
        (D(:,t2_idx_mean(q))' * D(:,t2_idx_mean(q)));
end

% Signal model vs. Data
% Note: no flip angle compensation
TE_model = 0:200;
Smean = NaN(nROI, length(TE_model));
for r = 1:nROI
    Smean(r,:) = SE_fun_v1(M0_sig_mean(r), T2_sig_mean(r), kap, flip_ex, TE_model, is_mag);
end

figure; hold on;
ColorOrd = get(gca, 'ColorOrder');
for r = 1:nROI
    errorbar(TE, abs(y_mean(r,:)), abs(y_std(r,:)), '*', 'Color', ColorOrd(r,:));
    plot(TE_model, Smean(r,:), 'Color', ColorOrd(r,:));
end
title('GE 3T 2D SE Data (d) vs. Model (m)');
xlabel('Echo time (TE)'); ylabel('Magnitude Signal (a.u.)');
legend('WM-d', 'WM-m', 'GM-d', 'GM-m', 'CSF-d', 'CSF-m', 'Location', 'NE');
if (pr), print('-depsc', strcat('sig_', suffix)); end;