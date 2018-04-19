%% Header Script for Regularized M0s/T1 Estimation from 2D SE IR Data at Variable TI
%       05.15.2015 -- Pulse sequences prepared
%       06.26.2015 -- Brain Data acquired (J. A. Fessler)
%       07.28.2015 -- BS .wav file found to be 8kHz, not 4kHz; b1 corrected
%       08.04.2015 -- Choose reg parameters carefully
%       08.05.2015 -- Added T1_ml print version w/ ROIs labeled
%                     Matched dictionary T1/T2 ranges to NIST 06/11 data
%       05.12.2016 -- Added manual roi selection capability
%                     Scaled kap_scale by fudge_factor
%       05.16.2016 -- Manual registration of ROI to fast scans
%
% Collected at two optimized TI times
% Written by: Gopal Nataraj

%% Raw data extraction
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    curdir = pwd; cd ../../../irt; setup(); cd(curdir);
end
addpath('../../data/Brain_06,26,15/');
addpath('../../model/ir/');
addpath('../../map/b1/');
addpath('../../map/sense/');
addpath('../../map/t1-t2/');

% Imaging parameters
nx = 256; ny = 256; nz = 1; nc = 32;         
FOVx = 240; FOVy = 240; FOVz = 5;           % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

TI = [650 900]';                            % ms
nTI = length(TI);                           
TR = 3000 * ones(nTI, 1);                   % ms
TE = 14 * ones(nTI, 1);                     % ms    Min allowed on scanner
flip_inv = pi;                              % rad   Nominal inversion flip
flip_ex = pi/2;                             % rad   Nominal excitation flip
wf = 0;                                     % Hz    Assume zero off-resonance
sl = 2;                                     %       Central slice (only for B1)
% rs = 0;                                     %       Toggle ROI selection prompts on/off
% pr = 0;                                     %       Toggle printing on/off


%% B1 mapping from Bloch-Siegert Data
if (~exist('kap_reg', 'var'))
    try
        addpath('../../data/Brain_06,26,15/BS');
        load('B1_reg');
    catch
        % Compute Bloch-Siegert constant
        % 07.28.2015: found that frequency offset is twice what is listed
        [~, mag_p, ph_p] = readwav('bs_8000Hz_12-Jun-2015.wav');
        [~, mag_m, ph_m] = readwav('bs_-8000Hz_12-Jun-2015.wav');
        w_rf = 8000;                            % Hz        (meas. 8.0645kHz off-res.)
        t = 4e-6 * (0:length(mag_p)-1)';        % s         (4 us sampling)
        gam = 4.258e3;                          % Hz/Gauss
        w_bs = (gam*mag_p).^2/(2*w_rf);         % Hz
        phi_bs = trapz(t, (2*pi*w_bs));         % rad
        K_bs = phi_bs / max(mag_p.^2);          % rad/Gauss^2

        % Load (3D) Bloch-Siegert Data
        % For, now only use single slice of interest
        dir = pwd; cd('/Volumes/a2/gnataraj/Data/2015,06,26_brain_32ch/BS');
        load('bs_p4000Hz.mat'); ybs_p_coil = permute(flipdim(ims(:,:,sl,:), 1), [2 1 4 3]);
        load('bs_m4000Hz.mat'); ybs_m_coil = permute(flipdim(ims(:,:,sl,:), 1), [2 1 4 3]);
        cd(dir);
        
        % Regularized B1 map estimation and flip angle correction
        [b1_init, b1_reg, ~] = mri_bsb1map(ybs_p_coil, ybs_m_coil,...
            'Kbs', K_bs, 'log2b_b1', 2, 'iter_b1', 1000);
        b1_peak = 0.075;                                % Gauss, for 90deg flip
        kap_reg = b1_reg / b1_peak;                     % Flip angle scale factor
        
        % Save regularized B1 maps for future use
        dir = pwd; cd('../../data/Brain_06,26,15/BS');
        save('B1_reg', 'b1_peak', 'b1_init', 'b1_reg', 'kap_reg'); cd(dir);
    end
    
    % DEBUG: try fudging the scale factor
%     kap_reg = ones(nx, ny);
%     fudge_factor = 1.0;
    kap_reg = fudge_factor * kap_reg;

    figure; im(b1_init, [0 b1_peak], 'cbar');
    figure, im(b1_reg,  [0 b1_peak], 'cbar');
    figure; im(kap_reg, [0 1], 'cbar');
end


%% Coil-combine the SE IR raw data
if (~exist('y_im', 'var'))
    try
        % First try to load SE IR coil-combined data
        addpath('../../data/Brain_06,26,15/IR/');
        load('ims_coil_comb_IR.mat');
    catch
        % Grab the SE IR coil data
        dir = pwd; cd('/Volumes/a2/gnataraj/Data/2015,06,26_brain_32ch/IR');
        load('im_ir_coil.mat'); cd(dir);
        
        % Use mri_multidata_coil_combine.m to merge nc-channel coil data
        y_coil = permute(im_ir, [2 1 4 3]);                 % [nx ny nTI nc]
        [y_im, smap, y_im_ssos, s_ml, coil_cost] = ...      % [nx ny nTI]
            mri_multidata_coil_combine(y_coil, 'thresh', 0.001, 'log2b', 10);
        
        % Save coil-combined data for future use
        dir = pwd; cd('../../data/Brain_06,26,15/IR/');
        save('ims_coil_comb_IR.mat', 'y_im', 'smap'); cd(dir);
    end
end

% For inversion recovery, using complex data is important: is_mag off
% Note that this is inconsistent with SPGR and DESS recons
is_mag = 0;

% Crop images for display
x_crop  = [38:217];  nx = length(x_crop);
% y_crop  = [8:227];   ny = length(y_crop);
y_crop = [10:229];  ny = length(y_crop);
y_im   = y_im(x_crop, y_crop, :);
kap_reg = kap_reg(x_crop, y_crop);

% Trick: normalize data by median of non-background data values
% so that the effective regularization beta is scale-invariant
tmp = abs(y_im);
tmp = median(tmp(tmp > 0.1 * max(tmp(:))));
y_im = y_im / tmp;

% Create tight and loose masks
tight_mask = imfill(squeeze(abs(y_im(:,:,ceil(nTI/2))) >...
    0.03*max(col(abs(y_im(:,:,ceil(nTI/2)))))), 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 10));


%% Dictionary creation (could be precomputed)
% Estimate a constant flip-angle scaling factor
kap_scale = fudge_factor * mean(kap_reg(tight_mask));       % Hopefully approx 1
% kap_scale = mean(kap_reg(roi_mask));                        % Hopefully approx 1

% Compute dictionary signal pathways
% T1 = logspace(1, 4.5, 1000);
T1 = logspace(2, 4, 300);
D = NaN(nTI, length(T1));
K = length(T1);
for t1 = 1:length(T1);
    for i = 1:nTI
        D(i, t1) = IR_fun_v1(1, T1(t1), kap_scale, flip_inv, flip_ex,...
            TR(i), TI(i), TE(i), wf, is_mag);
    end
end


%% Maximum-likelihood Estimation
% Dictionary-based estimation via variable-projection method
weights = ones(nTI, 1);
W = spdiags(weights, 0, nTI, nTI);      % Weighting matrix
y = reshape(permute(y_im, [3 1 2]), [nTI nx*ny]);

tic; 
maxProd = zeros(1, nx*ny);
t1_idx = ones(1, nx*ny);
for k = 1:K
    % Compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;
    
    % If the kth inner product is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    t1_idx(update) = k;
end
time_ml = toc;

% Extract indices for maximum-likelihood maps
T1_ml = reshape(T1(t1_idx), [nx ny]);

% Assume negligible off-resonance effects (short TE)
wf_ml = zeros(nx, ny);

% % M02 initial guess (from dictionary)
% M02_ml = NaN(nx*ny, 1);
% for q = 1:length(t1_idx)
%     M02_ml(q) = (D(:,t1_idx(q))' * y(:,q)) ./ (D(:,t1_idx(q))' * D(:,t1_idx(q)));
% end
% M02_ml = reshape(M02_ml, [nx ny]);


%% Preprocessing and Masking
% Project Images to within range
T1max = 3000;       T1_ml = min(T1_ml, T1max);
T1min = 10;         T1_ml = max(T1_ml, T1min); 

% Set voxels inside loose mask but outside tight mask to mean
T1_ml(~tight_mask & loose_mask) = mean(col(T1_ml));

% Set voxels outside both tight and loose masks to zero
T1_ml(~loose_mask) = 0;
wf_ml(~loose_mask) = 0;

% Median filtering
T1_med = medfilt2(T1_ml);


%% Regularized, Joint M02, T1 Estimation
% Define iteration parameters
n_outer = 100;
niterM = 100;
niter1 = 100;

tolM = 10^-7;
tol1 = 10^-7;
disp = 0;

% Define regularizers, Rm and R1
betaM = nTI * 10^-9;                        % nTI * 2^-10;    
beta1 = nTI * 10^-7;                        % nTI * 2^-9;

% Delta values on the order of one standard deviation
deltaM = 10^-0.5;                           % 10^4;
delta1 = 10^1.5;                            % 10^-1;

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);

% Reshape data inputs to work with 3D implementation
y_3D = permute(y_im, [1 2 4 3]);

% Regularized reconstruction
tic; 
[M02_ml, M02_med, M02_reg, T1_reg, cost] = ...
    mri_irse_m02t1_map(T1_med, kap_reg, y_3D, TI, TR, loose_mask, T1max, T1min,...
    Rm, R1, n_outer, niterM, niter1, tolM, tol1, is_mag, disp);
time_reg = toc;

% Postprocessing for display
% For real data, can use this mask
% reg_mask = imfill(abs(M02_reg) >= 0.02*abs(max(M02_reg(:))), 'holes');
reg_mask = tight_mask;

% % Remove pixels outside reg_mask
M02_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  
M02_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; 
M02_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; 


%% Images, RMSE, and Comparisons
suffix = sprintf('%uIR_brain', nTI);

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
            roi_mask(:,:,r) = roipoly(T1_ml/T1max);
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
        error('Unrecognized input');
end

figure;
im(roi_mask);

% % ROI rectangles (x,y,w,h)
% posWM  = [ctrX(1)-rad(1) ctrY(1)-rad(1) 2*rad(1)+1 2*rad(1)+1];
% posGM  = [ctrX(2)-rad(2) ctrY(2)-rad(2) 2*rad(2)+1 2*rad(2)+1];
% posSF =  [ctrX(3)-rad(3) ctrY(3)-rad(3) 2*rad(1)+1 2*rad(3)+1];

% ROI boundaries
roi_bound.wm = bwboundaries(roi_mask(:,:,1), 'noholes');
roi_bound.gm = bwboundaries(roi_mask(:,:,2), 'noholes');

% M02 Images 
m0smin = 0; m0smax = 16; m0s_rng = [m0smin m0smax];
figure; im('notick', abs(M02_ml), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M02_ml_', suffix, '.eps')), end;
figure; im('notick', abs(M02_med), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M02_med_', suffix, '.eps')), end;
figure; im('notick', abs(M02_reg), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M02_reg_', suffix, '.eps')), end;

% T1 Images
t1min = 500; t1max = 1500; t1_rng = [t1min t1max]; 
figure; im('notick', T1_ml, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_ml_', suffix, '.eps')), end;
figure; im('notick', T1_med, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_med_', suffix, '.eps')), end;
figure; im('notick', T1_reg, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_reg_', suffix, '.eps')), end;

% T1 Images w/ ROIs
figure; im('notick', T1_ml, t1_rng, 'cbar', ' '); hold on;...
    for i = 1:length(roi_bound.wm)
        bound_wm = roi_bound.wm{i};
        plot(bound_wm(:,1), bound_wm(:,2), 'm', 'LineWidth', 1);
    end
    for i = 1:length(roi_bound.gm)
        bound_gm = roi_bound.gm{i};
        plot(bound_gm(:,1), bound_gm(:,2), 'c', 'LineWidth', 1);
    end
    hold off;
if (pr), print('-depsc', strcat('T1_ml_ROI_', suffix, '.eps')); end;

% Cost vs. Iteration
figure; hold on;
scatter(1:2:2*n_outer, cost(2:2:end), 'bo');
scatter(2:2:2*n_outer, cost(3:2:end), 'ro');
plot(0:2*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M02 update', 'T1 update');
%print('-depsc','cost_vs_iter.eps');


%% Sanity check: does the data fit the model?
% Note: plots only make sense if TRs and TEs are fixed across scans

% Extract data, M02, T1 means and standard deviations from ROIs
y_re_mean = NaN(nROI, nTI);     y_re_std = NaN(nROI, nTI);
y_im_mean = NaN(nROI, nTI);     y_im_std = NaN(nROI, nTI);

M02_re_ml_mean = NaN(nROI, 1);  M02_re_ml_std = NaN(nROI, 1);
M02_im_ml_mean = NaN(nROI, 1);  M02_im_ml_std = NaN(nROI, 1);
M02_re_reg_mean = NaN(nROI, 1); M02_re_reg_std = NaN(nROI, 1);
M02_im_reg_mean = NaN(nROI, 1); M02_im_reg_std = NaN(nROI, 1);

T1_ml_mean = NaN(nROI, 1);      T1_ml_std = NaN(nROI, 1);
T1_reg_mean = NaN(nROI, 1);     T1_reg_std = NaN(nROI, 1);

for r = 1:nROI
    for i = 1:nTI
        y_re_mean(r,i)  = mean(masker(real(y_im(:,:,i)), roi_mask(:,:,r)));
        y_im_mean(r,i)  = mean(masker(imag(y_im(:,:,i)), roi_mask(:,:,r)));
        y_re_std(r,i)   = std(masker(real(y_im(:,:,i)), roi_mask(:,:,r)));
        y_im_std(r,i)   = std(masker(imag(y_im(:,:,i)), roi_mask(:,:,r)));
    end
    
    M02_re_ml_mean(r)   = mean(masker(real(M02_ml), roi_mask(:,:,r)));
    M02_im_ml_mean(r)   = mean(masker(imag(M02_ml), roi_mask(:,:,r)));
    M02_re_ml_std(r)    = std(masker(real(M02_ml), roi_mask(:,:,r)));
    M02_im_ml_std(r)    = std(masker(imag(M02_ml), roi_mask(:,:,r)));
    
    M02_re_reg_mean(r)  = mean(masker(real(M02_reg), roi_mask(:,:,r)));
    M02_im_reg_mean(r)  = mean(masker(imag(M02_reg), roi_mask(:,:,r)));
    M02_re_reg_std(r)   = std(masker(real(M02_reg), roi_mask(:,:,r)));
    M02_im_reg_std(r)   = std(masker(imag(M02_reg), roi_mask(:,:,r)));
    
    T1_ml_mean(r)       = mean(masker(T1_ml, roi_mask(:,:,r)));
    T1_ml_std(r)        = std(masker(T1_ml, roi_mask(:,:,r)));
    
    T1_reg_mean(r)      = mean(masker(T1_reg, roi_mask(:,:,r)));
    T1_reg_std(r)       = std(masker(T1_reg, roi_mask(:,:,r)));
    
%     for i = 1:nTI
%         [y_re_mean(r, i), y_re_std(r, i)]...
%             = multiMeans(real(y_im(:,:,i)), [ctrX(r) ctrY(r)], rad(r));
%         [y_im_mean(r, i), y_im_std(r, i)]...
%             = multiMeans(imag(y_im(:,:,i)), [ctrX(r) ctrY(r)], rad(r));
%     end
%     
%     [M02_re_ml_mean(r), M02_re_ml_std(r)] ...
%         = multiMeans(real(M02_ml), [ctrX(r) ctrY(r)], rad(r));
%     [M02_im_ml_mean(r), M02_im_ml_std(r)] ...
%         = multiMeans(imag(M02_ml), [ctrX(r) ctrY(r)], rad(r));
%     [M02_re_reg_mean(r), M02_re_reg_std(r)] ...
%         = multiMeans(real(M02_reg), [ctrX(r) ctrY(r)], rad(r));
%     [M02_im_reg_mean(r), M02_im_reg_std(r)] ...
%         = multiMeans(imag(M02_reg), [ctrX(r) ctrY(r)], rad(r));
%     
%     [T1_ml_mean(r), T1_ml_std(r)] ...
%         = multiMeans(T1_ml, [ctrX(r) ctrY(r)], rad(r));
%     [T1_reg_mean(r), T1_reg_std(r)] ...
%         = multiMeans(T1_reg, [ctrX(r) ctrY(r)], rad(r));
end

% Print Statistics (previous contents overwritten!)
if (pr)
    fid = fopen(strcat('Summary_06,26,15_', suffix), 'w');
    fprintf(fid, strcat('Estimation Statistics for 2 RF SE IR Scans'));
    fprintf(fid, sprintf('\n\tML run time: %0.2f seconds', time_ml));
    fprintf(fid, sprintf('\n\tReg run time: %0.2f seconds', time_reg));
    
    fprintf(fid, '\n\nRe[M02] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M02_re_ml_mean(r),   char(177),  M02_re_ml_std(r),...
            M02_re_reg_mean(r),  char(177),  M02_re_reg_std(r));
    end
    
    fprintf(fid, '\n\nIm[M02] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M02_im_ml_mean(r),   char(177),  M02_im_ml_std(r),...
            M02_im_reg_mean(r),  char(177),  M02_im_reg_std(r));
    end

    fprintf(fid, '\n\nT1 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            T1_ml_mean(r), char(177), T1_ml_std(r),...
            T1_reg_mean(r), char(177), T1_reg_std(r));
    end
    fclose(fid);
end

% Use average values in ROIs to get corresponding ML parameters
% This prevents data noise from being amplified in estimating curves
maxProd_mean = zeros(1, nROI); 
t1_idx_mean = zeros(1, nROI);
y_mean = complex(y_re_mean, y_im_mean);
for k = 1:K
    % Compute kth inner product
    hess_mean = abs(D(:,k)' * W * D(:,k));
    ytild_mean = D(:,k)' * W * y_mean' / sqrt(hess_mean);
    newProd_mean = abs(ytild_mean).^2;
    
    % If the kth inner product is largest, save k
    update_mean = newProd_mean > maxProd_mean;
    maxProd_mean(update_mean) = newProd_mean(update_mean);
    t1_idx_mean(update_mean) = k;
end

% Extract indices for method-of-moments maps
T1_sig_mean = T1(t1_idx_mean);

% M0 initial guess (from dictionary)
M02_sig_mean = NaN(nROI, 1);
for q = 1:length(t1_idx_mean)
    M02_sig_mean(q) = (D(:,t1_idx_mean(q))' * y_mean(q,:).') ./...
        (D(:,t1_idx_mean(q))' * D(:,t1_idx_mean(q)));
end

% Signal model vs. Data
% Note: flip angle compensation using average value, kap_scale
TI_model = 0:1500;
Smean = NaN(nROI, length(TI_model));
for r = 1:nROI
    Smean(r,:) = IR_fun(M02_sig_mean(r), T1_sig_mean(r), kap_scale, flip_inv,...
        flip_ex, mean(TR), TI_model, mean(TE), wf, is_mag);
end

% Real components
figure; hold on;
ColorOrd = get(gca, 'ColorOrder');
for r = 1:nROI
    errorbar(TI, y_re_mean(r,:), y_re_std(r,:), '*', 'Color', ColorOrd(r,:));
    plot(TI_model, real(Smean(r,:)), 'Color', ColorOrd(r,:));
end
title('GE 3T 2D IR SE Re[Data] (d) vs. Re[Model] (m)');
xlabel('Inversion time (TI)'); ylabel('Real Signal (a.u.)');
legend('WM-d', 'WM-m', 'GM-d', 'GM-m', 'CSF-d', 'CSF-m', 'Location', 'NE');
if (pr), print('-depsc', strcat('sig_real_', suffix)); end;

% Imaginary components
figure; hold on;
ColorOrd = get(gca, 'ColorOrder');
for r = 1:nROI
    errorbar(TI, y_im_mean(r,:), y_im_std(r,:), '*', 'Color', ColorOrd(r,:));
    plot(TI_model, imag(Smean(r,:)), 'Color', ColorOrd(r,:));
end
title('GE 3T 2D IR SE Im[Data] (d) vs. Im[Model] (m): Odd ROIs');
xlabel('Inversion time (TI)'); ylabel('Imaginary Signal (a.u.)');
legend('WM-d', 'WM-m', 'GM-d', 'GM-m', 'CSF-d', 'CSF-m', 'Location', 'NE');
if (pr), print('-depsc', strcat('sig_imag_', suffix)); end;
