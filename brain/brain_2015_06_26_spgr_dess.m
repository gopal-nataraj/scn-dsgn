%% Header Script for Regularized M0s/T1/T2 Recon from SPGR/DESS
%       05.15.2015 -- Pulse sequences prepared
%       06.26.2015 -- Brain Data acquired (J. A. Fessler)
%       07.28.2015 -- BS .wav file found to be 8kHz, not 4kHz; b1 corrected
%       08.04.2015 -- Choose reg parameters carefully
%                     Matched dictionary T1/T2 ranges to NIST 06/12 data
%       05.13.2016 -- Added manual roi selection capability
%
% Collected for comparing optimized (2,1), (1,1), and (0,2) scans
% Written by: Gopal Nataraj

%% Raw data extraction
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    curdir = pwd; cd ../../../irt; setup(); cd(curdir);
end
addpath('../../data/Brain_06,26,15/');
addpath('../../model/dess/');
addpath('../../model/spgr/');
addpath('../../map/b1/');
addpath('../../map/sense/');
addpath('../../map/t1-t2/');

% Global imaging parameters
nx = 256; ny = 256; nz = 6; nc = 32;
FOVx = 240; FOVy = 240; FOVz = 30;          % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

TE_global = 4.67;                           % ms
ncycles = 2;                                % Number of spoiling cycles
wf = 0;                                     % Assume zero off-resonance
sl = 2;                                     % Central slice selected
% rs = 0;                                     %       Toggle ROI selection prompts on/off
% pr = 0;                                     % Toggle printing on/off


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
        load('bs_p4000Hz.mat'); ybs_p_coil = permute(flipdim(implus(:,:,sl,:),  1), [2 1 4 3]);
        load('bs_m4000Hz.mat'); ybs_m_coil = permute(flipdim(imminus(:,:,sl,:), 1), [2 1 4 3]);
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


%% Coil-combine the relevant SPGR and DESS data
% Choose a scan profile: (0,2), (1,1), or (2,1)
% n_spgr = 2;                                 % Number of SPGR scans
% n_dess = 1;                                 % Number of DESS scans
profile = sprintf('%u SPGR, %u DESS', n_spgr, n_dess);
clear('ys_im', 'yp_im', 'ym_im');           % If profile changed, reload data

% Set profile-specific imaging parameters
switch profile
    case '0 SPGR, 2 DESS'
        flip_s_deg = []';                   % degrees
        flip_d_deg = [35, 10]';             % degrees
        TRs = []';                          % ms
        TRd = [24.4, 17.5]';                % ms
    case '1 SPGR, 1 DESS'
        flip_s_deg = [15]';                 % degrees
        flip_d_deg = [10]';                 % degrees
        TRs = [13.9]';                      % ms
        TRd = [28.0]';                      % ms
    case '2 SPGR, 1 DESS'
        flip_s_deg = [15, 5]';              % degrees
        flip_d_deg = [30]';                 % degrees
        TRs = [12.2, 12.2]';                % ms
        TRd = [17.5]';                      % ms
    otherwise 
        warning('Unexpected scan profile!');
end
nfs = length(TRs);                          % Number of SPGR scans
nfd = length(TRd);                          % Number of DESS scans
M = nfs + 2*nfd;                            % Number of datasets
flip_s = flip_s_deg * (pi/180);             % radians
flip_d = flip_d_deg * (pi/180);             % radians
TEs = TE_global * ones(nfs, 1);             % ms
TEd = TE_global * ones(nfd, 1);             % ms

if (~exist('yp_im', 'var'))
    try 
        % First try to load SPGR/DESS coil-combined data
        addpath(sprintf('../../data/Brain_06,26,15/%uSPGR,%uDESS', n_spgr, n_dess));
        load(sprintf('ims_coil_comb_%uSPGR%uDESS_slice%u.mat', n_spgr, n_dess, sl));
    catch
        ys_coil = NaN(nx, ny, nfs, nc);
        yp_coil = NaN(nx, ny, nfd, nc);
        ym_coil = NaN(nx, ny, nfd, nc);
        
        % Grab one mediocoronal slice of SPGR data
        dir = pwd; cd('/Volumes/a2/gnataraj/Data/2015,06,26_brain_32ch/SPGR');
        for l = 1:nfs
            file = sprintf('SPGR_TR%2.1f_flip%02d', TRs(l), flip_s_deg(l));
            load(strcat(strrep(file, '.', 'p'), '.mat'));
            ys_coil(:,:,l,:) = ims(:,:,sl,:);                   % [nx ny nc]
        end
        cd(dir);
        
        % Grab one mediocoronal slice of DESS data
        dir = pwd; cd('/Volumes/a2/gnataraj/Data/2015,06,26_brain_32ch/DESS');
        for l = 1:nfd
            file = sprintf('DESS_TR%2.1f_flip%02d_echo1', TRd(l), flip_d_deg(l));
            load(strcat(strrep(file, '.', 'p'), '.mat'));
            yp_coil(:,:,l,:) = implus(:,:,sl,:);                % [nx ny nc]
            
            file = sprintf('DESS_TR%2.1f_flip%02d_echo2', TRd(l), flip_d_deg(l));
            load(strcat(strrep(file, '.', 'p'), '.mat'));
            ym_coil(:,:,l,:) = imminus(:,:,sl,:);               % [nx ny nc]
        end
        cd(dir);
        
        % Take transpose to orient anterior --> posterior
        % Flip image to compensate DESS second echo inverted in readout direction
        ys_coil = permute(flipdim(ys_coil, 1), [2 1 3 4]);
        yp_coil = permute(flipdim(yp_coil, 1), [2 1 3 4]);
        ym_coil = permute(flipdim(ym_coil, 1), [2 1 3 4]);
        
        % Use mri_multidata_coil_combine.m to merge nc-channel coil data
        % Heavy smap regularization required to avoid errors in y_im
        y_coil = cat(3, ys_coil, yp_coil, ym_coil);         % [nx ny M nc]
        [y_im, smap, y_im_ssos, s_ml, coil_cost] = ...
            mri_multidata_coil_combine(y_coil, 'log2b', 20);% [nx ny M]
        ys_im = y_im(:,:,1:nfs);                            % [nx ny nfs]
        yp_im = y_im(:,:,nfs+1:nfs+nfd);                    % [nx ny nfd]
        ym_im = y_im(:,:,nfs+nfd+1:end);                    % [nx ny nfd]
        
        % Save coil-combined data for future use
        dir = pwd; cd(sprintf('../../data/Brain_06,26,15/%uSPGR,%uDESS', n_spgr, n_dess));
        save(sprintf('ims_coil_comb_%uSPGR%uDESS_slice%u.mat', n_spgr, n_dess, sl),...
            'ys_im', 'yp_im', 'ym_im', 'smap'); cd(dir);
    end
end
  
% Even though the data is complex, take the magnitude and turn is_mag on
%   1) Large linear phase observed in main field direction; as this phase
%      is not the same for both echoes, would cause phase errors. 
%   2) Susceptibility effects could cause large off-resonance, distorting phase
ys_im = abs(ys_im);
yp_im = abs(yp_im);
ym_im = abs(ym_im);
is_mag = 1;

% Crop images for display
x_crop  = [38:217];  nx = length(x_crop);
y_crop  = [8:227];   ny = length(y_crop);
ys_im   = ys_im(x_crop, y_crop, :);
yp_im   = yp_im(x_crop, y_crop, :);
ym_im   = ym_im(x_crop, y_crop, :);
kap_reg = kap_reg(x_crop, y_crop);

% Create tight and loose masks
tight_mask = imfill(squeeze(abs(yp_im(:,:,ceil(nfd/2))) >...
    0.05*max(col(abs(yp_im(:,:,ceil(nfd/2)))))), 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 10)); 
        

%% Dictionary creation (could be precomputed)
% Estimate a constant flip-angle scaling factor
kap_scale = fudge_factor * mean(kap_reg(tight_mask));       % Hopefully approx 1 
% kap_scale = mean(kap_reg(roi_mask));                        % Hopefully approx 1

% Compute dictionary signal pathways
T1 = logspace(2, 4, 300);
T2 = logspace(1, 3, 300);
D = NaN(M, length(T1), length(T2));
K = length(T1) * length(T2);
for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        % SPGR dictionary component
        for a = 1:nfs
            D(a, t1, t2) = spgr_fun(1, T1(t1), kap_scale,...
                flip_s(a), TRs(a), TEs(a), wf, 1);
        end

        % DESS dictionary component
        for a = 1:nfd
            [D(nfs+a, t1, t2), D(nfs+nfd+a, t1, t2)] = ...
                dess_fun(1, T1(t1), T2(t2), kap_scale,...
                flip_d(a), TRd(a), TEd(a), wf, 1);
        end
    end
end
D = reshape(D, [M, K]);


%% Maximum-likelihood Estimation
% Dictionary-based estimation via variable-projection method
weights = ones(M, 1);
W = spdiags(weights, 0, M, M);      % Weighting Matrix
y = reshape(permute(cat(3, ys_im, yp_im, ym_im), [3 1 2]), [M nx*ny]);

tic;
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
time_ml = toc;

% Extract indices for maximum-likelihood maps
[t1_idx, t2_idx] = ind2sub([length(T1) length(T2)], idx);
T1_ml = reshape(T1(t1_idx), [nx ny]); 
T2_ml = reshape(T2(t2_idx), [nx ny]);        

% % M0s initial guess
% M0s_ml = NaN(nx*ny, 1);
% for q = 1:length(idx)
%     M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
% end
% M0s_ml = reshape(M0s_ml, [nx ny]);

% Assume negligible off-resonance effects (short TE)
wf_ml = zeros(nx, ny);        
        

%% Preprocessing and Masking
% Project Images to within range
T1max = 5000;       T1_ml = min(T1_ml, T1max);
T1min = 5;          T1_ml = max(T1_ml, T1min);
T2max = 2000;       T2_ml = min(T2_ml, T2max);
T2min = 5;          T2_ml = max(T2_ml, T2min);

% Set voxels inside loose mask but outside tight mask to mean
T1_ml(~tight_mask & loose_mask) = mean(col(T1_ml));   
T2_ml(~tight_mask & loose_mask) = mean(col(T2_ml));   
wf_ml(~tight_mask & loose_mask) = mean(col(wf_ml)); 

% Set voxels outside both tight and loose masks to zero
T1_ml(~loose_mask) = 0;  
T2_ml(~loose_mask) = 0;  
wf_ml(~loose_mask) = 0;  

% Median filtering
T1_med = medfilt2(T1_ml);  
T2_med = medfilt2(T2_ml); 


%% Regularized, Joint M0s,T1,T2,kap Estimation
% Define maximum iteration parameters
n_outer = 100;
niterM = 100;   
niter1 = 100;   
niter2 = 100;   
niterK = 0;                                 % kap not updated

tolM = 10^-7;   
tol1 = 10^-7;  
tol2 = 10^-7;  
tolk = 10^-7;                               % irrelevant
disp = 0; 

% Define regularizers, Rm, R1, R2, and Rk
betaM = M * 10^-9;                          
beta1 = M * 10^-7;                          
beta2 = M * 10^-8;                          
betaK = M * 10^-1;                          % irrelevant

% Delta values on the order of one standard deviation
deltaM = 10^-0.5;
delta1 = 10^1.5;
delta2 = 10^0.5;

Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM, 'type_penal', 'mat');
R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);
Rk = Reg1(loose_mask, 'pot_arg', {'quad'}, 'beta', betaK);

% Reshape data inputs to work with 3D implementation
ys_3D = permute(ys_im, [1 2 4 3]);
yp_3D = permute(yp_im, [1 2 4 3]);
ym_3D = permute(ym_im, [1 2 4 3]);

% % Method One: Sequential Estimation (old code, needs edit)
% % Need a M0s initial guess
% M0s_ml = NaN(nx*ny, 1);
% for q = 1:length(idx)
%     M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
% end
% M0s_ml = reshape(M0s_ml, [nx ny]);
% M0s_med = M0s_ml;
% 
% % M0s, T1 regularized reconstruction from SPGR
% [M0s_reg, T1_reg, ~] = mri_despot1_m0t1map(M0s_ml, E1s_med,...
%     ys_3D, flip_s_3D, loose_mask, Rm, R1, T1max, TR_s, n_outer,...
%     niterM, niter1, tolM, tol1, disp);
% 
% % T2 regularized reconstruction from DESS
% E1d_reg = exp(-TR_d ./ T1_reg);
% [M0s_reg, T1_reg, T2_reg, wf_reg, ~] = ...
%     mri_dess_t2map_nest83_M0s(M0s_reg, E1d_reg, E2d_med, wf_ml, ...
%     flip_d_3D, yp_3D, ym_3D, loose_mask, T2max, T2min,...
%     TR_d, TE, R2, n_outer, niter2, tol2, disp);

% Method Two: Joint Estimation
% Assume we know the ground truth flip angle scaling
tic; 
[M0s_ml, M0s_med, M0s_reg, T1_reg, T2_reg, kap_reg, wf_reg, cost] = ...
    mri_spgr_dess_m0st1t2kap_map(T1_med, T2_med, kap_reg, wf_ml, flip_s,...
    flip_d, ys_3D, yp_3D, ym_3D, loose_mask, T1max, T1min, T2max, T2min,...
    TRs, TRd, TEs, TEd, Rm, R1, R2, Rk, n_outer, niterM, niter1,...
    niter2, niterK, tolM, tol1, tol2, tolk, is_mag, disp);
time_reg = toc; 

% Postprocessing for display
% For real data, can use this mask
% reg_mask = imfill(abs(M0s_reg) >= 0.1*abs(max(M0s_reg(:))), 'holes');
reg_mask = tight_mask;

% Remove pixels outside reg_mask
M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  
M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; 
M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0; 

% %% For comparison, Method-of-Moments T2 Estimation (biased)
% E2sqrd = NaN(nx, ny);
% for xx = 1:nx
%     for yy = 1:ny
%         Y = squeeze(ym_im(xx, yy, :));
%         A = squeeze(yp_im(xx, yy, :));
%         C = [A Y];
%         sig_small = min(svds(C'*C));
%         E2sqrd(xx, yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
%     end
% end
% T2_mom = -2*(mean(TRd)-TE_global) ./ log(E2sqrd);
% T2_mom(~reg_mask) = 0;


%% Images, RMSE, and Comparisons
suffix = sprintf('%uSPGR%uDESS_brain', nfs, nfd);

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

% M0s Images 
m0smin = 0; m0smax = 16; m0s_rng = [m0smin m0smax];
figure; im('notick', abs(M0s_ml), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_ml_', suffix, '.eps')), end;
figure; im('notick', abs(M0s_med), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_med_', suffix, '.eps')), end;
figure; im('notick', abs(M0s_reg), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_reg_', suffix, '.eps')), end;

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

% T2 Images
t2min = 20; t2max = 120; t2_rng = [t2min t2max]; 
% figure; im('notick', T2_mom, t2_rng, 'cbar', ' ');
% if (pr), print('-deps', strcat('T2_mom_', suffix, '.eps')), end;
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

% kap Images
kapmin = 0; kapmax = 2;
figure; im('notick', kap_reg, [kapmin kapmax], 'cbar', ' ');
if (pr), print('-deps', strcat('kap_reg_', suffix)), end;

% Cost vs. Iteration
figure; hold on;
scatter(1:4:4*n_outer, cost(2:4:end), 'bo');
scatter(2:4:4*n_outer, cost(3:4:end), 'ro');
scatter(3:4:4*n_outer, cost(4:4:end), 'go');
scatter(4:4:4*n_outer, cost(5:4:end), 'mo');
plot(0:4*n_outer, cost, 'k'); hold off;
title('Cost vs. iteration');
legend('M0s update', 'T1 update', 'T2 update', 'kap update');
% print('-depsc', strcat('cost_', suffix, '.eps'));


%% Sanity check: does the data fit the model?
% Note: plots only make sense if TRs are fixed across scans

% Extract data means and standard deviations from ROI
ys_mean = NaN(nROI, nfs); ys_std = NaN(nROI, nfs);
yp_mean = NaN(nROI, nfd); yp_std = NaN(nROI, nfd);
ym_mean = NaN(nROI, nfd); ym_std = NaN(nROI, nfd);

M0s_re_ml_mean = NaN(nROI, 1);  M0s_re_ml_std = NaN(nROI, 1);
M0s_im_ml_mean = NaN(nROI, 1);  M0s_im_ml_std = NaN(nROI, 1);
M0s_re_reg_mean = NaN(nROI, 1); M0s_re_reg_std = NaN(nROI, 1);
M0s_im_reg_mean = NaN(nROI, 1); M0s_im_reg_std = NaN(nROI, 1);

T1_ml_mean = NaN(nROI, 1);      T1_ml_std = NaN(nROI, 1);
T1_reg_mean = NaN(nROI, 1);     T1_reg_std = NaN(nROI, 1);

T2_ml_mean = NaN(nROI, 1);      T2_ml_std = NaN(nROI, 1);
T2_reg_mean = NaN(nROI, 1);     T2_reg_std = NaN(nROI, 1);
% T2_mom_mean = NaN(nROI, 1);     T2_mom_std = NaN(nROI, 1);

kap_reg_mean = NaN(nROI, 1);    kap_reg_std = NaN(nROI, 1);

for r = 1:nROI
    for l = 1:nfs
        ys_mean(r,l)     = mean(masker(ys_im(:,:,l), roi_mask(:,:,r)));
        ys_std(r,l)      = std(masker(ys_im(:,:,l), roi_mask(:,:,r)));
    end
    
    for l = 1:nfd
        yp_mean(r,l)     = mean(masker(yp_im(:,:,l), roi_mask(:,:,r)));
        yp_std(r,l)      = std(masker(yp_im(:,:,l), roi_mask(:,:,r)));
        
        ym_mean(r,l)     = mean(masker(ym_im(:,:,l), roi_mask(:,:,r)));
        ym_std(r,l)      = std(masker(ym_im(:,:,l), roi_mask(:,:,r)));
    end
    
    T1_ml_mean(r)       = mean(masker(T1_ml, roi_mask(:,:,r)));
    T1_ml_std(r)        = std(masker(T1_ml, roi_mask(:,:,r)));
    T1_reg_mean(r)      = mean(masker(T1_reg, roi_mask(:,:,r)));
    T1_reg_std(r)       = std(masker(T1_reg, roi_mask(:,:,r)));
    
    T2_ml_mean(r)       = mean(masker(T2_ml, roi_mask(:,:,r)));
    T2_ml_std(r)        = std(masker(T2_ml, roi_mask(:,:,r)));
    T2_reg_mean(r)      = mean(masker(T2_reg, roi_mask(:,:,r)));
    T2_reg_std(r)       = std(masker(T2_reg, roi_mask(:,:,r)));
    
    kap_reg_mean(r)     = mean(masker(kap_reg, roi_mask(:,:,r)));
    kap_reg_std(r)      = std(masker(kap_reg, roi_mask(:,:,r)));
end

% Print Statistics (previous contents overwritten!)
if (pr)
    fid = fopen(strcat('Summary_06,26,15_', suffix), 'w');
    fprintf(fid, sprintf('Estimation Statistics for (%u,%u) Profile)',...
        n_spgr, n_dess));
    fprintf(fid, sprintf('\n\tML run time: %0.2f seconds', time_ml));
    fprintf(fid, sprintf('\n\tReg run time: %0.2f seconds', time_reg));
    
    fprintf(fid, '\n\nRe[M0s] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M0s_re_ml_mean(r),   char(177),  M0s_re_ml_std(r),...
            M0s_re_reg_mean(r),  char(177),  M0s_re_reg_std(r));
    end
    
    fprintf(fid, '\n\nIm[M0s] Summary Statistics\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n', roi_labels{r},...
            M0s_im_ml_mean(r),   char(177),  M0s_im_ml_std(r),...
            M0s_im_reg_mean(r),  char(177),  M0s_im_reg_std(r));
    end

    fprintf(fid, '\n\nT1 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\n', roi_labels{r},...
            T1_ml_mean(r), char(177), T1_ml_std(r),...
            T1_reg_mean(r), char(177), T1_reg_std(r));
    end

    fprintf(fid, '\n\nT2 Summary Statistics:\n');
    fprintf(fid, '\tML\t\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.2f\t%c%0.2f\t%0.2f\t%c%0.2f\n', roi_labels{r},...
            T2_ml_mean(r), char(177), T2_ml_std(r),...
            T2_reg_mean(r), char(177), T2_reg_std(r));
    end
        
    fprintf(fid, '\n\nKappa Summary Statistics:\n');
    fprintf(fid, '\tReg\n');
    for r = 1:nROI
        fprintf(fid, '%s:\t%0.2f\t%c%0.2f\n', roi_labels{r},...
            kap_reg_mean(r), char(177), kap_reg_std(r));
    end
    fclose(fid);
end

% Use average values in ROIs to get corresponding ML parameters
% This prevents data noise from being amplified in estimating curves
y_mean = [ys_mean yp_mean ym_mean]'; 
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

% Extract indices for ML maps (kap_sig_mean constant)
[t1_idx_mean, t2_idx_mean] ...
    = ind2sub([length(T1) length(T2)], idx_mean);
T1_sig_mean = T1(t1_idx_mean);
T2_sig_mean = T2(t2_idx_mean);
kap_sig_mean = kap_reg_mean';
% kap_sig_mean = kap_scale * ones(size(idx_mean));

% M0s initial guess (from dictionary)
M0s_sig_mean = NaN(nROI, 1);
for q = 1:length(idx_mean)
    M0s_sig_mean(q) = (D(:,idx_mean(q))' * y_mean(:,q)) ./...
        (D(:,idx_mean(q))' * D(:,idx_mean(q)));
end

% Signal model vs. Data
for r = 1:nROI
    flip_model = [0:90] * (pi/180) * kap_sig_mean(r);
    flips_data = flip_s * kap_sig_mean(r);
    flipd_data = flip_d * kap_sig_mean(r);
    
    Ss_mean = spgr_fun(M0s_sig_mean(r), T1_sig_mean(r), 1,...
        flip_model, mean(TRs), TE_global, wf, 1);
    [Sp_mean, Sm_mean] = dess_fun(M0s_sig_mean(r), T1_sig_mean(r),...
        T2_sig_mean(r), 1, flip_model, mean(TRd), TE_global, wf, 1);
    
    figure; hold on;
    if (nfs > 0)
        errorbar(flips_data * (180/pi), abs(ys_mean(r,:)), ys_std(r,:), 'g*');
        plot(flip_model * (180/pi), abs(Ss_mean), 'y');
    end
    if (nfd > 0)
        errorbar(flipd_data * (180/pi), abs(yp_mean(r,:)), yp_std(r,:), 'b*');
        plot(flip_model * (180/pi), abs(Sp_mean), 'c');

        errorbar(flipd_data * (180/pi), abs(ym_mean(r,:)), ym_std(r,:), 'r*');
        plot(flip_model * (180/pi), abs(Sm_mean), 'm');
    end
    hold off;

    if (nfs==0)
        legend('SSFP-FID Data', 'SSFP-FID ML Model', 'SSFP-Echo Data', 'SSFP-Echo ML Model');
    elseif (nfd==0)
        legend('SPGR Data', 'SPGR ML Model');
    else
        legend('SPGR Data', 'SPGR ML Model', 'SSFP-FID Data',...
            'SSFP-FID ML Model', 'SSFP-Echo Data', 'SSFP-Echo ML Model');
    end
    tit = sprintf('Brain %s: DESS Model vs. Data, (T1,T2,kap) = (%0.1fms,%0.2fms,%0.2f)',...
        roi_labels{r}, T1_reg_mean(r), T2_reg_mean(r), kap_reg_mean(r)); title(tit); 
    xlabel('(Compensated) flip angle (deg)'); 
    ylabel('Magnitude Signal (a.u.)');
    if (pr), print('-depsc', strcat(sprintf('sig_%s_', roi_labels{r}), suffix, '.eps')), end;
end             