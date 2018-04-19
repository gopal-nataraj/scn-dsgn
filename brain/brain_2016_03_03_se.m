% script brain_2016_03_03_se.m
% compares t2 estimates as echo times are varied
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2016-02-24      protocol prepared
%   v0.2    2016-03-03      data acquired
%   v1.1    2016-03-04      data processed; original

%% raw data extraction
% header files and IRT Setup
if (~exist('irtdir', 'var'))
    cd ../irt; setup(); cd ../matlab;
end
addpath('Brain_03,03,16/');
addpath('Calculations/SE/');
addpath('Calculations/SENSE/');
addpath('Helper/');

% imaging parameters
nx = 256; ny = 256; nz = 1; nc = 32;        
FOVx = 240; FOVy = 240; FOVz = 5;           % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

TE = [10 30 150]';                          % ms
nTE = length(TE);                           
TR = 3000 * ones(nTE, 1);                   % ms
flip_ex = pi/2;                             % rad   Nominal excitation flip
wf = 0;                                     % Hz    Assume zero off-resonance
kap = 1;                                    %       Flip angle scaling
pr = 0;                                     % Toggle printing on/off

%% coil-combine the SE coil data
if (~exist('y_im', 'var'))
    try
        % first try to load SE coil-combined data
        addpath('Brain_03,03,16/');
        load('ims_coil_comb_SE.mat');
    catch
        % grab the SE coil data
        path = '/Volumes/General Storage/Documents/school/michigan/research/data/2016,03,03_brain_32ch/';
        dir = pwd; cd(path);
        load('im_se_coil.mat'); cd(dir);
        
        % use mri_multidata_coil_combine.m to merge nc-channel coil data
        y_coil = permute(im_se, [2 1 4 3]);                 % [nx ny nTE nc]
        [y_im, smap] = mri_multidata_coil_combine(y_coil);  % [nx ny nTE]
        
        % save coil-combined data for future use
        dir = pwd; cd('Brain_03,03,16/');
        save('ims_coil_comb_SE.mat', 'y_im', 'smap'); cd(dir);
    end
end

% even though the data is complex, take the magnitude and turn is_mag on
% this is to remain consistent with the SPGR and DESS recons
y_im = abs(y_im);
is_mag = 1;

% trick: normalize data by median of non-background data values
%   so that the effective regularization beta is scale-invariant
tmp = abs(y_im);
tmp = median(tmp(tmp > 0.1 * max(tmp(:))));
y_im = y_im / tmp;

% crop images for display
x_crop  = [38:217];  nx = length(x_crop);
y_crop  = [31:256];  ny = length(y_crop);
y_im    = y_im(x_crop, y_crop, :);

% create tight and loose masks
tight_mask = imfill(squeeze(abs(y_im(:,:,ceil(nTE/2))) >...
    0.01*max(col(abs(y_im(:,:,ceil(nTE/2)))))), 'holes');
tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
loose_mask = imdilate(tight_mask, strel('disk', 10));


%% dictionary creation (could be precomputed)
% T2 = logspace(0, 3.5, 1000);
T2 = logspace(1, 3, 300);
D = NaN(nTE, length(T2));
K = length(T2);

for t2 = 1:length(T2)
    for i = 1:nTE
        D(i, t2) = SE_fun(1, T2(t2), kap, flip_ex, TE(i), is_mag);
    end
end

%% maximum-likelihood estimates from different TE combinations
n_comb = 4;
T2_ml_part = NaN(nx, ny, n_comb);
time_ml_part = NaN(n_comb, 1);

for n = 1:n_comb
    switch n
        case 1  % te = [10,30] ms
            idx = [1; 2]; 
            fprintf('T2 ML from TE = (%u,%u)ms...', TE(idx));
        case 2  % te = [10,150] ms
            idx = [1; 3];
            fprintf('T2 ML from TE = (%u,%u)ms...', TE(idx));
        case 3  % te = [30, 150] ms
            idx = [2; 3];
            fprintf('T2 ML from TE = (%u,%u)ms...', TE(idx));
        case 4  % te = [10,30,150] ms
            idx = [1; 2; 3];
            fprintf('T2 ML from TE = (%u,%u,%u)ms...', TE(idx));
    end
    
    % extract relevant data and dictionary
    y_im_part = y_im(:,:,idx);
    D_part = D(idx, :);
    nTE_part = length(idx);
    
    % dictionary-based estimation via variable projection method
    weights = ones(nTE_part, 1);
    W_part = spdiags(weights, 0, nTE_part, nTE_part);
    y_part = reshape(permute(y_im_part, [3 1 2]), [nTE_part nx*ny]);
    
    tic;
    maxProd = zeros(1, nx*ny);
    t2_idx = ones(1, nx*ny);
    for k = 1:K
        % Compute kth inner product
        hess = abs(D_part(:,k)' * W_part * D_part(:,k));
        ytild = D_part(:,k)' * W_part * y_part / sqrt(hess);
        newProd = abs(ytild).^2;

        % If the kth inner product is largest, save k
        update = newProd > maxProd;
        maxProd(update) = newProd(update);
        t2_idx(update) = k;
    end
    time_ml_part(n) = toc;
    fprintf('done in %3.0fms.\n', time_ml_part(n)*1000);
    
    % extract indices for maximum-likelihood maps
    T2_ml_part(:,:,n) = reshape(T2(t2_idx), [nx ny]);
end

% assume negligible off-resonance effects and flip angle variation
wf_ml = zeros(nx, ny);
kap_ml = ones(nx, ny);


%% postprocessing and masking
% project images to within range
T2max = 2000;       T2_ml_part = min(T2_ml_part, T2max);
T2min = 1;          T2_ml_part = max(T2_ml_part, T2min);

% Set voxels inside loose mask but outside tight mask to mean
tight_mask = repmat(tight_mask, [1 1 n_comb]);
loose_mask = repmat(loose_mask, [1 1 n_comb]);
T2_ml_part(~tight_mask & loose_mask) = mean(col(T2_ml_part));  

% Set voxels outside both tight and loose masks to zero
T2_ml_part(~loose_mask) = 0;  


%% images, rmse, and comparisons
% t2 ranges
t2min = 20;
t2max = 120;
t2_rng = [t2min t2max];

% define ROIs (WM, GM)
ctrX = [113 79];
ctrY = [60  48];
rad  = [4   2];
nROI = length(rad);
roi_labels = {'WM', 'GM'};

T2_ml_mean  = NaN(nROI, n_comb);
T2_ml_std   = NaN(nROI, n_comb);

% summary statistics
if (pr)
    fid = fopen('2016-03-03,summary', 'w');
    fprintf(fid, 'SE scan ROI statistics, varying echo times');
    fprintf(fid, '\n\nT2 Summary Statistics:\n');
    fprintf(fid, '\t\t\tWM\t\tGM\n');
end

for n = 1:n_comb
    switch n
        case 1
            suffix = '10,30ms';
        case 2
            suffix = '10,150ms';
        case 3
            suffix = '30,150ms';
        case 4
            suffix = '10,30,150ms';
    end
    
    % t2 images
    figure; im('notick', T2_ml_part(:, :, n), t2_rng, 'cbar', ' ');
    if (pr), print('-deps', strcat('T2_ml_SE_TE_', suffix, '.eps')); end;
    
    % roi mean and std. dev.
    for r = 1:nROI
        [T2_ml_mean(r, n), T2_ml_std(r, n)] ...
            = multiMeans(T2_ml_part(:,:,n), [ctrX(r) ctrY(r)], rad(r));
    end
    
    % summary statistics, cont'd.
    if (pr)
        str = '%s:\t%0.3f\t%c%0.3f\t%0.3f\t%c%0.3f\n';
        fprintf(fid, str, suffix,...
            T2_ml_mean(1, n), char(177), T2_ml_std(1, n),...
            T2_ml_mean(2, n), char(177), T2_ml_std(2, n));
    end
end

if (pr)
    fclose(fid);
end

