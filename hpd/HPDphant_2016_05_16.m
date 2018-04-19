% script HPDphant_2016_05_16.m
% compare t1/t2 estimates from optimized pairs of (2,1), (1,1), (0,2)
%   (spgr, dess) scans, from several independent scan repetitions
% warning: some dess coil data is clipped, and weighted zero in coil combination
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   2016-05-16      hpd phantom data acquired
%   2016-05-17      original recon script
%   2016-05-18      masking varied across profile, but now fixed across reps
%   2016-05-19      profile comparison and performance plots vs. published nist numbers

%% initial setup
% header files and irt setup
if (~exist('irtdir', 'var'))
    cd ../irt; setup(); cd ../matlab;
end
addpath('HPD_Phantom_05,16,16/');
addpath('Calculations/BS');
addpath('Calculations/SPGR/');
addpath('Calculations/DESS/');
addpath('Calculations/SENSE/');
addpath('Calculations/Contrib/');
addpath('Helper/');

% global imaging parameters
nx = 256; ny = 256; nz = 8;
FOVx = 240; FOVy = 240; FOVz = 40;          % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

% coil.wghts = [1 0 1 1 1 1 0 1]';            % coil weights (0-wght coils may be clipped)
coil.wghts = [0 0 0 1 1 0 0 0]';            % coil weights (0-wght coils may be clipped)
coil.iter = 10;                             % coil combination outer iteration count
coil.log2b = -5;                            % log2(coil reg parameter)
coil.arg = {...
    'coilweights', coil.wghts,...
    'nouter', coil.iter,...
    'log2b', coil.log2b};

nc = 8;                                     % number of receive coils
nr = 11;                                    % number of repetitions
nprof = 3;                                  % number of profiles to test
ncycles = 2;                                % number of dess spoiling cycles

TE_global = 4.67;                           % ms
wf = 0;                                     % off-resonance

sl = 3;                                     % 5th slice selected ([3:6] only extracted)
shift = [0 0];                              % [x y] pixels to circshift
vi = 0;                                     % toggle viewing reg images on/off
pr = 0;                                     % toggle printing on/off 
sv = 0;                                     % toggle data saving on/off


%% b1 mapping from bloch-siegert data
if (~exist('kap_reg', 'var'))
    try
        addpath('HPD_Phantom_05,16,16/BS/');
        load('B1_reg');
    catch
        % compute bloch-siegert constant
        [~, mag_p, ph_p] = readwav('bsp.wav');
        [~, mag_m, ph_m] = readwav('bsm.wav');
        w_rf = 4000;                            % Hz        (meas. ??? kHz off-res.)
        t = 4e-6 * (0:length(mag_p)-1)';        % s         (4 us sampling)
        gam = 4.258e3;                          % Hz/Gauss
        w_bs = (gam*mag_p).^2/(2*w_rf);         % Hz
        phi_bs = trapz(t, (2*pi*w_bs));         % rad
        K_bs = phi_bs / max(mag_p.^2);          % rad/Gauss^2
        
        % load (3d) bloch-siegert data from all repetitions
        % for now, only use single slice of interest
        ybs_p_coil = NaN(nx, ny, nc, nr);
        ybs_m_coil = NaN(nx, ny, nc, nr);
        dir = pwd; 
        cd('/Volumes/a2/gnataraj/Data/2016,05,16_hpd_8ch/');
        for ir = 1:nr
            load(sprintf('im,hpd-phant,bs-spgr-dess,pfreq,rep%u.mat', ir));
            ybs_p_coil(:,:,:,ir) = ims{1}.echo1(:,:,sl,:);
            ybs_m_coil(:,:,:,ir) = ims{2}.echo1(:,:,sl,:);
        end
        clear('ims');
        cd(dir)
        
        % circshift to center object in fov
        ybs_p_coil = circshift(ybs_p_coil, [shift 0 0]);
        ybs_m_coil = circshift(ybs_m_coil, [shift 0 0]);
        
        % estimate separate b1 maps for each repetition
        b1_init = NaN(nx, ny, nr);
        b1_reg = NaN(nx, ny, nr);
        for ir = 1:nr
            % regularized b1 map estimation and flip angle correction
            [b1_init(:,:,ir), b1_reg(:,:,ir), ~] = ...
                mri_bsb1map(ybs_p_coil(:,:,:,ir), ybs_m_coil(:,:,:,ir),...
                'Kbs', K_bs, 'log2b_b1', 2, 'iter_b1', 1000, 'coilOpt', coil.arg);
        end
        b1_peak = 0.050;                        % Gauss, for 90deg flip
        kap_reg = b1_reg / b1_peak;             % flip angle scale factor
        
        figure; im(b1_init, [0 b1_peak], 'cbar');
        figure, im(b1_reg,  [0 b1_peak], 'cbar');
        figure; im(kap_reg, [0 1], 'cbar');
        
        % Save regularized B1 maps for future use
        dir = pwd; cd('HPD_Phantom_05,16,16/BS/');
        save('B1_reg', 'b1_peak', 'b1_init', 'b1_reg', 'kap_reg'); cd(dir);
    end
    
    % DEBUG: try fudging the scale factor (default = 1)
    fudge_factor = 1.016/0.913;
    kap_reg = fudge_factor * kap_reg;
end


%% global constant instantiation
% define rois
ctrX = [ 74  84 111 144 170 181 170 143 110  84 106 106 148 148];
ctrY = [130 161 180 179 160 128  97  77  78  99 108 150 150 107];
rad  = [  5   5   5   5   5   5   5   5   5   5   5   5   5   5];
nROI = length(rad);

% roi_masks
roi_mask = false(nx, ny, nROI);
[xx, yy] = ndgrid(1:nx, 1:ny);
for r = 1:nROI
    tmp = (xx-ctrX(r)).^2 + (yy-ctrY(r)).^2 <= rad(r)^2;
    roi_mask(:,:,r) = tmp;
end
roi_mask_all = sum(roi_mask, 3) ~= 0;

% ml variable instantiation
time_ml = NaN(nr, nprof);
M0s_ml = NaN(nx, ny, nr, nprof);
T1_ml = NaN(nx, ny, nr, nprof);
T2_ml = NaN(nx, ny, nr, nprof);
wf_ml = zeros(nx, ny, nr, nprof);

% ml postprocessing limits
T1max = 5000;
T1min = 5;
T2max = 2000;
T2min = 5;

% median-filtered variable instantiation
M0s_med = NaN(nx, ny, nr, nprof);
T1_med = NaN(nx, ny, nr, nprof);
T2_med = NaN(nx, ny, nr, nprof);

% max iterations
n_outer = 100;
niterM = 100;   
niter1 = 100;   
niter2 = 100;   
niterK = 0;                                 % kap not updated

% stopping tolerance
tolM = 10^-7;   
tol1 = 10^-7;  
tol2 = 10^-7;  
tolk = 10^-7;                               % irrelevant

% regularizer curvatures
deltaM = 10^-0.5;
delta1 = 10^1.5;
delta2 = 10^0.5;

% regularized variable instantiation
time_reg = NaN(nr, nprof);
M0s_reg = NaN(nx, ny, nr, nprof);
T1_reg = NaN(nx, ny, nr, nprof);
T2_reg = NaN(nx, ny, nr, nprof);
wf_reg = NaN(nx, ny, nr, nprof);

try
    % first try to load maps
    fprintf('\nTrying to load parameter maps...');
    load(sprintf('maps_slice%u.mat', sl+2));
    fprintf('success!\n');
catch
    %% if unavailable, make maps for: (0,2), (1,1), or (2,1)
    fprintf('failed.\n');
    for ip = 1:nprof
        switch ip
            case 1   
                n_spgr = 0;
                n_dess = 2;
                flip_s_deg = []';                   % degrees
                flip_d_deg = [35, 10]';             % degrees
                TRs = []';                          % ms
                TRd = [24.4, 17.5]';                % ms
                mask.ir = 8;                        % mask design rep
                mask.thresh = 0.04;                 % mask design threshold
            case 2
                n_spgr = 1;
                n_dess = 1;
                flip_s_deg = [15]';                 % degrees
                flip_d_deg = [10]';                 % degrees
                TRs = [13.9]';                      % ms
                TRd = [28.0]';                      % ms
                mask.ir = 3;                        % mask design rep
                mask.thresh = 0.05;                 % mask design threshold
            case 3
                n_spgr = 2;
                n_dess = 1;
                flip_s_deg = [15, 5]';              % degrees
                flip_d_deg = [30]';                 % degrees
                TRs = [12.2, 12.2]';                % ms
                TRd = [17.5]';                      % ms
                mask.ir = 6;                        % mask design rep
                mask.thresh = 0.05;                 % mask design threshold
            otherwise 
                warning('Unexpected scan profile!');   
        end
        profile = sprintf('%u SPGR, %u DESS', n_spgr, n_dess);
        nfs = length(TRs);                          % number of SPGR scans
        nfd = length(TRd);                          % number of DESS scans
        M = nfs + 2*nfd;                            % number of datasets
        flip_s = flip_s_deg * (pi/180);             % radians
        flip_d = flip_d_deg * (pi/180);             % radians
        TEs = TE_global * ones(nfs, 1);             % ms
        TEd = TE_global * ones(nfd, 1);             % ms       
        
        try
            % first try to load spgr/dess coil-combined data
            fprintf('\nTrying to load coil-combined images for (%u,%u) profile...',...
                n_spgr, n_dess);
            addpath(sprintf('HPD_Phantom_05,16,16/%uSPGR,%uDESS', n_spgr, n_dess));
            load(sprintf('ims_coil_comb_%uSPGR%uDESS_slice%u.mat', n_spgr, n_dess, sl+2));
            fprintf('success!\n');
        catch
            fprintf('failed.\n');
            ys_coil = NaN(nx, ny, nfs, nc, nr);
            yp_coil = NaN(nx, ny, nfd, nc, nr);
            ym_coil = NaN(nx, ny, nfd, nc, nr);
            
            % if not found extract data from all repetitions
            dir = pwd;
            cd('/Volumes/a2/gnataraj/Data/2016,05,16_hpd_8ch/');
            for ir = 1:nr
                fprintf('Loading coil data for (%u,%u) profile, repetition %u.\n',...
                    n_spgr, n_dess, ir);
                load(sprintf('im,hpd-phant,bs-spgr-dess,pfreq,rep%u.mat', ir));
                switch ip
                    case 1
                        yp_coil(:,:,1,:,ir) = ims{8}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{8}.echo2(:,:,sl,:);
                        yp_coil(:,:,2,:,ir) = ims{9}.echo1(:,:,sl,:);
                        ym_coil(:,:,2,:,ir) = ims{9}.echo2(:,:,sl,:);
                    case 2
                        ys_coil(:,:,1,:,ir) = ims{5}.echo1(:,:,sl,:);
                        yp_coil(:,:,1,:,ir) = ims{7}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{7}.echo2(:,:,sl,:);
                    case 3
                        ys_coil(:,:,1,:,ir) = ims{3}.echo1(:,:,sl,:);
                        ys_coil(:,:,2,:,ir) = ims{4}.echo1(:,:,sl,:);
                        yp_coil(:,:,1,:,ir) = ims{6}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{6}.echo2(:,:,sl,:);
                    otherwise
                        warning('Unexpected scan profile!');
                end
            end
            clear('ims');
            cd(dir);
            
            % make sure coil arrays are fully populated
            if sum(isnan(col(cat(3, ys_coil, yp_coil, ym_coil)))) ~= 0
                error('Incomplete loading of data!');
            end
            
            % circshift to center object in fov
            ys_coil = circshift(ys_coil, [shift 0 0 0]);
            yp_coil = circshift(yp_coil, [shift 0 0 0]);
            ym_coil = circshift(ym_coil, [shift 0 0 0]);
            
            % coil-combine each repetition separately
            y_coil = cat(3, ys_coil, yp_coil, ym_coil);     % [nx ny M nc nr]
            y_im = NaN(nx, ny, M, nr);                      % [nx ny M nr]
            smap = NaN(nx, ny, nc, nr);                     % [nx ny nc nr]
            for ir = 1:nr
                [y_im(:,:,:,ir), smap(:,:,:,ir)] = ...      % [nx ny M], [nx ny nc]
                    mri_multidata_coil_combine(y_coil(:,:,:,:,ir), coil.arg{:});
            end

            % separate out the coil-combined data types
            ys_im = y_im(:,:,1:nfs,:);                      % [nx ny nfs nr]
            yp_im = y_im(:,:,nfs+1:nfs+nfd,:);              % [nx ny nfd nr]
            ym_im = y_im(:,:,nfs+nfd+1:end,:);              % [nx ny nfd nr]
            
            % save coil-combined data for future use
            dir = pwd;
            cd(sprintf('HPD_Phantom_05,16,16/%uSPGR,%uDESS', n_spgr, n_dess));
            save(sprintf('ims_coil_comb_%uSPGR%uDESS_slice%u.mat', n_spgr, n_dess, sl+2),...
                'ys_im', 'yp_im', 'ym_im', 'smap'); 
            cd(dir);
        end
        
        % Even though the data is complex, take the magnitude and turn is_mag on
        %   1) Large linear phase observed in main field direction; as this phase
        %      is not the same for both echoes, would cause phase errors. 
        %   2) Susceptibility effects could cause large off-resonance, distorting phase
        ys_im = abs(ys_im);
        yp_im = abs(yp_im);
        ym_im = abs(ym_im);
        is_mag = 1;
        
        
        %% dictionary creation
        % estimate a constant flip-angle scaling factor
        % to ensure dictionary constant, use all repetitions here
        kap_scale = mean(kap_reg(repmat(roi_mask_all, [1 1 nr])));

        % compute dictionary signal pathways
        T1 = logspace(2, 4, 300);
        T2 = logspace(1, 3, 300);
        D = NaN(M, length(T1), length(T2));
        K = length(T1) * length(T2);
        for t1 = 1:length(T1)
            for t2 = 1:length(T2)
                % SPGR dictionary component
                for a = 1:nfs
                    D(a, t1, t2) = spgr_fun(1, T1(t1), kap_scale,...
                        flip_s(a), TRs(a), TEs(a), wf, is_mag);
                end

                % DESS dictionary component
                for a = 1:nfd
                    [D(nfs+a, t1, t2), D(nfs+nfd+a, t1, t2)] = ...
                        dess_fun(1, T1(t1), T2(t2), kap_scale,...
                        flip_d(a), TRd(a), TEd(a), wf, is_mag);
                end
            end
        end
        D = reshape(D, [M, K]);
        
        
        %% t1,t2 mapping by repetition
        % weighting matrix
        weights = ones(M, 1);
        W = spdiags(weights, 0, M, M);
        
        % fix tight/loose masks for all repetitions
        tight_mask = imfill(squeeze(abs(yp_im(:,:,ceil(nfd/2),mask.ir)) >...
            mask.thresh * max(col(abs(yp_im(:,:,ceil(nfd/2),mask.ir))))), 'holes');
        tight_mask = imdilate(~imdilate(~tight_mask, strel('disk', 5)), strel('disk', 5));
        loose_mask = imdilate(tight_mask, strel('disk', 10));
        figure; im(tight_mask);
        
        % regularization parameters
        betaM = M * 10^-9;                          
        beta1 = M * 10^-7;                          
        beta2 = M * 10^-8;                          
        betaK = M * 10^-1;                          % irrelevant

        % regularizers
        Rm = Reg1(loose_mask, 'pot_arg', {'hyper3', deltaM}, 'beta', betaM,...
            'type_penal', 'mat');
        R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
        R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);
        Rk = Reg1(loose_mask, 'pot_arg', {'quad'}, 'beta', betaK);
        
        for ir = 1:nr
            fprintf('\nPerforming ML est for repetition %u, profile %u.\n', ir, ip);
            
            %% ml estimation
            y = reshape(permute(cat(3, ys_im(:,:,:,ir), yp_im(:,:,:,ir),...
                ym_im(:,:,:,ir)), [3 1 2]), [M nx*ny]);

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
            time_ml(ir, ip) = toc;

            % extract indices for ml maps
            [t1_idx, t2_idx] = ind2sub([length(T1) length(T2)], idx);
            tmp1 = reshape(T1(t1_idx), [nx ny]);    % placeholder for cur rep
            tmp2 = reshape(T2(t2_idx), [nx ny]);    % placeholder for cur rep

            % project ml maps to within range
            tmp1 = min(tmp1, T1max);
            tmp1 = max(tmp1, T1min);
            tmp2 = min(tmp2, T2max);
            tmp2 = max(tmp2, T2min);

            % set voxels inside loose mask but outside tight mask to mean
            tmp1(~tight_mask & loose_mask) = mean(col(tmp1));
            tmp2(~tight_mask & loose_mask) = mean(col(tmp2));

            % set voxels outside both tight and loose masks to zero
            tmp1(~loose_mask) = 0;
            tmp2(~loose_mask) = 0;

            % store ml estimates
            T1_ml(:,:,ir,ip) = tmp1;
            T2_ml(:,:,ir,ip) = tmp2;

            %% regularized estimation
            fprintf('Performing reg est for repetition %u, profile %u.\n', ir, ip);
            
            % median-filtered initialization
            T1_med(:,:,ir,ip) = medfilt2(tmp1);
            T2_med(:,:,ir,ip) = medfilt2(tmp2);

            % reshape data inputs to work with 3d implementation
            ys_3d = permute(ys_im(:,:,:,ir), [1 2 4 3]);
            yp_3d = permute(yp_im(:,:,:,ir), [1 2 4 3]);
            ym_3d = permute(ym_im(:,:,:,ir), [1 2 4 3]);

            % joint estimation
            tic;
            [M0s_ml(:,:,ir,ip), M0s_med(:,:,ir,ip), M0s_reg(:,:,ir,ip),...
                T1_reg(:,:,ir,ip), T2_reg(:,:,ir,ip), ~, wf_reg(:,:,ir,ip), ~] = ...
                mri_spgr_dess_m0st1t2kap_map(...
                T1_med(:,:,ir,ip), T2_med(:,:,ir,ip),...
                kap_reg(:,:,ir), wf_ml(:,:,ir,ip),...
                flip_s, flip_d, ys_3d, yp_3d, ym_3d, loose_mask,...
                T1max, T1min, T2max, T2min, TRs, TRd, TEs, TEd,...
                Rm, R1, R2, Rk, n_outer, niterM, niter1, niter2, niterK,...
                tolM, tol1, tol2, tolk, is_mag, vi);
            time_reg(ir,ip) = toc;
        end
    end
    
    % bulk post-processing for display
    reg_mask = repmat(tight_mask, [1 1 nr nprof]);
    M0s_ml(~reg_mask) = 0;  T1_ml(~reg_mask) = 0;  T2_ml(~reg_mask) = 0;  
    M0s_med(~reg_mask) = 0; T1_med(~reg_mask) = 0; T2_med(~reg_mask) = 0; 
    M0s_reg(~reg_mask) = 0; T1_reg(~reg_mask) = 0; T2_reg(~reg_mask) = 0;
    
    % save maps for future use
    if (sv)
        dir = pwd;
        cd('HPD_Phantom_05,16,16');
        save(sprintf('maps_slice%u.mat', sl+2), 'M0s_*', 'T1_*', 'T2_*',...
            'fudge_factor', 'kap_scale', 'time_*');
        cd(dir);
    end
end


%% display maps
suffix = 'HPD';

% m0s images 
m0smin = 0; m0smax = 8; 
m0s_rng = [m0smin m0smax];
figure; im('notick', 'row', nprof, 'col', nr, abs(M0s_ml), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_ml_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, abs(M0s_med), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_med_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, abs(M0s_reg), m0s_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('M0s_reg_', suffix, '.eps')), end;

% t1 images
% t1min = 0; t1max = 2000; 
t1min = 500; t1max = 1500;
t1_rng = [t1min t1max];
figure; im('notick', 'row', nprof, 'col', nr, T1_ml, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_ml_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, T1_med, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_med_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, T1_reg, t1_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T1_reg_', suffix, '.eps')), end;

% t2 images
% t2min = 0; t2max = 500; 
t2min = 20; t2max = 120;
t2_rng = [t2min t2max]; 
figure; im('notick', 'row', nprof, 'col', nr, T2_ml, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_ml_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, T2_med, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_med_', suffix, '.eps')), end;
figure; im('notick', 'row', nprof, 'col', nr, T2_reg, t2_rng, 'cbar', ' ');
if (pr), print('-deps', strcat('T2_reg_', suffix, '.eps')), end;

% kap images
kapmin = 0; kapmax = 1;
figure; im('notick', 'row', 1, 'col', nr, kap_reg, [kapmin kapmax], 'cbar', ' ');
if (pr), print('-deps', strcat('kap_reg_', suffix)), end;


%% compute first and second order statistics
% sample mean maps
T1_ml_mean  = squeeze(mean(T1_ml, 3));
T1_med_mean = squeeze(mean(T1_med, 3));
T1_reg_mean = squeeze(mean(T1_reg, 3));

T2_ml_mean  = squeeze(mean(T2_ml, 3));
T2_med_mean = squeeze(mean(T2_med, 3));
T2_reg_mean = squeeze(mean(T2_reg, 3));

% sample standard deviation maps
T1_ml_std   = squeeze(std(T1_ml, 0, 3));
T1_med_std  = squeeze(std(T1_med, 0, 3));
T1_reg_std  = squeeze(std(T1_reg, 0, 3));

T2_ml_std   = squeeze(std(T2_ml, 0, 3));
T2_med_std  = squeeze(std(T2_med, 0, 3));
T2_reg_std  = squeeze(std(T2_reg, 0, 3));

% log10 sample mean maps
logT1_ml_mean   = squeeze(mean(log10(T1_ml), 3));
logT1_med_mean  = squeeze(mean(log10(T1_med), 3));
logT1_reg_mean  = squeeze(mean(log10(T1_reg), 3));

logT2_ml_mean   = squeeze(mean(log10(T2_ml), 3));
logT2_med_mean  = squeeze(mean(log10(T2_med), 3));
logT2_reg_mean  = squeeze(mean(log10(T2_reg), 3));

% log10 sample standard deviation maps
logT1_ml_std    = squeeze(std(log10(T1_ml), 0, 3));
logT1_med_std   = squeeze(std(log10(T1_med), 0, 3));
logT1_reg_std   = squeeze(std(log10(T1_reg), 0, 3));

logT2_ml_std    = squeeze(std(log10(T2_ml), 0, 3));
logT2_med_std   = squeeze(std(log10(T2_med), 0, 3));
logT2_reg_std   = squeeze(std(log10(T2_reg), 0, 3));

% pool sample statistics over rois
T1_ml_mean_roi  = NaN(nROI, nprof);
T1_med_mean_roi = NaN(nROI, nprof);
T1_reg_mean_roi = NaN(nROI, nprof);

T2_ml_mean_roi  = NaN(nROI, nprof);
T2_med_mean_roi = NaN(nROI, nprof);
T2_reg_mean_roi = NaN(nROI, nprof);

T1_ml_std_roi   = NaN(nROI, nprof);
T1_med_std_roi  = NaN(nROI, nprof);
T1_reg_std_roi  = NaN(nROI, nprof);

T2_ml_std_roi   = NaN(nROI, nprof);
T2_med_std_roi  = NaN(nROI, nprof);
T2_reg_std_roi  = NaN(nROI, nprof);

logT1_ml_mean_roi  = NaN(nROI, nprof);
logT1_med_mean_roi = NaN(nROI, nprof);
logT1_reg_mean_roi = NaN(nROI, nprof);

logT2_ml_mean_roi  = NaN(nROI, nprof);
logT2_med_mean_roi = NaN(nROI, nprof);
logT2_reg_mean_roi = NaN(nROI, nprof);

logT1_ml_std_roi   = NaN(nROI, nprof);
logT1_med_std_roi  = NaN(nROI, nprof);
logT1_reg_std_roi  = NaN(nROI, nprof);

logT2_ml_std_roi   = NaN(nROI, nprof);
logT2_med_std_roi  = NaN(nROI, nprof);
logT2_reg_std_roi  = NaN(nROI, nprof);

for ir = 1:nROI
    for ip = 1:nprof
        T1_ml_mean_roi(ir, ip)  = mean(masker(T1_ml_mean(:,:,ip), roi_mask(:,:,ir)));
        T1_med_mean_roi(ir, ip) = mean(masker(T1_med_mean(:,:,ip), roi_mask(:,:,ir)));
        T1_reg_mean_roi(ir, ip) = mean(masker(T1_reg_mean(:,:,ip), roi_mask(:,:,ir)));
        
        T2_ml_mean_roi(ir, ip)  = mean(masker(T2_ml_mean(:,:,ip), roi_mask(:,:,ir)));
        T2_med_mean_roi(ir, ip) = mean(masker(T2_med_mean(:,:,ip), roi_mask(:,:,ir)));
        T2_reg_mean_roi(ir, ip) = mean(masker(T2_reg_mean(:,:,ip), roi_mask(:,:,ir)));
        
        T1_ml_std_roi(ir, ip)   = mean(masker(T1_ml_std(:,:,ip), roi_mask(:,:,ir)));
        T1_med_std_roi(ir, ip)  = mean(masker(T1_med_std(:,:,ip), roi_mask(:,:,ir)));
        T1_reg_std_roi(ir, ip)  = mean(masker(T1_reg_std(:,:,ip), roi_mask(:,:,ir)));
        
        T2_ml_std_roi(ir, ip)   = mean(masker(T2_ml_std(:,:,ip), roi_mask(:,:,ir)));
        T2_med_std_roi(ir, ip)  = mean(masker(T2_med_std(:,:,ip), roi_mask(:,:,ir)));
        T2_reg_std_roi(ir, ip)  = mean(masker(T2_reg_std(:,:,ip), roi_mask(:,:,ir)));
        
        logT1_ml_mean_roi(ir, ip)   = mean(masker(logT1_ml_mean(:,:,ip), roi_mask(:,:,ir)));
        logT1_med_mean_roi(ir, ip)  = mean(masker(logT1_med_mean(:,:,ip), roi_mask(:,:,ir)));
        logT1_reg_mean_roi(ir, ip)  = mean(masker(logT1_reg_mean(:,:,ip), roi_mask(:,:,ir)));
        
        logT2_ml_mean_roi(ir, ip)   = mean(masker(logT2_ml_mean(:,:,ip), roi_mask(:,:,ir)));
        logT2_med_mean_roi(ir, ip)  = mean(masker(logT2_med_mean(:,:,ip), roi_mask(:,:,ir)));
        logT2_reg_mean_roi(ir, ip)  = mean(masker(logT2_reg_mean(:,:,ip), roi_mask(:,:,ir)));
        
        logT1_ml_std_roi(ir, ip)    = mean(masker(logT1_ml_std(:,:,ip), roi_mask(:,:,ir)));
        logT1_med_std_roi(ir, ip)   = mean(masker(logT1_med_std(:,:,ip), roi_mask(:,:,ir)));
        logT1_reg_std_roi(ir, ip)   = mean(masker(logT1_reg_std(:,:,ip), roi_mask(:,:,ir)));
        
        logT2_ml_std_roi(ir, ip)    = mean(masker(logT2_ml_std(:,:,ip), roi_mask(:,:,ir)));
        logT2_med_std_roi(ir, ip)   = mean(masker(logT2_med_std(:,:,ip), roi_mask(:,:,ir)));
        logT2_reg_std_roi(ir, ip)   = mean(masker(logT2_reg_std(:,:,ip), roi_mask(:,:,ir)));
    end
end


%% worst-case sample ml estimator standard deviations over vials of interest
% define vials of interest
roi_t = 5:7;
roi_b = 4:8;
vials = 1:nROI;

% relative scaling of parameters
scale.t1 = 0.1;
scale.t2 = 1;
Psi_ml_std_roi = scale.t1 * T1_ml_std_roi + scale.t2 * T2_ml_std_roi;

% tight-range worst-case values
sigwt_t1 = max(T1_ml_std_roi(roi_t,:), [], 1);
sigwt_t2 = max(T2_ml_std_roi(roi_t,:), [], 1);
Psiwt = max(Psi_ml_std_roi(roi_t,:), [], 1);

% broad-range worst-case values
sigwb_t1 = max(T1_ml_std_roi(roi_b,:), [], 1);
sigwb_t2 = max(T2_ml_std_roi(roi_b,:), [], 1);
Psiwb = max(Psi_ml_std_roi(roi_b,:), [], 1);


%% performance summary 
% print statistics to display
fprintf('\n\nT1 ML Summary Statistics\n');
for ir = 1:nROI
    fprintf('Vial%u:\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\n', ir,...
        T1_ml_mean_roi(ir, 1), char(177), T1_ml_std_roi(ir, 1),...
        T1_ml_mean_roi(ir, 2), char(177), T1_ml_std_roi(ir, 2),...
        T1_ml_mean_roi(ir, 3), char(177), T1_ml_std_roi(ir, 3));
end
            
fprintf('\n\nT1 RLS Summary Statistics\n');
for ir = 1:nROI
    fprintf('Vial%u:\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\n', ir,...
        T1_reg_mean_roi(ir, 1), char(177), T1_reg_std_roi(ir, 1),...
        T1_reg_mean_roi(ir, 2), char(177), T1_reg_std_roi(ir, 2),...
        T1_reg_mean_roi(ir, 3), char(177), T1_reg_std_roi(ir, 3));
end
    
fprintf('\n\nT2 ML Summary Statistics\n');
for ir = 1:nROI
    fprintf('Vial%u:\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\n', ir,...
        T2_ml_mean_roi(ir, 1), char(177), T2_ml_std_roi(ir, 1),...
        T2_ml_mean_roi(ir, 2), char(177), T2_ml_std_roi(ir, 2),...
        T2_ml_mean_roi(ir, 3), char(177), T2_ml_std_roi(ir, 3));
end
            
fprintf('\n\nT2 RLS Summary Statistics\n');
for ir = 1:nROI
    fprintf('Vial%u:\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\t%7.2f\t%c%5.2f\n', ir,...
        T2_reg_mean_roi(ir, 1), char(177), T2_reg_std_roi(ir, 1),...
        T2_reg_mean_roi(ir, 2), char(177), T2_reg_std_roi(ir, 2),...
        T2_reg_mean_roi(ir, 3), char(177), T2_reg_std_roi(ir, 3));
end

% print profile comparison
fprintf('\n\nWorst-case performance summary over tight t1/t2 ranges\n');
fprintf('sigwt1\t%7.2f\t%7.2f\t%7.2f\n', sigwt_t1);
fprintf('sigwt2\t%7.2f\t%7.2f\t%7.2f\n', sigwt_t2);
fprintf('Psiwt\t%7.2f\t%7.2f\t%7.2f\n', Psiwt);

fprintf('\n\nWorst-case performance summary over broad t1/t2 ranges\n');
fprintf('sigwb1\t%7.2f\t%7.2f\t%7.2f\n', sigwb_t1);
fprintf('sigwb2\t%7.2f\t%7.2f\t%7.2f\n', sigwb_t2);
fprintf('Psiwb\t%7.2f\t%7.2f\t%7.2f\n', Psiwb);

%% comparison with nist values
T1ideal = 0:10000;
T2ideal = 0:1000;
labels = cellstr(num2str([1:nROI]'));

% nist t1/t2 mean/std dev
T1_nist_mean_roi = ...
    [2480   2173   1907   1604    1332    1044    801.7  608.6  458.4  336.5  244.2   176.6   126.9   90.9]';
T1_nist_std_roi = ...
    [  10.8   14.7   10.3    7.2     0.8     3.2    1.70   1.03   0.33   0.18   0.09    0.09    0.03    0.05]'; 

T2_nist_mean_roi = ...
    [ 581.3  403.5  278.1  190.94  133.27   96.89  64.07  46.42  31.97  22.56  15.813  11.237   7.911  5.592]';
T2_nist_std_roi = ...
    [   0.39   0.55   0.28   0.011   0.073   0.049  0.034  0.014  0.083  0.012  0.0061  0.0057  0.0037 0.0055]';

% log-transformation 
% not quite right... better to rescale data but impossible for nist
T1_nist_cfvar       = T1_nist_mean_roi ./ T1_nist_std_roi;
logT1_nist_mean_roi = log10(T1_nist_mean_roi);
logT1_nist_std_roi  = logT1_nist_mean_roi ./ T1_nist_cfvar;

T2_nist_cfvar       = T2_nist_mean_roi ./ T2_nist_std_roi;
logT2_nist_mean_roi = log10(T2_nist_mean_roi);
logT2_nist_std_roi  = logT2_nist_mean_roi ./ T2_nist_cfvar;

% t1/t2 roi fill areas
T1min_roi = 800;   
T1max_roi = 1400;      
t1_xbox_roi = log10([T1min_roi T1min_roi T1max_roi T1max_roi]);
t1_ybox_roi = log10([T1min_roi T1max_roi T1max_roi T1min_roi]);

T2min_roi = 50;    
T2max_roi = 120;
t2_xbox_roi = log10([T2min_roi T2min_roi T2max_roi T2max_roi]);
t2_ybox_roi = log10([T2min_roi T2max_roi T2max_roi T2min_roi]);

% t1/t2 robust-range fill areas
T1min_rob = 400;   
T1max_rob = 2000;      
t1_xbox_rob = log10([T1min_rob T1min_rob T1max_rob T1max_rob]);
t1_ybox_rob = log10([T1min_rob T1max_rob T1max_rob T1min_rob]);

T2min_rob = 40;    
T2max_rob = 200;
t2_xbox_rob = log10([T2min_rob T2min_rob T2max_rob T2max_rob]);
t2_ybox_rob = log10([T2min_rob T2max_rob T2max_rob T2min_rob]);

% define colors
o = [1 0.6 0.2];        % Orange
f = [34 139 34]/256;    % ForestGreen

% figure 1: t1 ml from spgr/dess vs. nist values
figure; hold on;
fill(t1_xbox_rob, t1_ybox_rob, 'y');
fill(t1_xbox_roi, t1_ybox_roi, o);
t1_ideal = plot(log10(T1ideal), log10(T1ideal), 'k--', 'LineWidth', 1);
t1_ml_21 = errorbarxy(logT1_nist_mean_roi, logT1_ml_mean_roi(:,1),...
    logT1_nist_std_roi, logT1_ml_std_roi(:,1), 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'c');
t1_ml_11 = errorbarxy(logT1_nist_mean_roi, logT1_ml_mean_roi(:,2),...
    logT1_nist_std_roi, logT1_ml_std_roi(:,2), 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'v', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'm');
t1_ml_02 = errorbarxy(logT1_nist_mean_roi, logT1_ml_mean_roi(:,3),...
    logT1_nist_std_roi, logT1_ml_std_roi(:,3), 'Color', f,...
    'LineStyle', 'none', 'Marker', 'd', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'g');
txton = 1:nROI;
text(logT1_nist_mean_roi(txton), logT1_ml_mean_roi(txton,1),...
    labels(txton), 'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'right', 'FontSize', 16);

hold off; axis tight; axis square; grid on;
set(gca, 'XLim', [1.9 3.5]);
set(gca, 'XTick', [1.9:0.2:3.5]);
set(gca, 'YLim', [1.9 3.5]);
set(gca, 'YTick', [1.9:0.2:3.5]);
xlabel('NIST log_1_0(T1) Estimates (log_1_0 ms)', 'FontSize', 16);
xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
ylabel('SPGR/DESS log_1_0(T1) ML Estimates (log_1_0 ms)', 'FontSize', 16);

leg_t1ml = legend([t1_ideal, t1_ml_21(1), t1_ml_11(1), t1_ml_02(1)],...
    'Ideal', '(0,2)', '(1,1)', '(2,1)', 'Location', 'NW');
set(leg_t1ml, 'FontSize', 18);
if (pr), print('-depsc', 't1_ml_compare_HPD.eps'); end;

% figure 2: t2 ml from spgr/dess vs. nist values
figure; hold on;
fill(t2_xbox_rob, t2_ybox_rob, 'y');
fill(t2_xbox_roi, t2_ybox_roi, o);
t2_ideal = plot(log10(T2ideal), log10(T2ideal), 'k--', 'LineWidth', 1);
t2_ml_21 = errorbarxy(logT2_nist_mean_roi, logT2_ml_mean_roi(:,1),...
    logT2_nist_std_roi, logT2_ml_std_roi(:,1), 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'c');
t2_ml_11 = errorbarxy(logT2_nist_mean_roi, logT2_ml_mean_roi(:,2),...
    logT2_nist_std_roi, logT2_ml_std_roi(:,2), 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'v', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'm');
t2_ml_02 = errorbarxy(logT2_nist_mean_roi, logT2_ml_mean_roi(:,3),...
    logT2_nist_std_roi, logT2_ml_std_roi(:,3), 'Color', f,...
    'LineStyle', 'none', 'Marker', 'd', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'g');
txton = 1:12;
text(logT2_nist_mean_roi(txton), logT2_ml_mean_roi(txton,1),...
    labels(txton), 'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'right', 'FontSize', 16);

hold off; axis tight; axis square; grid on;
set(gca, 'XLim', [1 2.8]);
set(gca, 'XTick', [1:0.3:2.8]);
set(gca, 'YLim', [1 2.8]);
set(gca, 'YTick', [1:0.3:2.8]);
xlabel('NIST log_1_0(T2) Estimates (log_1_0 ms)', 'FontSize', 16);
xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
ylabel('SPGR/DESS log_1_0(T2) ML Estimates (log_1_0 ms)', 'FontSize', 16);

leg_t2ml = legend([t2_ideal, t2_ml_21(1), t2_ml_11(1), t2_ml_02(1)],...
    'Ideal', '(0,2)', '(1,1)', '(2,1)', 'Location', 'NW');
set(leg_t2ml, 'FontSize', 18);
if (pr), print('-depsc', 't2_ml_compare_HPD.eps'); end;

% figure 3: t1 rls from spgr/dess vs. nist values
figure; hold on;
fill(t1_xbox_rob, t1_ybox_rob, 'y');
fill(t1_xbox_roi, t1_ybox_roi, o);
t1_ideal = plot(log10(T1ideal), log10(T1ideal), 'k--', 'LineWidth', 1);
t1_reg_21 = errorbarxy(logT1_nist_mean_roi, logT1_reg_mean_roi(:,1),...
    logT1_nist_std_roi, logT1_reg_std_roi(:,1), 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'c');
t1_reg_11 = errorbarxy(logT1_nist_mean_roi, logT1_reg_mean_roi(:,2),...
    logT1_nist_std_roi, logT1_reg_std_roi(:,2), 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'v', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'm');
t1_reg_02 = errorbarxy(logT1_nist_mean_roi, logT1_reg_mean_roi(:,3),...
    logT1_nist_std_roi, logT1_reg_std_roi(:,3), 'Color', f,...
    'LineStyle', 'none', 'Marker', 'd', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'g');
txton = 1:nROI;
text(logT1_nist_mean_roi(txton), logT1_reg_mean_roi(txton,1),...
    labels(txton), 'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'right', 'FontSize', 16);

hold off; axis tight; axis square; grid on;
set(gca, 'XLim', [1.9 3.5]);
set(gca, 'XTick', [1.9:0.2:3.5]);
set(gca, 'YLim', [1.9 3.5]);
set(gca, 'YTick', [1.9:0.2:3.5]);
xlabel('NIST log_1_0(T1) Estimates (log_1_0 ms)', 'FontSize', 16);
xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
ylabel('SPGR/DESS log_1_0(T1) RLS Estimates (log_1_0 ms)', 'FontSize', 16);

leg_t1reg = legend([t1_ideal, t1_reg_21(1), t1_reg_11(1), t1_reg_02(1)],...
    'Ideal', '(0,2)', '(1,1)', '(2,1)', 'Location', 'NW');
set(leg_t1reg, 'FontSize', 18);
if (pr), print('-depsc', 't1_reg_compare_HPD.eps'); end;

% figure 4: t2 rls from spgr/dess vs. nist values
figure; hold on;
fill(t2_xbox_rob, t2_ybox_rob, 'y');
fill(t2_xbox_roi, t2_ybox_roi, o);
t2_ideal = plot(log10(T2ideal), log10(T2ideal), 'k--', 'LineWidth', 1);
t2_reg_21 = errorbarxy(logT2_nist_mean_roi, logT2_reg_mean_roi(:,1),...
    logT2_nist_std_roi, logT2_reg_std_roi(:,1), 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'c');
t2_reg_11 = errorbarxy(logT2_nist_mean_roi, logT2_reg_mean_roi(:,2),...
    logT2_nist_std_roi, logT2_reg_std_roi(:,2), 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'v', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'm');
t2_reg_02 = errorbarxy(logT2_nist_mean_roi, logT2_reg_mean_roi(:,3),...
    logT2_nist_std_roi, logT2_reg_std_roi(:,3), 'Color', f,...
    'LineStyle', 'none', 'Marker', 'd', 'LineWidth', 1.0,...
    'MarkerSize', 8, 'MarkerFaceColor', 'g');
txton = 1:12;
text(logT2_nist_mean_roi(txton), logT2_reg_mean_roi(txton,1),...
    labels(txton), 'VerticalAlignment', 'bottom',...
    'HorizontalAlignment', 'right', 'FontSize', 16);

hold off; axis tight; axis square; grid on;
set(gca, 'XLim', [1 2.8]);
set(gca, 'XTick', [1:0.3:2.8]);
set(gca, 'YLim', [1 2.8]);
set(gca, 'YTick', [1:0.3:2.8]);
xlabel('NIST log_1_0(T2) Estimates (log_1_0 ms)', 'FontSize', 16);
xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
ylabel('SPGR/DESS log_1_0(T2) RLS Estimates (log_1_0 ms)', 'FontSize', 16);

leg_t2reg = legend([t2_ideal, t2_reg_21(1), t2_reg_11(1), t2_reg_02(1)],...
    'Ideal', '(0,2)', '(1,1)', '(2,1)', 'Location', 'NW');
set(leg_t2reg, 'FontSize', 18);
if (pr), print('-depsc', 't2_reg_compare_HPD.eps'); end;
