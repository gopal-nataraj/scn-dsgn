% script HPDphant_2016_05_20_v1.m
% compare t1/t2 estimates from optimized pairs of (2,1), (1,1), (0,2)
%   (spgr, dess) scans, from several independent scan repetitions
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   2016-05-20      hpd phantom data acquired, with lower receive gains
%   2016-05-21      original recon script 
%   2016-05-23      added code to examine noise statistics

%% initial setup
% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

addpath('../../data/HPD_Phantom_05,20,16/');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../crb');
addpath('../../etc');

% global imaging parameters
nx = 256; ny = 256; nz = 8;
FOVx = 240; FOVy = 240; FOVz = 40;          % mm
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz;                        % mm^3

coil.wghts = [1 1 1 1 1 1 1 1]';            % coil weights (0-wght coils may be clipped)
coil.iter = 10;                             % coil combination outer iteration count
coil.log2b = -5;                            % log2(coil reg parameter)
coil.arg = {...
    'coilweights', coil.wghts,...
    'nouter', coil.iter,...
    'log2b', coil.log2b};

nc = 8;                                     % number of receive coils
nr = 10;                                    % number of repetitions
nprof = 3;                                  % number of profiles to test
ncycles = 2;                                % number of dess spoiling cycles

TE_global = 4.67;                           % ms
wf = 0;                                     % off-resonance

sl = 2;                                     % 4th slice selected ([3:6] only extracted)
shift = [0 0];                              % [x y] pixels to circshift
vi = 0;                                     % toggle viewing reg images on/off
ni = 0;                                     % toggle viewing noise images on/off;
pr = 0;                                     % toggle printing on/off 
sv = 0;                                     % toggle data saving on/off

%% b1 mapping from bloch-siegert data
if (~exist('kap_reg', 'var'))
    try
        addpath('../../data/HPD_Phantom_05,20,16/BS/');
        load('B1_reg');
    catch
        % compute bloch-siegert constant
        [~, mag_p, ph_p] = readwav('bsp.wav');
        [~, mag_m, ph_m] = readwav('bsm.wav');
        w_rf = 8000;                            % Hz        (now set using rffreq instead of iath)
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
        cd('/Volumes/a2/gnataraj/Data/2016,05,20_hpd_8ch/');
        for ir = 1:nr
            load(sprintf('im,hpd-phant,bs-spgr-dess,pfreq,rep%u.mat', ir));
            ybs_p_coil(:,:,:,ir) = ims{8}.echo1(:,:,sl,:);
            ybs_m_coil(:,:,:,ir) = ims{9}.echo1(:,:,sl,:);
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
        b1_peak = 0.075;                        % Gauss, for 90deg flip
        kap_reg = b1_reg / b1_peak;             % flip angle scale factor
        
        figure; im(b1_init, [0 b1_peak], 'cbar');
        figure, im(b1_reg,  [0 b1_peak], 'cbar');
        figure; im(kap_reg, [0 1], 'cbar');
        
        % Save regularized B1 maps for future use
        dir = pwd; cd('HPD_Phantom_05,20,16/BS/');
        save('B1_reg', 'b1_peak', 'b1_init', 'b1_reg', 'kap_reg'); cd(dir);
    end
    
    % DEBUG: try fudging the scale factor (default = 1)
    fudge_factor = 1.016/0.913;
    kap_reg = fudge_factor * kap_reg;
end


%% global constant instantiation
% % define rois
ctrX = [ 74  84 111 144 171 181 172 145 112  85 107 106 149 149];
ctrY = [130 161 181 181 162 131 100  80  79  99 109 151 152 109];
nROI = length(ctrX);
rad = 3 * ones(1, nROI);
% rad = 5 * ones(1, nROI);

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

% noise statistics
noise = cell(nprof,1);

try
    % first try to load maps
    fprintf('\nTrying to load parameter maps...');
    load(sprintf('maps_slice%u.mat', sl+2));
    fprintf('success!\n');
catch
    %% if unavailable, make maps for: (2,1), (1,1), or (0,2)
    fprintf('failed.\n');
    for ip = 1:nprof
        switch ip
            case 1
                n_spgr = 2;
                n_dess = 1;
                flip_s_deg = [15, 5]';              % degrees
                flip_d_deg = [30]';                 % degrees
                TRs = [12.2, 12.2]';                % ms
                TRd = [17.5]';                      % ms
                mask.ir = 6;                        % mask design rep
                mask.thresh = 0.05;                 % mask design threshold
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
                n_spgr = 0;
                n_dess = 2;
                flip_s_deg = []';                   % degrees
                flip_d_deg = [35, 10]';             % degrees
                TRs = []';                          % ms
                TRd = [24.4, 17.5]';                % ms
                mask.ir = 8;                        % mask design rep
                mask.thresh = 0.04;                 % mask design threshold
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
            addpath(sprintf('HPD_Phantom_05,20,16/%uSPGR,%uDESS', n_spgr, n_dess));
            load(sprintf('ims_coil_comb_%uSPGR%uDESS_slice%u.mat', n_spgr, n_dess, sl+2));
            fprintf('success!\n');
        catch
            fprintf('failed.\n');
            ys_coil = NaN(nx, ny, nfs, nc, nr);
            yp_coil = NaN(nx, ny, nfd, nc, nr);
            ym_coil = NaN(nx, ny, nfd, nc, nr);
            
            % if not found extract data from all repetitions
            dir = pwd;
            cd('/Volumes/a2/gnataraj/Data/2016,05,20_hpd_8ch/');
            for ir = 1:nr
                fprintf('Loading coil data for (%u,%u) profile, repetition %u.\n',...
                    n_spgr, n_dess, ir);
                load(sprintf('im,hpd-phant,bs-spgr-dess,pfreq,rep%u.mat', ir));
                switch ip
                    case 1
                        ys_coil(:,:,1,:,ir) = ims{5}.echo1(:,:,sl,:);
                        ys_coil(:,:,2,:,ir) = ims{6}.echo1(:,:,sl,:);
                        yp_coil(:,:,1,:,ir) = ims{2}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{2}.echo2(:,:,sl,:);
                    case 2
                        ys_coil(:,:,1,:,ir) = ims{7}.echo1(:,:,sl,:);
                        yp_coil(:,:,1,:,ir) = ims{3}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{3}.echo2(:,:,sl,:);
                    case 3
                        yp_coil(:,:,1,:,ir) = ims{1}.echo1(:,:,sl,:);
                        ym_coil(:,:,1,:,ir) = ims{1}.echo2(:,:,sl,:);
                        yp_coil(:,:,2,:,ir) = ims{4}.echo1(:,:,sl,:);
                        ym_coil(:,:,2,:,ir) = ims{4}.echo2(:,:,sl,:);
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
            cd(sprintf('HPD_Phantom_05,20,16/%uSPGR,%uDESS', n_spgr, n_dess));
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
        T1.d = logspace(2, 4, 300);
        T2.d = logspace(1, 3, 300);
        D = NaN(M, length(T1.d), length(T2.d));
        K = length(T1.d) * length(T2.d);
        for t1 = 1:length(T1.d)
            for t2 = 1:length(T2.d)
                % SPGR dictionary component
                for a = 1:nfs
                    D(a, t1, t2) = spgr_fun(1, T1.d(t1), kap_scale,...
                        flip_s(a), TRs(a), TEs(a), wf, is_mag);
                end

                % DESS dictionary component
                for a = 1:nfd
                    [D(nfs+a, t1, t2), D(nfs+nfd+a, t1, t2)] = ...
                        dess_fun(1, T1.d(t1), T2.d(t2), kap_scale,...
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
        
        % sanity check: coil-combined data noise maps
        if nfs>0
            noise{ip}.s.mean    = embed(mean(masker(ys_im, ~loose_mask), 3), ~loose_mask);      % [nx ny nfs]
            noise{ip}.s.std     = embed(std(masker(ys_im, ~loose_mask), 0, 3), ~loose_mask);    % [nx ny nfs]
        end
        if nfd>0
            noise{ip}.p.mean    = embed(mean(masker(yp_im, ~loose_mask), 3), ~loose_mask);      % [nx ny nfd]
            noise{ip}.p.std     = embed(std(masker(yp_im, ~loose_mask), 0, 3), ~loose_mask);    % [nx ny nfd]
            
            noise{ip}.m.mean    = embed(mean(masker(ym_im, ~loose_mask), 3), ~loose_mask);      % [nx ny nfd]
            noise{ip}.m.std     = embed(std(masker(ym_im, ~loose_mask), 0, 3), ~loose_mask);    % [nx ny nfd]
        end
        
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
            [t1_idx, t2_idx] = ind2sub([length(T1.d) length(T2.d)], idx);
            tmp1 = reshape(T1.d(t1_idx), [nx ny]);    % placeholder for cur rep
            tmp2 = reshape(T2.d(t2_idx), [nx ny]);    % placeholder for cur rep

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
        cd('HPD_Phantom_05,20,16');
        save(sprintf('maps_slice%u.mat', sl+2), 'M0s_*', 'T1_*', 'T2_*',...
            'fudge_factor', 'kap_scale', 'time_*', '*mask', 'noise');
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

% noise images
if (ni)
    n.min = 0; n.max = 0.05;
    figure; im('notick', noise{1}.s.mean, [n.min n.max], 'cbar', '(2,1) spgr data noise mean');
    figure; im('notick', noise{2}.s.mean, [n.min n.max], 'cbar', '(1,1) spgr data noise mean');

    figure; im('notick', noise{1}.p.mean, [n.min n.max], 'cbar', '(2,1) dess echo1 data noise mean');
    figure; im('notick', noise{2}.p.mean, [n.min n.max], 'cbar', '(1,1) dess echo1 data noise mean');
    figure; im('notick', noise{3}.p.mean, [n.min n.max], 'cbar', '(0,2) dess echo1 data noise mean');

    figure; im('notick', noise{1}.m.mean, [n.min n.max], 'cbar', '(2,1) dess echo2 data noise mean');
    figure; im('notick', noise{2}.m.mean, [n.min n.max], 'cbar', '(1,1) dess echo2 data noise mean');
    figure; im('notick', noise{3}.m.mean, [n.min n.max], 'cbar', '(0,2) dess echo2 data noise mean');

    figure; im('notick', noise{1}.s.std, [n.min n.max], 'cbar', '(2,1) spgr data noise std dev');
    figure; im('notick', noise{2}.s.std, [n.min n.max], 'cbar', '(1,1) spgr data noise std dev');

    figure; im('notick', noise{1}.p.std, [n.min n.max], 'cbar', '(2,1) dess echo1 data noise std dev');
    figure; im('notick', noise{2}.p.std, [n.min n.max], 'cbar', '(1,1) dess echo1 data noise std dev');
    figure; im('notick', noise{3}.p.std, [n.min n.max], 'cbar', '(0,2) dess echo1 data noise std dev');

    figure; im('notick', noise{1}.m.std, [n.min n.max], 'cbar', '(2,1) dess echo2 data noise std dev');
    figure; im('notick', noise{2}.m.std, [n.min n.max], 'cbar', '(1,1) dess echo2 data noise std dev');
    figure; im('notick', noise{3}.m.std, [n.min n.max], 'cbar', '(0,2) dess echo2 data noise std dev');
end


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
roi.t1.t = 5:6;
roi.t2.t = 6:7;
[roi.psi.t, iroi.t1.t, iroi.t2.t] = intersect(roi.t1.t, roi.t2.t);
roi.t1.b = 3:9;
roi.t2.b = 4:8;
[roi.psi.b, iroi.t1.b, iroi.t2.b] = intersect(roi.t1.b, roi.t2.b);

% relative scaling of parameters
scale.t1 = 0.1;
scale.t2 = 1;
Psi_ml_std_roi = scale.t1 * T1_ml_std_roi + scale.t2 * T2_ml_std_roi;

% tight-range worst-case values
sigwt_t1 = max(T1_ml_std_roi(roi.t1.t,:), [], 1);
sigwt_t2 = max(T2_ml_std_roi(roi.t2.t,:), [], 1);
Psiwt = max(Psi_ml_std_roi(roi.psi.t,:), [], 1);

% broad-range worst-case values
sigwb_t1 = max(T1_ml_std_roi(roi.t1.b,:), [], 1);
sigwb_t2 = max(T2_ml_std_roi(roi.t2.b,:), [], 1);
Psiwb = max(Psi_ml_std_roi(roi.psi.b,:), [], 1);


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
fprintf('\n\nEmpirical performance summary over tight t1/t2 ranges\n');
fprintf('sigwt1\t%7.2f\t%7.2f\t%7.2f\n', sigwt_t1);
fprintf('sigwt2\t%7.2f\t%7.2f\t%7.2f\n', sigwt_t2);
fprintf('Psiwt\t%7.2f\t%7.2f\t%7.2f\n', Psiwt);

fprintf('\n\nEmpirical performance summary over broad t1/t2 ranges\n');
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
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
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
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
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
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
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
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
set(leg_t2reg, 'FontSize', 18);
if (pr), print('-depsc', 't2_reg_compare_HPD.eps'); end;


%% cramer rao bound predictions
% conservative noise mask
noise_mask = ~imdilate(tight_mask, strel('disk', 20));

Sig_inv = cell(nprof,1);
sigwt1 = cell(nprof,1);
sigwt2 = cell(nprof,1);
Psi = cell(nprof,1);

for ip = 1:nprof
  % define profile
  switch ip
    case 1
        nfs = 2;
        nfd = 1;
        flip_s = [15, 5]' * (pi/180);       % rad
        flip_d = [30]' * (pi/180);          % rad
        TRs = [12.2, 12.2]';                % ms
        TRd = [17.5]';                      % ms
    case 2
        nfs = 1;
        nfd = 1;
        flip_s = [15]' * (pi/180);          % rad
        flip_d = [10]' * (pi/180);          % rad
        TRs = [13.9]';                      % ms
        TRd = [28.0]';                      % ms
    case 3
        nfs = 0;
        nfd = 2;
        flip_s = []' * (pi/180);            % rad
        flip_d = [35, 10]' * (pi/180);      % rad
        TRs = []';                          % ms
        TRd = [24.4, 17.5]';                % ms
  end
  
  % scale flip angles
  flip_s_sc = kap_scale * flip_s;
  flip_d_sc = kap_scale * flip_d;
  
  % estimate noise covariance matrix
  if nfs>0
    var.s = mean(masker(noise{ip}.s.std.^2, noise_mask), 1);
  else
    var.s = [];
  end
  if nfd>0
    var.p = mean(masker(noise{ip}.p.std.^2, noise_mask), 1);
    var.m = mean(masker(noise{ip}.m.std.^2, noise_mask), 1);
  else
    var.p = [];
    var.m = [];
  end
%   Sig_inv{ip} = diag(1 ./ [var.s var.p var.m]');
  Sig_inv{ip} = diag((1./1.49e-7) * ones(nfs+2*nfd,1));

  % find worst case t1/t2 std dev and cost over tight range
  sigt1.t = NaN(length(roi.t1.t), length(roi.t1.t));
  sigt2.t = NaN(length(roi.t2.t), length(roi.t2.t));
  T1.t = T1_nist_mean_roi(roi.t1.t);
  T2.t = T2_nist_mean_roi(roi.t2.t);
  
  for it1 = 1:length(T1.t)
    for it2 = 1:length(T2.t)
      [~, ~, sigt1.t(it1,it2), sigt2.t(it1,it2)] = ...
        norm_crlb_dess_3parm(T1.t(it1), T2.t(it2), 0,...
        flip_s_sc, flip_d_sc, TRs, TRd, TE_global, Sig_inv{ip}, 0);
    end
  end
  
  sigwt1{ip}.t = max(col(sigt1.t));
  sigwt2{ip}.t = max(col(sigt2.t));
  Psi{ip}.t = max(col(scale.t1 * sigt1.t(iroi.t1.t, iroi.t2.t) ...
    + scale.t2 * sigt2.t(iroi.t1.t, iroi.t2.t)));
  
  % first worst case t1/t2 std dev and cost over broad range
  sigt1.b = NaN(length(roi.t1.b), length(roi.t1.b));
  sigt2.b = NaN(length(roi.t2.b), length(roi.t2.b));
  T1.b = T1_nist_mean_roi(roi.t1.b);
  T2.b = T2_nist_mean_roi(roi.t2.b);
  
  for it1 = 1:length(T1.b)
    for it2 = 1:length(T2.b)
      [~, ~, sigt1.b(it1,it2), sigt2.b(it1,it2)] = ...
        norm_crlb_dess_3parm(T1.b(it1), T2.b(it2), 0,...
        flip_s_sc, flip_d_sc, TRs, TRd, TE_global, Sig_inv{ip}, 0);
    end
  end
  
  sigwt1{ip}.b = max(col(sigt1.b));
  sigwt2{ip}.b = max(col(sigt2.b));
  Psi{ip}.b = max(col(scale.t1 * sigt1.b(iroi.t1.b, iroi.t2.b) ...
    + scale.t2 * sigt2.b(iroi.t1.b, iroi.t2.b)));
end

% print profile comparison
fprintf('\n\nPredicted performance summary over tight t1/t2 ranges\n');
fprintf('sigwt1\t%7.2f\t%7.2f\t%7.2f\n', sigwt1{1}.t, sigwt1{2}.t, sigwt1{3}.t);
fprintf('sigwt2\t%7.2f\t%7.2f\t%7.2f\n', sigwt2{1}.t, sigwt2{2}.t, sigwt2{3}.t);
fprintf('Psiwt\t%7.2f\t%7.2f\t%7.2f\n', Psi{1}.t, Psi{2}.t, Psi{3}.t);

fprintf('\n\nPredicted performance summary over broad t1/t2 ranges\n');
fprintf('sigwb1\t%7.2f\t%7.2f\t%7.2f\n', sigwt1{1}.b, sigwt1{2}.b, sigwt1{3}.b);
fprintf('sigwb2\t%7.2f\t%7.2f\t%7.2f\n', sigwt2{1}.b, sigwt2{2}.b, sigwt2{3}.b);
fprintf('Psiwb\t%7.2f\t%7.2f\t%7.2f\n', Psi{1}.b, Psi{2}.b, Psi{3}.b);
  