%% Regularized M0/T1 Recon from Synthetic Data: Compare over flips & SNR

%% Synthesis of Phantom SPGR Data
% Load digital phantom
f.filename = 'phantom_1.0mm_msles2_crisp.fld';
cd DigitalPhantom; labels = fld_read(f.filename); cd ..;

% Get true values
slice = 95; [nx, ny, nz] = size(labels);
M0_true = double(squeeze(labels(:,:,slice)));
T1_true = double(squeeze(labels(:,:,slice)));
T2s_true = double(squeeze(labels(:,:,slice)));
for idx = 0:10
    f = mri_brainweb_params(idx);
    M0_true(M0_true == idx) = f.pd;
    T1_true(T1_true == idx) = f.t1;
    T2s_true(T2s_true == idx) = f.t2s;
end

% NRMSE Formula
nrmse = @(tst, tru) sqrt(sum((abs(tst(:))-tru(:)).^2) ./ sum(abs(tru(:)).^2));
fullflip = [5 30]' * pi/180;
fullSNR = [30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]';
err_cube = NaN(4, length(fullflip)-1, length(fullSNR));

% Varying flip angle
for ff = 2:length(fullflip)
    
    % Imaging Parameters
    flip = fullflip(1:ff);
    nf = length(flip);
    TR = 22.1;                              % ms
    TE = 5.1;                               % ms
    E1_true = exp(-TR./T1_true);

    % Forward model: make true data
    y_true = NaN(nx, ny, nf);
    for a = 1:nf
        y_true(:,:,a) = fft2(M0_true .* exp(-TE./T2s_true) .* (1-E1_true)...
            .* sin(flip(a)) ./ (1 - E1_true .* cos(flip(a))));
    end
    
    % Varying SNR
    for ss = 1:length(fullSNR)

        % Add complex white gaussian noise
        SNR = fullSNR(ss);                  % dB
        sigma = exp(-SNR/20) * norm(y_true(:)) / sqrt(2*numel(y_true));
        y = y_true + sigma * (randn(size(y_true)) + 1i * randn(size(y_true)));
        % sigma = exp(-SNR/20) * norm(y_true(:)) / sqrt(numel(y_true));
        % y = y_true + sigma * randn(size(y_true));
        printm('snr = %g dB', 20*log(norm(y_true(:)) / norm(col(y_true-y))));

        % Take noisy data back to image domain
        % Since simulating only body coil, sqrt(SOS) is just the magnitude
        y_im = abs(ifft2(y));

        %% Method-of-Moments Estimation
        E1_map_mom = NaN(nx, ny);
        intercept_map_mom = NaN(nx, ny);
        flipmap = permute(repmat(flip, [1 nx ny]), [2 3 1]);

        for xx = 1:nx
            for yy = 1:ny
                Y = squeeze(y_im(xx,yy,:) ./ sin(flipmap(xx,yy,:)));
                X = squeeze(y_im(xx,yy,:) ./ tan(flipmap(xx,yy,:)));
                W = diag(squeeze((1-cos(flipmap(xx,yy,:)))./sin(flipmap(xx,yy,:)))); 
                A = [X ones(length(flip), 1)];
                regression = (A' * W * A) \ (A' * W * Y);
                E1_map_mom(xx,yy) = regression(1);
                intercept_map_mom(xx,yy) = regression(2);
            end
        end

        M0_map_mom = intercept_map_mom ./ (1-E1_map_mom);
        T1_map_mom = -TR ./ log(abs(E1_map_mom));     % Should be pure real

        % Project Images to within range 
        T1_max = 5000; 
        E1_max = exp(-TR ./ T1_max);
        M0_max = 2;

        T1_map_mom = max(T1_map_mom, 0);
        T1_map_mom = min(T1_map_mom, T1_max);
        M0_map_mom = max(M0_map_mom, 0);
        M0_map_mom = min(M0_map_mom, M0_max);

        %% Preprocessing and Masking
        % Masking
        tight_mask = (M0_map_mom > 0.1) & (T1_map_mom < T1_max) & ...
            (T1_map_mom > 0) & ~isnan(T1_map_mom); 
        loose_mask = imdilate(tight_mask, strel('disk', 5));

        % Set voxels inside loose mask but outside tight mask to mean
        M0_map_mom(~tight_mask & loose_mask) = mean(col(M0_map_mom));
        T1_map_mom(~tight_mask & loose_mask) = mean(col(T1_map_mom));

        M0_map_mom(~loose_mask) = 0;
        T1_map_mom(~loose_mask) = 0;

        % Update the E1 map
        E1_map_mom = exp(-TR ./ T1_map_mom);

        %% Iterative Reconstruction
        % Define iteration parameters
        n_outer_iter = 20;
        niter1 = 200;
        niter2 = 200; 
        tol_1 = 10^-6;
        tol_2 = 10^-6;
        disp = 0;

        % Define regularizers, R1 and R2
        beta1 = 2^-1;%2^-3;
        beta2 = 2^7;%2^5;
        delta1 = 10^-6;%10^-4
        delta2 = 10^-5;
        R1 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
        R2 = Reg1(loose_mask, 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

        % Regularized reconstruction
        y_tmp = permute(y_im, [1 2 4 3]);
        flipmap_tmp = permute(flipmap, [1 2 4 3]);
        mask_tmp = permute(permute(loose_mask, [3 2 1]), [3 2 1]);

        [M0_map, T1_map, cost] = mri_despot1_m0t1map(M0_map_mom, E1_map_mom,...
            y_tmp, flipmap_tmp, loose_mask, R1, R2, T1_max, TR, n_outer_iter, ...
            niter1, niter2, tol_1, tol_2, disp);

        %% Postprocessing for display and analysis
        
        % T1 values in voxels with M0 = 0 should also be 0
        thresh = 0.1;
        T1_map_mom(M0_map_mom < thresh * max(col(M0_map_mom))) = 0;
        T1_map(M0_map < thresh * max(col(M0_map))) = 0;
        
        % Only examine relevant voxels; mask out rest
        post_mask = labels(:,:,slice) > 0; 
        M0_map_mom(~post_mask) = 0; M0_map(~post_mask) = 0; 
        T1_map_mom(~post_mask) = 0; T1_map(~post_mask) = 0;

%         % Images
%         figure; im(M0_map_mom, [0 1], 'cbar'); 
%         figure; im(M0_map, [0 1], 'cbar');
%         figure; im(M0_true, [0 1], 'cbar');
%         figure; im(T1_map_mom, [0 4000], 'cbar');
%         figure; im(T1_map, [0 4000], 'cbar');
%         figure; im(T1_true, [0 4000], 'cbar');

        %% Analysis
        err_cube(1,ff-1,ss) = nrmse(M0_map_mom, M0_true);
        err_cube(2,ff-1,ss) = nrmse(M0_map, M0_true);
        err_cube(3,ff-1,ss) = nrmse(T1_map_mom, T1_true);
        err_cube(4,ff-1,ss) = nrmse(T1_map, T1_true);
    end
end
save('nrmse_digiphant_moreSNR', 'err_cube', 'fullflip', 'fullSNR');

figure(1); hold on;
plot(fullSNR, squeeze(err_cube(1,1,:)), 'b--', 'Linewidth', 2);
plot(fullSNR, squeeze(err_cube(2,1,:)), 'b-', 'Linewidth', 2);
plot(fullSNR, squeeze(err_cube(3,1,:)), 'r--', 'Linewidth', 2);
plot(fullSNR, squeeze(err_cube(4,1,:)), 'r-', 'Linewidth', 2);
hold off; grid on; 
title('NRMSE of M0 and T1 vs. SNR','FontSize',12); 
xlabel('SNR (dB)','FontSize',12); ylabel('NRMSE','FontSize',12);
legend('M0*, conventional', 'M0*, MBET1', 'T1, conventional', ...
    'T1, MBET1','FontSize',12);
% print -depsc nrmse_vs_snr.eps