% Sanity Check Noise Analysis on Synthetic Phantom - Four Parameters
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

%% Raw data creation 
% Header files and IRT Setup
if (~exist('irtdir', 'var'))
    cd ../../IRT; setup(); cd ../Scripts/CRLB_Analysis;
end
addpath('../../Scripts/');
addpath('../Calculations/DESS/');

% Imaging Parameters
nx = 64; ny = 64; nz = 1; nc = 9;
% opflip = [3 5 8 11 15 20 29 37:10:77]'; % degrees
opflip = [11 37 77]';                   % degrees
flip = opflip * pi/180; 
nf = length(flip);
TR = 15.1;                              % ms
TE = 5.04;                              % ms
wf = 0;                                 % Zero off-resonance
rcoil = 100;                          % mm

% Create object mask
phant_rad = nx/2;
[xx, yy] = ndgrid(-nx/2:nx/2-1, -ny/2:ny/2-1);
mask = (xx.^2 + yy.^2) <= phant_rad.^2;

% True parameters 
T1_true = 550;                          % ms
T2_true = 40;                           % ms
kap_true = 1.3;                         % a.u.
M0s_true = 0.7*mask;                    % a.u. (sqrt-SoS scaled up by sqrt(nc))

% % Create (single coil) noiseless phantom data
% yp_true = NaN(nx, ny, nf);
% ym_true = NaN(nx, ny, nf);
% for a = 1:nf
%     [S1true, S2true] = dess_fun(TR, TE, M0s_true, T1_true,...
%         T2_true, kap_true*flip(a), wf);
%     yp_true(:,:,a) = fft2(S1true);
%     ym_true(:,:,a) = fft2(S2true);
% end

% For multicoil data, create sensitivity maps
smap = mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', nc, 'rcoil', rcoil);

% Create multi-coil noiseless phantom data
yp_true = NaN(nx, ny, nc, nf);
ym_true = NaN(nx, ny, nc, nf);
for a = 1:nf
    [S1true, S2true] = dess_fun(TR, TE, M0s_true, T1_true,...
        T2_true, kap_true*flip(a), wf);
    for c = 1:nc
        yp_true(:,:,c,a) = fft2(smap(:,:,c) .* S1true);
        ym_true(:,:,c,a) = fft2(smap(:,:,c) .* S2true);
    end
end

% Add complex white gaussian noise (compensate for FFT operation)
var_im = 1.31e-7; sigma_im = sqrt(var_im); 
var_k = var_im * (nx*ny); sigma_k = sqrt(var_k);
yp = yp_true + sigma_k * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
ym = ym_true + sigma_k * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
printm('snr_p = %g dB', 20*log(norm(yp_true(:)) / norm(col(yp_true-yp))));
printm('snr_m = %g dB', 20*log(norm(ym_true(:)) / norm(col(ym_true-ym))));

% Take noisy data back to image domain
% Since simulating only body coil, sqrt(SOS) is just the magnitude
yp_im = ifft2(yp);
ym_im = ifft2(ym);

% Compute the square-root sum-of-squares raw images
yp_sos = squeeze(sqrt(sum(abs(yp_im).^2, 3)));
ym_sos = squeeze(sqrt(sum(abs(ym_im).^2, 3)));

% %% PART ONE: Compute empirical noise variance
% % Choose a noise-only region
% noise_msk = repmat(~imdilate(mask, strel('disk', 5)), [1 1 nf]);
% 
% % Extract and vectorize noise
% noise_p = reshape(yp_im(noise_msk), [length(yp_im(noise_msk))/nf nf]);
% noise_m = reshape(ym_im(noise_msk), [length(ym_im(noise_msk))/nf nf]);
% 
% % Covariance matrix estimate, of size (M x M); hopefully off-diagonals are small
% Sig_hat = cov([noise_p noise_m]);
% 
% % Compute the closest scaled-identity matrix to Sig_hat
% % This is simply given by scaling identity by the mean of the eigenvalues
% mean_var_hat = mean(eig(Sig_hat));      
% Sig_hat_diag = diag(repmat(mean_var_hat, [2*nf 1]));

%% PART ONE: Empirical vs. Theoretical comparisons
% 1a. Noise regions 
% Extract and vectorize noise
noise_msk = repmat(~imdilate(mask, strel('disk', 2)), [1 1 nf]);
noise_p = reshape(yp_sos(noise_msk), [length(yp_sos(noise_msk))/nf nf]);
noise_m = reshape(ym_sos(noise_msk), [length(ym_sos(noise_msk))/nf nf]);

% Noise mean empirical and theoretical estimates 
mu_hat  = mean(col([noise_p noise_m]));
mu_theo = (sqrt(2)*gamma((2*nc+1)/2)/gamma(nc)) * sigma_im; 

% Covariance matrix estimate, of size (M x M)
% Noise variance empirical and theoretical estimates
Sigma_hat = cov([noise_p noise_m]);
var_hat  = mean(eig(Sigma_hat));      
var_theo = (2*nc - (sqrt(2)*gamma((2*nc+1)/2)/gamma(nc))^2) * var_im; 

% % Compute the closest scaled-identity matrix to Sig_hat
% % This is simply given by scaling identity by the mean of the eigenvalues
% Sig_hat_diag = diag(repmat(var_hat, [2*nf 1]));
% figure; im(cat(3, Sig_hat, diag(diag(Sig_hat)), diag(diag(Sig_hat)) - Sig_hat,...
%     Sig_hat_diag, Sig_hat_diag - Sig_hat), 'cbar');

% Empirical vs. theoretical distributions
data = col([noise_p noise_m]); Xd = linspace(min(data), max(data), 100);...
    Yd = 2^(1-nc) * (Xd/sqrt(var_im)).^(2*nc-1) .* exp(-((Xd/sqrt(var_im)).^2)/2) / gamma(nc);...
    [Nd, Bd] = hist(data, 15);...
    figure; bar(Bd, Nd/max(Nd)*max(Yd)); hold on;...
    plot(Xd, Yd, 'r', 'LineWidth', 2); hold off; title('Sqrt-SoS Noise Distribution');

% 1b. Signal regions
% Extract and vectorize data
sig_mask = repmat((xx.^2 + yy.^2) <= (nx/16).^2, [1 1 nf]);
signal_p = reshape(yp_sos(sig_mask), [length(yp_sos(sig_mask))/nf nf]);
signal_m = reshape(ym_sos(sig_mask), [length(ym_sos(sig_mask))/nf nf]);

% Covariance matrix estimate, of size (M x M)
% In the high-SNR regime, var(Sqrt-SoS) \approx 0.5*var(single coil)
sig_Sigma_hat = cov([signal_p signal_m]);
sig_var_hat  = mean(diag(sig_Sigma_hat));
sig_var_theo = var_im;  

%% PART TWO: Four-Parameter VarPro for computing empirical (M0s/T1/T2/kappa) est. variance
% High-resolution dictionary creation
T1 = logspace(log10(500), log10(600), 50);
T2 = logspace(log10(38), log10(42), 50);
kappa = 2 .^ linspace(0.3, 0.5, 50);
D = NaN(2*length(flip), length(T1), length(T2), length(kappa));
for t1 = 1:length(T1)
    for t2 = 1:length(T2)
        for s = 1:length(kappa)
            for a = 1:nf
                [D(a, t1, t2, s), D(nf+a, t1, t2, s)] = ...
                    dess_fun(TR, TE, 1, T1(t1), T2(t2), flip(a)*kappa(s), wf);
            end
        end
    end
end
D = reshape(D, [2*length(flip) length(T1)*length(T2)*length(kappa)]);
D = abs(D);                             % Use magnitude dictionary for SoS

% (ML) Dictionary-based estimation via variable-projection method
weights = ones(2*nf, 1);
W = spdiags(weights, 0, 2*nf, 2*nf);    % Weighting Matrix
y = reshape(permute(cat(3, yp_sos, ym_sos), [3 1 2]), [2*nf nx*ny]);

maxProd = zeros(1, nx*ny);
idx = zeros(1, nx*ny);
for k = 1:(length(T1)*length(T2)*length(kappa))
    % Compute kth inner product
    hess = abs(D(:,k)' * W * D(:,k));
    ytild = D(:,k)' * W * y / sqrt(hess);
    newProd = abs(ytild).^2;

    % If the kth inner product is largest, save k
    update = newProd > maxProd;
    maxProd(update) = newProd(update);
    idx(update) = k;
end

% Extract indices for method-of-moments maps
[t1_idx, t2_idx, s_idx] = ind2sub([length(T1) length(T2) length(kappa)], idx);
T1_ml = reshape(T1(t1_idx), [nx ny]); 
T2_ml = reshape(T2(t2_idx), [nx ny]);
kappa_ml = reshape(kappa(s_idx), [nx ny]);

% M0s initial guess
M0s_ml = NaN(nx*ny, 1);
for q = 1:length(idx)
    M0s_ml(q) = (D(:,idx(q))' * y(:,q)) ./ (D(:,idx(q))' * D(:,idx(q)));
end
M0s_ml = reshape(M0s_ml, [nx ny]);

% % Mask out irrelevant voxels
% T1_ml(~mask) = 0; 
% T2_ml(~mask) = 0; 
% kappa_ml(~mask) = 0;
% M0s_ml(~mask) = 0;

% Empirical means and standard deviations from ROI
center = [32 32];
radius = 20;
roi_msk = false(size(mask));
roi_msk(center(1)-radius:center(1)+radius, center(2)-radius:center(2)+radius) = true; 
[M0s_ml_mean, M0s_ml_sd]      = multiMeans(M0s_ml, center, radius);
[T1_ml_mean, T1_ml_sd]        = multiMeans(T1_ml, center, radius);
[T2_ml_mean, T2_ml_sd]        = multiMeans(T2_ml, center, radius);
[kappa_ml_mean, kappa_ml_sd]  = multiMeans(kappa_ml, center, radius);

%% PART THREE: CRLB computation to obtain theoretical (M0s/T1/T2/kappa) est. variance
% Note: norm_crlb_dess_4parm() yields standard deviations for
% abs(M0s)*sqrt(SoS(sensitivities)) = 1. We must scale this CRLB 
% to match the signal level of the data. 

% Four-parameter time-uncompensated CRLB
TR_dess = TR * ones(nf, 1);
Sigma_inv = (1 ./ sig_var_theo) * speye(2*nf);
time_comp = 0;
[~, M0s_theory_sd, T1_theory_norm, T2_theory_norm, kappa_theory_norm]...
    = norm_crlb_dess_4parm(T1_ml_mean, T2_ml_mean, kappa_ml_mean, wf,...
    []', flip, []', TR_dess, TE, Sigma_inv, time_comp);

% Scale the normalized CRLB to match signal level of the data
% Note M0s_theory_sd is the std. dev. of M0s*sqrt(SoS(sensitivities)).
% M0s_theory_sd / sqrt(SoS(sens)) would be LATENT std. dev. of M0s
% M0s_theory_norm = M0s_theory_sd / M0s_ml_mean is std. dev. of M0s = 1.
T1_theory_sd  = T1_theory_norm  / M0s_ml_mean;
T2_theory_sd  = T2_theory_norm  / M0s_ml_mean;
kappa_theory_sd = kappa_theory_norm / M0s_ml_mean; 

% % Note: crlb_dess_kappa() incorporates the noise covariance automatically
% Sigma_inv = (1 ./ mean_var_hat) * speye(2*nf);
% TR_dess = TR * ones(nf, 1);
% scan_time = sum(TR_dess);
% 
% % Four-parameter time-compensated CRLB
% [~, M0s_theory_tcsd, T1_theory_tcsd, T2_theory_tcsd, kappa_theory_tcsd]...
%     = crlb_dess_kappa(M0s_ml_mean, T1_ml_mean, T2_ml_mean,...
%     kappa_ml_mean, wf, []', flip, []', TR_dess, TE, Sigma_inv);
% 
% % Undo time-compensation for comparison with empirical values
% M0s_theory_sd   = M0s_theory_tcsd / sqrt(scan_time);
% T1_theory_sd    = T1_theory_tcsd / sqrt(scan_time);
% T2_theory_sd    = T2_theory_tcsd / sqrt(scan_time);
% kappa_theory_sd = kappa_theory_tcsd / sqrt(scan_time);

%% PART FOUR: Graphical outputs
% Histogram of ML estimates to observe Gaussianity of estimator
Xm = linspace(min(M0s_ml(roi_msk)), max(M0s_ml(roi_msk)), 100);...
    Ym = normpdf(Xm, M0s_ml_mean, M0s_theory_sd);... 
    [Nm, Bm] = hist(M0s_ml(roi_msk), 20);... 
    figure; bar(Bm, Nm/max(Nm)*max(Ym)); hold on;...
    plot(Xm, Ym, 'r', 'LineWidth', 2); hold off;...
    title('Histogram of M0s Estimates');...
    legend('Empirical', 'Gaussian w/ CRLB std. dev.');
%     print -depsc 'M0s_hist_4parm,2flip.eps';
 
X1 = linspace(min(T1_ml(roi_msk)), max(T1_ml(roi_msk)), 100);...
    Y1 = normpdf(X1, T1_ml_mean, T1_theory_sd);...
    [N1, B1] = hist(T1_ml(roi_msk), 20); hold on;...
    figure; bar(B1, N1/max(N1)*max(Y1)); hold on;...
    plot(X1, Y1, 'r', 'LineWidth', 2); hold off;...
    title('Histogram of T1 Estimates');...
    legend('Empirical', 'Gaussian w/ CRLB std. dev.');
%     print -depsc 'T1_hist_4parm,2flip.eps';

X2 = linspace(min(T2_ml(roi_msk)), max(T2_ml(roi_msk)), 100);...
    Y2 = normpdf(X2, T2_ml_mean, T2_theory_sd);...
    [N2, B2] = hist(T2_ml(roi_msk), 20); hold on;...
    figure; bar(B2, N2/max(N2)*max(Y2)); hold on;...
    plot(X2, Y2, 'r', 'LineWidth', 2); hold off;...
    title('Histogram of T2 Estimates');...
    legend('Empirical', 'Gaussian w/ CRLB std. dev.');
%     print -depsc 'T2_hist_4parm,2flip.eps';

Xk = linspace(min(kappa_ml(roi_msk)), max(kappa_ml(roi_msk)), 100);...
    Yk = normpdf(Xk, kappa_ml_mean, kappa_theory_sd);...
    [Nk, Bk] = hist(kappa_ml(roi_msk), 20); hold on;...
    figure; bar(Bk, Nk/max(Nk)*max(Yk)); hold on;...
    plot(Xk, Yk, 'r', 'LineWidth', 2); hold off;...
    title('Histogarm of kappa Estimates');...
    legend('Empirical', 'Gaussian w/ CRLB std. dev.');
%     print -depsc 'kap_hist_4parm,2flip.eps';