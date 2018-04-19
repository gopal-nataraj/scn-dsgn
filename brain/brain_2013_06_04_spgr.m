%% Header Script for Regularized M0/T1 Recon from 06/05/13 Data

% Load Data
cd InVivo060413; loaddata_invivo(); cd ..;

% Set parameters
TR = 22.1;                              % ms
TE = 5.7;                               % ms
Kbs = 6.3638e+03;                       % deg/Gauss^2
FOVx = 260; FOVy = 260; FOVz = 60;      % mm
flip_spgr = [10 20 30 60]' * pi/180;    % rad
flip_dess = [20 40 60]' * pi/180;       % rad

% Compute sum of squares
im_spgr_sos = squeeze(sqrt(sum(abs(im_spgr).^2, 4)));
[nx, ny, nz, nf] = size(im_spgr_sos);
dx = FOVx/nx; dy = FOVy/ny; dz = FOVz/nz;   % mm
voxelVol = dx*dy*dz; 

% Regularized B1 map estimation for flip maps
b1_init = col(b1map);
bs_data = cat(5, implus, imminus); bs_size = size(bs_data);
bs_data = permute(reshape(bs_data, [prod(bs_size(1:3)) ...
    bs_size(4) bs_size(5)]), [1 3 2]);
niter_bs = 500;
delta_bs = 0.001; 
beta_bs = 2^25; 
bs_mask = true(bs_size(1:3));
R_bs = Reg1(bs_mask, 'pot_arg', {'quad'}, 'beta', beta_bs);
[b1reg, b1gam] = regdbsb1(b1_init, bs_data, Kbs, R_bs, niter_bs);
b1reg = embed(b1reg(:,niter_bs), bs_mask);

% Expand scaleflip map to xdim x ydim x zdim
b1_nom = 0.05; scaleflip = b1reg/b1_nom; 
scaleflip_expanded = permute(imresize(permute(imresize(scaleflip,...
    [nx ny]), [2 3 1]), [ny nz]), [3 1 2]);

% Flip angle correction
flip_spgr_map = NaN(size(im_spgr_sos));
for ang = 1:nf
    flip_spgr_map(:,:,:,ang) = scaleflip_expanded * flip_spgr(ang);
end

%% MoM T1 MAP ESTIMATION 

% Linear Least-Squares Regression (DESPOT-1)
% For W, error covariance approximates E1 to 1
if (~exist('T1_map_mom','var'))
    try 
        S = load('T1map_wmom_invivo_trueflip');
%         E1_map_mom = S.E1_map_mom;
        T1_map_mom = S.T1_map_mom;
        M0_map_mom = S.M0_map_mom;
    catch
        E1_map_mom = NaN(nx, ny, nz);
        intercept_map_mom = NaN(nx, ny, nz);
        
        for xx = 1:nx
            for yy = 1:ny
                for zz = 1:nz
                    if (scaleflip_expanded(xx,yy,zz) > 0.03)
                        Y = squeeze(im_spgr_sos(xx,yy,zz,:) ./ ...
                            sin(flip_spgr_map(xx,yy,zz,:)));
                        X = squeeze(im_spgr_sos(xx,yy,zz,:) ./ ...
                            tan(flip_spgr_map(xx,yy,zz,:)));
                        W = diag(squeeze((1-cos(flip_spgr_map(xx,yy,zz,:))) ...
                            ./ sin(flip_spgr_map(xx,yy,zz,:)))); 
                        A = [X ones(length(flip_spgr), 1)];
                        regression = (A' * W * A) \ (A' * W * Y);
                        E1_map_mom(xx,yy,zz) = regression(1);
                        intercept_map_mom(xx,yy,zz) = regression(2);
                    end
                end
            end
        end
        
        T1_map_mom = -TR / log(abs(E1_map_mom));     % Should be pure real
        M0_map_mom = intercept_map_mom ./ (1-E1_map_mom);
        save('T1map_wmom_invivo_trueflip', 'T1_map_mom', 'E1_map_mom', 'M0_map_mom');
    end
end 

% Project Images to within range 
T1_max = 3000; 
E1_max = exp(-TR/T1_max);
M0_max = 5;

T1_map_mom = max(T1_map_mom, 0);
T1_map_mom = min(T1_map_mom, T1_max);
M0_map_mom = max(M0_map_mom, 0);
M0_map_mom = min(M0_map_mom, M0_max);

%% PREPROCESSING AND MASKING

% Masking
tight_mask = (M0_map_mom > 0.1) & (T1_map_mom < T1_max) & ...
    (T1_map_mom > 0) & ~isnan(T1_map_mom); 
loose_mask = imdilate(tight_mask, strel('disk', 15));

% Set voxels inside loose mask but outside tight mask to mean
M0_map_mom(~tight_mask & loose_mask) = mean(col(M0_map_mom));
T1_map_mom(~tight_mask & loose_mask) = mean(col(T1_map_mom));

M0_map_mom(~loose_mask) = 0;
T1_map_mom(~loose_mask) = 0;

% Compute Spin Density from M0
rho_map_mom = M0_map_mom / voxelVol;             % T/mm^3

% Update the E1 map
E1_map_mom = exp(-TR / T1_map_mom);

%% ITERATIVE RECONSTRUCTION

% Define iteration parameters
slice = 6;
n_outer_iter = 5;
niter1 = 100;
niter2 = 100; 
tol_1 = 10^-6;
tol_2 = 10^-6;

% Define regularizers, R1 and R2
beta1 = 2^-10;
beta2 = 2^-2;
delta1 = 10^-3;
delta2 = 10^-3;
R1 = Reg1(loose_mask(:,:,slice), 'pot_arg', {'hyper3', delta1}, 'beta', beta1);
R2 = Reg1(loose_mask(:,:,slice), 'pot_arg', {'hyper3', delta2}, 'beta', beta2);

% Regularized reconstruction
[M0map, T1map, cost] = mri_despot1_m0t1map(M0_map_mom(:,:,slice), ...
    E1_map_mom(:,:,slice), im_spgr_sos(:,:,slice,:), ...
    flip_spgr_map(:,:,slice,:), loose_mask(:,:,slice), R1, R2, ...
    T1_max, TR, n_outer_iter, niter1, niter2, tol_1, tol_2);
figure; im(M0map, [0 3], 'cbar');
figure; im(T1map, [0 1500], 'cbar');

% Plotting cost vs. iteration
figure; hold on;
scatter([1:2:2*n_outer_iter], cost(2:2:end), 'bo');
scatter([2:2:2*n_outer_iter], cost(3:2:end), 'ro');
plot([0:2*n_outer_iter], cost, 'k'); hold off;
title('Cost vs. iteration');
legend('After x1 (M0) update','After x2 (E1) update');

%% SANITY CHECKS AND SIGNAL COMPARISONS

% Test parameters
a_test = [0:0.1:90]'*pi/180;
radius = 2;
whm1 = [75 100];
csf1 = [93 133];

E1map = exp(-TR ./ T1map);
[whm1_mom_mod, whm1_dat, whm1_flip] = signal_model(squeeze(M0_map_mom(:,:,slice)),...
    squeeze(E1_map_mom(:,:,slice)), squeeze(im_spgr_sos(:,:,slice,:)), ...
    whm1, radius, a_test, squeeze(flip_spgr_map(:,:,slice,:))); 
[whm1_npls_mod, ~, ~] = signal_model(M0map, E1map, ...
    squeeze(im_spgr_sos(:,:,slice,:)), whm1, radius, a_test,...
    squeeze(flip_spgr_map(:,:,slice,:))); 

[csf1_mom_mod, csf1_dat, csf1_flip] = signal_model(squeeze(M0_map_mom(:,:,slice)),...
    squeeze(E1_map_mom(:,:,slice)), squeeze(im_spgr_sos(:,:,slice,:)),...
    csf1, radius, a_test, squeeze(flip_spgr_map(:,:,slice,:))); 
[csf1_npls_mod, ~, ~] = signal_model(M0map, E1map, ...
    squeeze(im_spgr_sos(:,:,slice,:)), csf1, radius, a_test,...
    squeeze(flip_spgr_map(:,:,slice,:))); 

% Raw signal vs. flip
figure; hold on;
plot(a_test*180/pi, whm1_mom_mod, 'b');
plot(a_test*180/pi, whm1_npls_mod, 'r');
errorbar(whm1_flip*180/pi, mean(whm1_dat, 1), std(whm1_dat, 1), 'k*');
hold off; grid on;
xlabel('Flip angle'); ylabel('SPGR signal'); title('WM');
legend('Method of Moments','Regularized','Raw data');

figure; hold on;
plot(a_test*180/pi, csf1_mom_mod, 'b');
plot(a_test*180/pi, csf1_npls_mod, 'r');
errorbar(csf1_flip*180/pi, mean(csf1_dat, 1), std(csf1_dat, 1), 'k*');
hold off; grid on;
xlabel('Flip angle'); ylabel('SPGR signal'); title('CSF');
legend('Method of Moments','Regularized','Raw data');

