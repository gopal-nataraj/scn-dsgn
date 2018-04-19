% script header_compare.m
% optimizing spgr scans for t1 estimation
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   1.1     2016-02-09      original
%   1.2     2016-07-18      add TR variation 

% constraints    
TR_tot = 1600;                        % ms    total TR time
TRs_res = 1;                          % ms    spgr TR resolution

flips_res = 1;                        % deg   spgr flip resolution
flips_min = 1;                        % deg   spgr flip minimum
flips_max = 120 - flips_min;          % deg   spgr flip maximum

% Test single-t1/kap vs. t1/kappa range
for i = 1:2
  switch i
    % wang et al
    case 1
      TRs_min = 800;                  % ms    min spgr TR  
      T1 = 1000;
      kap = 1;
      T1_rob = 1000;
      kap_rob = 1;
    % fleysher at al
    case 2
      TRs_min = 12.2;                 % ms    min spgr TR
      T1 = 1000;
      kap = 1;
      T1_rob = 1000;
      kap_rob = 1;
    % minmax, optimizing TR
    case 3
      TRs_min = 800;                  % ms    min spgr TR  
      T1 = linspace(800, 1400, 5);
      kap = linspace(0.9, 1.1, 5);
      T1_rob = linspace(400, 2000, 9);
      kap_rob = 2 .^ linspace(-1, 1, 5);
    % minmax, optimizing TR,flip
    case 4
      TRs_min = 12.2;                 % ms    min spgr TR  
      T1 = linspace(800, 1400, 5);
      kap = linspace(0.9, 1.1, 5);
      T1_rob = linspace(400, 2000, 9);
      kap_rob = 2 .^ linspace(-1, 1, 5);
  end

  T2 = 50;                            % typical constant, doesn't matter
  wf = 0;                             % rad/ms

  % Rescale measured noise variance to current voxel size
  oldVoxelVol = 1 * 1 * 5;                        % mm^3
  newVoxelVol = (240/256) * (240/256) * (30/6);   % mm^3
  oldvar_im_1coil = 2.62e-7;                      % a.u.; at 1x1x5 resolution 
  newvar_im_1coil = oldvar_im_1coil * (oldVoxelVol / newVoxelVol);

  % Other constant declarations
  TE = 4.67;                          % ms 
  noise_var_ssos = newvar_im_1coil/2; % High-SNR regime approximation

  a = 1;                              % T1 relative importance parameter
  b = 0;                              % T2 relative importance parameter

  tol = 0.01;                         % Global minima tolerance
  time_comp = 0;                      % Toggle time compensation on/off

  % Choose a scan profile: (2,0)
  n_spgr = 2;                         % Number of SPGR scans
  n_dess = 0;                         % Number of DESS scans
  profile = sprintf('%u SPGR, %u DESS', n_spgr, n_dess);

  % Introduce the experiment and start timer.
  fprintf('\n\nScan Profile: (%u SPGR, %u DESS)\n', n_spgr, n_dess);
  fprintf('Total TR constraint: %1.2f ms\n', TR_tot);
  fprintf('Weightings: (%1.1f,%1.1f,%1.1f)\n', 0, a, b);
  tic; 

  %% Scan sequence: 2 SPGR, 0 DESS
  % Controllable parameter grid: vary TRs, flips
  TRs_max = TR_tot - TRs_min;
%   TR1_spgr = [TRs_min : TRs_res : TRs_max];
  if TRs_min == TRs_max
    TR1_spgr = TRs_min;
  else
    TR1_spgr = logspace(log10(TRs_min), log10(TRs_max), 50);
  end
  spgr1 = [flips_min : flips_res : flips_max] * pi/180;     % rad
  spgr2 = [flips_min : flips_res : flips_max] * pi/180;     % rad

  % Grid definition
  sig_T1_20   = inf(length(TR1_spgr), length(spgr1), length(spgr2));
  Psi_20      = inf(length(TR1_spgr), length(spgr1), length(spgr2));

  % Covariance matrix: scaled identity with M = 2
  Sigma_inv = (1 ./ noise_var_ssos) * speye(2);

  %% Step one: min-max crlb analysis over parameter grid
  dess_20 = []';
  TRd_20 = []';
  for i_s1 = 1:length(spgr1)
    for i_s2 = 1:length(spgr2)
      spgr_20 = [spgr1(i_s1) spgr2(i_s2)]';

      for i_TRs1 = 1:length(TR1_spgr)
        % Second SPGR TR is already specified from first + constraint
        TRs_20 = [TR1_spgr(i_TRs1) TR_tot-TR1_spgr(i_TRs1)]';

        % Inner maximization: store max SDs over tight (T1,kap) range
        worst_sigT1 = NaN(length(T1), length(kap));

        for i_t1 = 1:length(T1)
          for i_kap = 1:length(kap)
            [~, ~, worst_sigT1(i_t1, i_kap)] = ...
              norm_crlb_dess_3parm(T1(i_t1), T2, wf,...
              kap(i_kap)*spgr_20, kap(i_kap)*dess_20,...
              TRs_20, TRd_20, TE, Sigma_inv, time_comp);
          end
        end

        % Store the worst-case sig_T1 value
        sig_T1_20(i_TRs1, i_s1, i_s2) = max(worst_sigT1(:));
        Psi_20(i_TRs1, i_s1, i_s2) = max(worst_sigT1(:));
      end
    end
  end

  % Find the indices of (multiple) Psi_20 minima, to within a tolerance
  Psimin_idx_20 = find( (Psi_20 - min(Psi_20(:))) ./ min(Psi_20(:)) <= tol);
  num_min_20 = length(Psimin_idx_20);

  %% Step two: select one minimum based on robustness over broad (T1, kap) range
  % Indices of each of the "num_min_20" total minima
  TRs1_idx_20 = NaN(num_min_20, 1);
  s1_idx_20 = NaN(num_min_20, 1);
  s2_idx_20 = NaN(num_min_20, 1);

  sigT1_worst_20  = NaN(num_min_20, length(T1_rob), length(kap_rob));

  for i_min = 1:num_min_20
    % Convert the 1D index to ND-grid indices
    [TRs1_idx_20(i_min), s1_idx_20(i_min), s2_idx_20(i_min)]...
      = ind2sub(size(Psi_20), Psimin_idx_20(i_min));

    % Store the minimizing parameters
    TRs1min_20  = TR1_spgr(TRs1_idx_20(i_min));
    TRs2min_20  = TR_tot - TRs1min_20;
    s1min_20    = spgr1(s1_idx_20(i_min));
    s2min_20    = spgr2(s2_idx_20(i_min));

    % Evalue CRLB at minimizers, over (T1_rob, kap_rob) ranges
    for i_t1 = 1:length(T1_rob)
      for i_kap = 1:length(kap_rob)
        [~, ~, sigT1_worst_20(i_min, i_t1, i_kap), ~] = ...
          norm_crlb_dess_3parm(T1_rob(i_t1), T2, wf,...
          kap_rob(i_kap)*[s1min_20 s2min_20]', kap_rob(i_kap)*[]',...
          [TRs1min_20 TRs2min_20]', []', TE, Sigma_inv, time_comp);
      end
    end
  end

  % Store the worst-case sigT1
  sigT1_min_20 = max(reshape(sigT1_worst_20,...
    [num_min_20 length(T1_rob)*length(kap_rob)]), [], 2);

  % Compute the worst-case Psi_20 value over the wider range
  Psi_min_20 = sigT1_worst_20;
  [Psi_maxmin_20, wide_idx_20] = max(reshape(Psi_min_20,...
    [num_min_20 length(T1_rob)*length(kap_rob)]), [], 2);

  % Extract the parameter indices that minimize Psi_maxmin_20 over minima
  % Psi_star - min(Psi(:)) is the degradation from narrow->wide range
  [Psi_star_20, i_star_20] = min(Psi_maxmin_20);
  sigT1_star_20 = sigT1_min_20(i_star_20);

  TRs1star_20     = TRs1_idx_20(i_star_20);
  s1star_20       = s1_idx_20(i_star_20);
  s2star_20       = s2_idx_20(i_star_20);

  % Display results 
  fprintf('Selected scan design: (%1.1f ms, %1.1f deg); (%1.1f ms, %1.1f deg)\n',...
    TR1_spgr(TRs1star_20), spgr1(s1star_20)*180/pi,...
    TR_tot-TR1_spgr(TRs1star_20), spgr2(s2star_20)*180/pi);
  fprintf('\t Worst-case sigT1: %6.2f\n',...
    sig_T1_20(TRs1star_20, s1star_20, s2star_20));
  fprintf('\t Robust-range sigT1: %6.2f\n', sigT1_star_20);
  fprintf('\t Worst-case Psi: (%6.1f * %6.2f) + (%6.1f * %6.3f) = %6.3f\n',...
    a, sig_T1_20(TRs1star_20, s1star_20, s2star_20),...
    b, 0,...
    Psi_20(TRs1star_20, s1star_20, s2star_20));
  fprintf('\t Robust-range Psi: %6.3f\n', Psi_star_20);
  fprintf('\t Selected from %d candidates, at a %d%% tolerance, in %5.0f seconds.\n\n',...
    num_min_20, 100*tol, toc);

  %% Graphical Output: 2D Heatmap for 2 SPGR, 0 DESS
  % Minimizer location(s)
  if (TR1_spgr(TRs1star_20) == TR_tot/2)
    xmin = [spgr1(s1star_20) spgr2(s2star_20)] * 180/pi;
    ymin = [spgr2(s2star_20) spgr1(s1star_20)] * 180/pi;
  else
    xmin = [spgr1(s1star_20)] * 180/pi;
    ymin = [spgr2(s2star_20)] * 180/pi;
  end

  % save 
  sv = 1;
  if sv
    switch i
      case 1
        save('sigt1-vs-s1s2,20,wang.mat',...
          'sig_T1_20', 'spgr*', 'TRs1star_20', 'xmin', 'ymin', 'i');
      case 2
        save('sigt1-vs-s1s2,20,fleysher.mat',...
          'sig_T1_20', 'spgr*', 'TRs1star_20', 'xmin', 'ymin', 'i');
      case 3
        save('sigt1-vs-s1s2,20,minmax,flip.mat',...
          'sig_T1_20', 'spgr*', 'TRs1star_20', 'xmin', 'ymin', 'i');
      case 4
        save('sigt1-vs-s1s2,20,minmax,flip-tr.mat',...
          'sig_T1_20', 'spgr*', 'TRs1star_20', 'xmin', 'ymin', 'i');
    end
  end
  
  % run from here down for plotting from saved heatmaps
  T1range = [3 7];                    % sigmaT1 plot range
  pr = 1;                             % Toggle print on/off
  
  switch i
    case 1
      tit = 'T1 Std. Dev. for nominal (T1,kap), w/o TR optimization';
      out = 'sigt1-vs-s1s2,20,wang.eps';
    case 2
      tit = 'T1 Std. Dev. for nominal (T1,kap), w/ TR optimization';
      out = 'sigt1-vs-s1s2,20,fleysher.eps';
    case 3
      tit = 'Max T1 Std. Dev. for GM/WM ROIs w/o TR optimization';
      out = 'sigwt1-vs-s1s2,20,minmax,flip.eps';
    case 4
      tit = 'Max T1 Std. Dev. for GM/WM ROIs, w/ TR optimization';
      out = 'sigwt1-vs-s1s2,20,minmax,flip-tr.eps';
  end
  
  tit = '';
  figure; 
  hold on; 
  imagesc(180*spgr1/pi, 180*spgr2/pi, ...
    squeeze(sig_T1_20(TRs1star_20, :, :))', T1range);...
  scatter(xmin, ymin, 250, 'p', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w');
  hold off;
  
  xlabel('SPGR \alpha_1 (deg)', 'FontSize', 14);
  ylabel('SPGR \alpha_2 (deg)', 'FontSize', 14);
  title(tit, 'FontSize', 14);
  axis xy square; 
  axis([0 120 0 120]);
  colormap('jet'); 
  colorbar;
  set(gca, 'Xtick', [0:10:120]);
  set(gca, 'Ytick', [0:10:120]);
  
  if pr
    print('-depsc', strcat(out));
  end
end

% % extra code to evaluate sigw at certain points
% % spgr_20 = [31 100]' * pi/180;
% spgr_20 = [24 102]' * pi/180;
% dess_20 = [];
% % TRs_20 = [1012.1 587.9]';
% TRs_20 = [874.8 725.2]';
% TRd_20 = [];
% TE = 4.67;                          % ms 
% 
% % T1 = linspace(800, 1400, 5);
% % kap = linspace(0.9, 1.1, 5);
% T1 = 1000;
% kap = 1;
% T2 = 50;                            % typical constant, doesn't matter
% wf = 0;                             % rad/ms
% 
% oldVoxelVol = 1 * 1 * 5;                        % mm^3
% newVoxelVol = (240/256) * (240/256) * (30/6);   % mm^3
% oldvar_im_1coil = 2.62e-7;                      % a.u.; at 1x1x5 resolution 
% newvar_im_1coil = oldvar_im_1coil * (oldVoxelVol / newVoxelVol);
% noise_var_ssos = newvar_im_1coil/2; % High-SNR regime approximation
% worst_sigT1 = NaN(length(T1), length(kap));
% Sigma_inv = (1 ./ noise_var_ssos) * speye(2);
% 
% time_comp = 0;                      % Toggle time compensation on/off
% 
% for i_t1 = 1:length(T1)
%   for i_kap = 1:length(kap)
%     [~, ~, worst_sigT1(i_t1, i_kap)] = ...
%       norm_crlb_dess_3parm(T1(i_t1), T2, wf,...
%       kap(i_kap)*spgr_20, kap(i_kap)*dess_20,...
%       TRs_20, TRd_20, TE, Sigma_inv, time_comp);
%   end
% end
% 
% max(worst_sigT1(:))

