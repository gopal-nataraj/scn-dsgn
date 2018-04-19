% script 2015_06_26_test.m
% test to check joint irse-se recon with old data
%
% gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2015-06-26    data acquired
%   v1.x    2015-07-28    initial processing of paper data
%   v2.1    2016-06-22    checking new model/workflow

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd;
  cd ../../../../irt;
  irtdir = pwd;
  setup();
  cd(curdir);
end

% add relevant directories
addpath('../../data/Brain_06,26,15');
addpath('../../data/Brain_06,26,15/BS');
addpath('../../data/Brain_06,26,15/IR');
addpath('../../data/Brain_06,26,15/SE');
addpath('../../map/t1-t2');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../etc');

% constant declarations
n.x = 256;
n.y = 256;
n.z = 6;
n.c = 32;
n.sl = 2;
n.p = 3;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

% acquisition parameters
P.ir.ti   = [650 900]';                 % ms
C.ir      = length(P.ir.ti);
E.ir      = 1;
P.ir.tr   = ones(C.ir,1) * 3000;        % ms
P.ir.te   = ones(C.ir,E.ir) * 14;       % ms
P.ir.aex  = ones(C.ir,1) * pi/2;        % rad
P.ir.aref = ones(C.ir,1) * pi;          % rad

P.se.te   = [10 150]';                  % ms
C.se      = length(P.se.te);
E.se      = 1;         
P.se.tr   = ones(C.se,1) * 3000;        % ms
P.se.aex  = ones(C.se,1) * pi/2;        % rad
P.se.aref = ones(C.se,1) * pi;          % rad

% header options
bool.sc   = true;
bool.im   = false;
bool.pr   = true;
bool.rs   = true;

% mask threshold
mask.thresh = 0.05;

% recon options
dict.t1.n = 200;                        % num t1 dict elements
dict.t2.n = 200;                        % num t2 dict elements
stop.iter = 10;                        % iteration stop criterion
bool.chat = true;                       % verbosity
bool.disp = false;                      % show image updates
bool.precon = true;                     % use preconditioner for rls estimation

% load nu
load('B1_reg');
if bool.sc
  scale = 1.016/0.913;
else
  scale = 1;
end
nu.b1.mom = b1_init * scale;
nu.b1.rls = b1_reg * scale;
nu.kap.rls = kap_reg * scale;
clear('b1_*', 'kap_reg');

% load coil-combined image data
load('ims_coil_comb_IR.mat');
ycc.ir = reshape(y_im, [n.x n.y C.ir E.ir]);
load('ims_coil_comb_SE.mat');
ycc.se = reshape(y_im, [n.x n.y C.se E.se]);
clear('y_im', 'smap');

% create masks from short te image
tmp = squeeze(abs(ycc.se(:,:,1)));
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

% test different scan profiles
for p = 1:n.p
  % set data weights based on profile
  switch p
    case 1                              % (4 ir, 4 se)
      wght.ir = [1 1 1 1]'  * ones(1,E.ir);     
      wght.se = [1 1 1 1]'  * ones(1,E.se);
    case 2                              % (4 ir, 4 se)
      wght.ir = [1 1 1 1]'  * ones(1,E.ir);     
      wght.se = [0 0 0 0]'  * ones(1,E.se);  
    case 3                              % (4 ir, 4 se)
      wght.ir = [0 0 0 0]'  * ones(1,E.ir);     
      wght.se = [1 1 1 1]'  * ones(1,E.se);
  end
  fprintf('\nTesting profile %d: (%u,%u) ir/se scans.\n',...
    p, sum(wght.ir(:,1)), sum(wght.se(:,1)));
  
  % parameter estimation
  mapArg = {...
    'mask.disp', mask.t,...
    'mask.est', mask.b,...
    'nu.kap', double(nu.kap.rls),...
    'wght.ir', wght.ir,...
    'wght.se', wght.se,...
    'dict.t1.n', dict.t1.n,...
    'dict.t2.n', dict.t2.n,...
    'stop.iter', stop.iter,...
    'bool.chat', bool.chat,...
    'bool.disp', bool.disp,...
    'bool.precon', bool.precon};
  [x(p), t(p)] = mri_m0t1t2_map(ycc, P, mapArg{:});
end

% roi selection
if bool.rs
  tmp = 'ROI selection: [d]efauls, [m]anual, or [s]aved? ';
  roi.sel = input(tmp, 's');
else
  roi.sel = 's';
end
roi.label = {'wm'; 'gm'; 'sf'};
n.roi = length(roi.label);

switch roi.sel
  case 'd'
    % define circular rois
    roi.x = [147 180 120];
    roi.y = [ 53 164  88];
    roi.r = [  4   2   2];
    
    roi.mask = false(n.x, n.y, n.roi);
    [roi.grid.x, roi.grid.y] = ndgrid(1:n.x, 1:n.y);
    for r = 1:n.roi
      tmp = false(n.x, n.y);
      tmp((roi.grid.x-roi.x(r)).^2 + (roi.grid.y-roi.y(r)).^2 <= roi.r(r)^2) = true;
      roi.mask(:,:,r) = tmp;
    end
    
  case 'm'
    roi.mask = false(n.x, n.y, n.roi);
    for r = 1:n.roi
      while true
        fprintf('Select %s polygon.\n', roi.label{r});
        tmp = roipoly(x(3).t2.rls/200);
        roi.mask(:,:,r) = roi.mask(:,:,r) | tmp;
        tmp = input('Add another polygon to this ROI [*/n]? ', 's');
        if strcmp(tmp,'n'), break; end;
      end
    end

    % save mask
    dir = pwd; 
    cd('../../data/Brain_06,26,15');
    if exist('roi_masks.mat', 'file')
      tmp = input('Overwrite previously-defined ROI mask [y/n]? ', 's');
      switch tmp
        case 'y'
          save('roi_masks.mat', 'roi'); 
          fprintf('roi_masks.mat overwritten.\n');
        case 'n'
          fprintf('New ROI mask not saved.\n');
        otherwise
          error('Unrecognized input.');
      end
    else
      save('roi_masks.mat', 'roi');
    end
    cd(dir);

  case 's'
    try
      addpath('../../data/Brain_06,26,15');
      load('roi_masks.mat');
      roi.mask = roi_mask; 
      clear('roi_mask');
    catch
      error('ROI masks not found!');
    end

  otherwise
    error('Unrecognized input!');
end

% boundaries
roi.bound.wm = bwboundaries(roi.mask(:,:,1), 'noholes');
roi.bound.gm = bwboundaries(roi.mask(:,:,2), 'noholes');

% images
rng.m0 = [0 n.c];
rng.t1 = [500 2000];
rng.t2 = [20 120];
if bool.im
  for p = 1:n.p
    figure; im('notick', x(p).t1.ml, rng.t1, 'cbar', ' ');...
      if bool.roi
        hold on;
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1), tmp(:,2), 'm', 'LineWidth', 1);
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1), tmp(:,2), 'c', 'LineWidth', 1);
        end
        hold off;
      end
      if bool.pr
        print('-depsc', sprintf('t1,ml,prof-%u.eps', p));
      end
      
    figure; im('notick', x(p).t1.rls, rng.t1, 'cbar', ' ');...
      if bool.roi
        hold on;
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1), tmp(:,2), 'm', 'LineWidth', 1);
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1), tmp(:,2), 'c', 'LineWidth', 1);
        end
        hold off;
      end
      if bool.pr
        print('-depsc', sprintf('t1,rls,prof-%u.eps', p));
      end
      
    figure; im('notick', x(p).t2.ml, rng.t2, 'cbar', ' ');...
      if bool.roi
        hold on;
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1), tmp(:,2), 'm', 'LineWidth', 1);
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1), tmp(:,2), 'c', 'LineWidth', 1);
        end
        hold off;
      end
      if bool.pr
        print('-depsc', sprintf('t2,ml,prof-%u.eps', p));
      end
      
    figure; im('notick', x(p).t2.rls, rng.t2, 'cbar', ' ');...
      if bool.roi
        hold on;
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1), tmp(:,2), 'm', 'LineWidth', 1);
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1), tmp(:,2), 'c', 'LineWidth', 1);
        end
        hold off;
      end
      if bool.pr
        print('-depsc', sprintf('t2,rls,prof-%u.eps', p));
      end
  end
else
  tmp = [x(:).m0];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.m0, 'cbar');
  tmp = [x(:).t1];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.t1, 'cbar');
  tmp = [x(:).t2];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.t2, 'cbar');
end

% summary statistics
stat.t1.ml.mean   = NaN(n.p, n.roi);
stat.t1.ml.std    = NaN(n.p, n.roi);
stat.t1.rls.mean  = NaN(n.p, n.roi);
stat.t1.rls.std   = NaN(n.p, n.roi);
stat.t2.ml.mean   = NaN(n.p, n.roi);
stat.t2.ml.std    = NaN(n.p, n.roi);
stat.t2.rls.mean  = NaN(n.p, n.roi);
stat.t2.rls.std   = NaN(n.p, n.roi);
for p = 1:n.p
  for r = 1:n.roi
    stat.t1.ml.mean(p,r)    = mean(masker(x(p).t1.ml, roi.mask(:,:,r)));
    stat.t1.ml.std(p,r)     = std(masker(x(p).t1.ml, roi.mask(:,:,r)));
    stat.t1.rls.mean(p,r)   = mean(masker(x(p).t1.rls, roi.mask(:,:,r)));
    stat.t1.rls.std(p,r)    = std(masker(x(p).t1.rls, roi.mask(:,:,r)));
    stat.t2.ml.mean(p,r)    = mean(masker(x(p).t2.ml, roi.mask(:,:,r)));
    stat.t2.ml.std(p,r)     = std(masker(x(p).t2.ml, roi.mask(:,:,r)));
    stat.t2.rls.mean(p,r)   = mean(masker(x(p).t2.rls, roi.mask(:,:,r)));
    stat.t2.rls.std(p,r)    = std(masker(x(p).t2.rls, roi.mask(:,:,r)));
  end
end

% print
if bool.pr
  fid = fopen('../../data/Brain_06,26,15/2015-06-26,stat', 'w');
  fprintf(fid, '2015-06-26 parameter estimate statistics\n');
  fprintf(fid, '\tcolumns are profiles\n');
  fprintf(fid, '\trows are ROIs\n');
  
  fprintf(fid, '\nT1 ML:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.ml.mean(p,r), char(177), stat.t1.ml.std(p,r));
    end
    fprintf(fid, '\n');
  end
  
  fprintf(fid, '\nT1 RLS:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.rls.mean(p,r), char(177), stat.t1.rls.std(p,r));
    end
    fprintf(fid, '\n');
  end
  
  fprintf(fid, '\nT2 ML:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.ml.mean(p,r), char(177), stat.t2.ml.std(p,r));
    end
    fprintf(fid, '\n');
  end
  
  fprintf(fid, '\nT2 RLS:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.rls.mean(p,r), char(177), stat.t2.rls.std(p,r));
    end
    fprintf(fid, '\n');
  end
  fclose(fid);
end