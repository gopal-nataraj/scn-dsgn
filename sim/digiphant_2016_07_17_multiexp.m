% script digiphant_2016_02_22_multiexp_test.m
% simulations to find whether model mismatch due to multiexponential T2
%
% Written by: Gopal Nataraj
% Copyright 2016, University of Michigan
%
% Version control
%   v1.1      2016-02-22      spin echo modeling only
%   v1.2      2016-07-16      minor tweaks for paper

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/DigitalPhantom');
addpath('../../map/t1-t2');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../etc');

% constant declarations
c.b0 = 3.0;
c.sl = 81;
c.sig.ir = 0;
c.sig.se = 0;
c.rng.inveff.tru = 0.00;
c.rng.kap.tru = 0.0;
c.rng.kap.est = 0.0;
c.rng.b0.tru  = 0.00;                   % kHz
c.rng.b0.est  = 0.00;                   % kHz
c.rng.r2p.tru = 0.00;                   % kHz
c.rng.r2p.est = 0.00;                   % kHz

% header options
bool.im   = true;                       % show images
bool.pr   = true;                       % print images and statistics

% compartmental water weights
%   rows are [my; ie; fr]
%   cols are [sf, gm, wm]
comp.wght = [0.00 0.00 0.15;...
             0.00 0.95 0.80;...
             1.00 0.05 0.05];
[n.comp, n.roi] = size(comp.wght);

% nominal compartmental values (my, ie, fr)
comp.xnom.t1  = [ 500 1000 3500];       % ms
comp.xnom.t2  = [  20   80  250];       % ms
comp.xnom.t2s = [  10   40   80];       % ms

% load digital phantom file
f.name = 'phantom_1.0mm_msles2_crisp.fld';
f.label = flipdim(fld_read(f.name, 'slice', c.sl), 2);
[n.x, n.y] = size(f.label);
n.p = 4;

% true multi-component maps (no exchange)
% assume each component perfectly inverted for simplicity
xtrue = struct('m0', [], 't1', [], 't2', [], 't2s', [], 'inveff', []);
[xx, yy] = ndgrid(linspace(-1,1,n.x), linspace(-1,1,n.y));
for i = 1:n.comp
  xtrue(i).m0 = NaN(n.x, n.y);
  xtrue(i).t1 = NaN(n.x, n.y);
  xtrue(i).t2 = NaN(n.x, n.y);
  xtrue(i).t2s = NaN(n.x, n.y);
  xtrue(i).inveff = (1.0+c.rng.inveff.tru/2) - c.rng.inveff.tru*(xx.^2 + yy.^2);
  for r = 0:10
    tmp = mri_brainweb_params(r, 'b0', c.b0);
    switch r
      case {1,2,3}                      % csf, gm, wm 
        xtrue(i).m0(f.label == r) = comp.wght(i,r) * tmp.pd;
        xtrue(i).t1(f.label == r) = comp.xnom.t1(i);
        xtrue(i).t2(f.label == r) = comp.xnom.t2(i);
        xtrue(i).t2s(f.label == r) = comp.xnom.t2s(i);
      otherwise   
        xtrue(i).m0(f.label == r) = tmp.pd/n.comp;
        xtrue(i).t1(f.label == r) = tmp.t1;
        xtrue(i).t2(f.label == r) = tmp.t2;
        xtrue(i).t2s(f.label == r) = tmp.t2s;
    end
  end
end

% make known parameter maps
nu.tru.kap = (1+c.rng.kap.tru/2) - c.rng.kap.tru*(xx.^2 + yy.^2);
nu.tru.b0  = (0+c.rng.b0.tru/2)  - c.rng.b0.tru *(xx.^2 + yy.^2);
nu.tru.r2p = (0+c.rng.r2p.tru/2) - c.rng.r2p.tru*(xx.^2 + yy.^2);
nu.est.kap = (1+c.rng.kap.est/2) - c.rng.kap.est*(xx.^2 + yy.^2);
nu.est.b0  = (0+c.rng.b0.est/2)  - c.rng.b0.est *(xx.^2 + yy.^2);
nu.est.r2p = (0+c.rng.r2p.est/2) - c.rng.r2p.est*(xx.^2 + yy.^2);

% acquisition parameters
P.ir.ti   = [50 150 450 1350]';         % ms
C.ir      = length(P.ir.ti);
E.ir      = 1;
P.ir.tr   = ones(C.ir,1) * 1400;        % ms
P.ir.te   = ones(C.ir,E.ir) * 10;       % ms
P.ir.ainv = ones(C.ir,1) * pi;          % rad
P.ir.aex  = ones(C.ir,1) * pi/2;        % rad
P.ir.aref = ones(C.ir,1) * pi;          % rad

P.se.te   = [10 30 60 150]';            % ms
C.se      = length(P.se.te);
E.se      = 1;         
P.se.tr   = ones(C.se,1) * 1000;        % ms
P.se.aex  = ones(C.se,1) * pi/2;        % rad
P.se.aref = ones(C.se,1) * pi;          % rad

% recon options
stop.iter = 0;                          % iteration stop criterion
bool.mag.ir = false;                    % use magnitude se-ir data/model
bool.mag.se = false;                    % use magnitude se data/model
bool.chat = true;                       % verbosity
bool.norm = false;                      % normalize coil-combined data before recon
bool.disp = false;                      % show image updates
bool.precon = true;                     % use preconditioner for rls estimation

% collect noiseless ir data
y.ir = zeros(n.x, n.y, C.ir, E.ir);
for d = 1:C.ir
  % add signals component-wise
  for i = 1:n.comp
    y.ir(:,:,d,1) = y.ir(:,:,d,1) + IR_fun_v4(...
      xtrue(i).m0, xtrue(i).t1, xtrue(i).t2,...
      P.ir.tr(d), P.ir.ti(d), P.ir.te(d,1),...
      'inveff', xtrue(i).inveff, 'kap', nu.tru.kap, 'wf', nu.tru.b0,...
      'flip_inv', P.ir.ainv(d), 'flip_ex', P.ir.aex(d), 'flip_ref', P.ir.aref(d));
  end
end

% collect noiseless se data
y.se = zeros(n.x, n.y, C.se, E.se);
for d = 1:C.se
  % add signals component-wise
  for i = 1:n.comp
    y.se(:,:,d,1) = y.se(:,:,d,1) + SE_fun_v4(...
      xtrue(i).m0, xtrue(i).t1, xtrue(i).t2,...
      P.se.tr(d), P.se.te(d,1),...
      'kap', nu.tru.kap, 'wf', nu.tru.b0, 'flip_ex', P.se.aex(d), 'flip_ref', P.se.aref(d));
  end
end

% add complex gaussian noise
n.ir = c.sig.ir * (randn(size(y.ir)) + 1i*randn(size(y.ir)));
n.se = c.sig.se * (randn(size(y.se)) + 1i*randn(size(y.se)));
y.ir = y.ir + n.ir;
y.se = y.se + n.se;

% set data weights based on profile
wght = struct([]);
for p = 1:n.p
  switch p
    case 1                              % (4 ir)
      wght(p).ir = [1 1 1 1]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
    case 2                              % (2 se), w/ te = (10,30) ms te
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [1 1 0 0]'   * ones(1,E.se);   
    case 3                              % (2 se), w/ te = (10,60) ms te
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [1 0 1 0]'   * ones(1,E.se);   
    case 4                              % (2 se), w/ te = (10,150) ms te
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [1 0 0 1]'   * ones(1,E.se);   
  end
end

% set dictionary parameters based on profile
dict = struct([]);
for p = 1:n.p
  dict(p).t1.n  = 300;
  dict(p).t1.ub = 10^3.5;
  dict(p).t1.lb = 10^1.5;
  dict(p).t1.sp = 'log';
  dict(p).t2.n  = 300;
  dict(p).t2.ub = 10^3.0;
  dict(p).t2.lb = 10^0.5;
  dict(p).t2.sp = 'log';
  dict(p).inveff.n = 100;
  switch p
    case 1
      dict(p).t2.n  = 1;
      dict(p).t2.ub = 10^6;
      dict(p).t2.lb = 10^6;
  end
end

% create masks from short te image
tmp = squeeze(abs(y.se(:,:,1)));
mask.thresh = 0.05;
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

% use magnitude data as appropriate
if bool.mag.ir
  y.ir = abs(y.ir);
end
if bool.mag.se
  y.se = abs(y.se);
end

% test different scan profiles
x = struct('m0', [], 't1', [], 't2', [], 'inveff', []);
t = struct('ml', [], 'rls', []);
for p = 1:n.p
  fprintf('\nTesting profile %d: (%u,%u) ir/se scans.\n',...
    p, sum(wght(p).ir(:,1)), sum(wght(p).se(:,1)));

  % set x0
  if p>1
    fprintf('  te = (%u,%u) ms.\n', P.se.te(wght(p).se>0));
    % for 2-se, set t1 from 4-ir ml est
    x0.t1 = x(1).t1.ml;
  else
    x0.t1 = [];
  end

  % parameter estimation
  opt.map = {...
    'mask.disp', mask.t,...
    'mask.est', mask.b,...
    'nu', nu.est,...
    'wght.ir', wght(p).ir,...
    'wght.se', wght(p).se,...
    'dict.t1', dict(p).t1,...
    'dict.t2', dict(p).t2,...
    'dict.inveff.n', dict(p).inveff.n,...
    'x0.t1', x0.t1,...
    'stop.iter', stop.iter,...
    'bool.mag.ir', bool.mag.ir,...
    'bool.mag.se', bool.mag.se,...
    'bool.chat', bool.chat,...
    'bool.norm', bool.norm,...
    'bool.disp', bool.disp,...
    'bool.precon', bool.precon};
  [tmp, t(p)] = mri_m0t1t2inveff_map(y, P, opt.map{:});

  % ensure similar structure across profiles
  x(p).m0 = tmp.m0;
  x(p).t1 = tmp.t1;
  x(p).t2 = tmp.t2;
  if isfield(tmp, 'inveff')
    x(p).inveff = tmp.inveff;
  else
    x(p).inveff.ml = ones(n.x,n.y);
    x(p).inveff.rls = ones(n.x,n.y);
  end
end

% omit profile 1 for analysis
if n.p == 4
  x(2).inveff = x(1).inveff;
  x(3).inveff = x(1).inveff;
  x(4).inveff = x(1).inveff;
  x(1) = [];
  n.p = 3;
end

% images
rng.m0 = [0 1];
rng.t1 = [500 1500];
rng.t2 = [20 120];
rng.inveff = [0.5 1];

opt.im.cmap = {'jet', 'gray'};
opt.im.prof = {'(10,30) ms', '(10,60) ms', '(10,150) ms'};
opt.im.t1 = {'T1 ML'};
opt.im.t2 = {'T2 ML'};

if bool.im
  for m = 1:length(opt.im.cmap)
    % t2 image
    tmp = [x(:).t2];
    figure; im('notick', 'row', 1, cat(3, tmp(:).ml), rng.t2, 'cbar', ' ');...
      colormap(opt.im.cmap{m});...
      tmp = colorbar; set(tmp, 'YTick', linspace(rng.t2(1), rng.t2(2), 6));
      hold on;
      text(col(n.x/2-1:n.x:5*n.x/2), zeros(n.p,1), col(opt.im.prof),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k');
      text(zeros(1,1), col(n.y/2-1:n.y:1*n.y/2), col(opt.im.t2),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
      hold off;
      if bool.pr
        tmp = sprintf('../../data/DigitalPhantom/2016-07-17,syn,multi-t2,%s.eps', opt.im.cmap{m});
        print('-depsc', tmp);
      end
  end
end

% % summary statistics
% stat.t1.ml.mean   = NaN(n.p, n.roi);
% stat.t1.ml.std    = NaN(n.p, n.roi);
% stat.t1.rls.mean  = NaN(n.p, n.roi);
% stat.t1.rls.std   = NaN(n.p, n.roi);
% stat.t2.ml.mean   = NaN(n.p, n.roi);
% stat.t2.ml.std    = NaN(n.p, n.roi);
% stat.t2.rls.mean  = NaN(n.p, n.roi);
% stat.t2.rls.std   = NaN(n.p, n.roi);
% for p = 1:n.p
%   for r = 1:n.roi
%     stat.t1.ml.mean(p,r)    = mean(masker(x(p).t1.ml, roi.mask(:,:,r)));
%     stat.t1.ml.std(p,r)     = std(masker(x(p).t1.ml, roi.mask(:,:,r)));
%     stat.t1.rls.mean(p,r)   = mean(masker(x(p).t1.rls, roi.mask(:,:,r)));
%     stat.t1.rls.std(p,r)    = std(masker(x(p).t1.rls, roi.mask(:,:,r)));
%     stat.t2.ml.mean(p,r)    = mean(masker(x(p).t2.ml, roi.mask(:,:,r)));
%     stat.t2.ml.std(p,r)     = std(masker(x(p).t2.ml, roi.mask(:,:,r)));
%     stat.t2.rls.mean(p,r)   = mean(masker(x(p).t2.rls, roi.mask(:,:,r)));
%     stat.t2.rls.std(p,r)    = std(masker(x(p).t2.rls, roi.mask(:,:,r)));
%   end
% end
% 
% % print
% if bool.pr
%   fid = fopen('../../data/Brain_05,31,16/2016-05-31,brain,multi-t2,stat', 'w');
%   fprintf(fid, '2016-05-31 parameter estimate statistics\n');
%   fprintf(fid, '\tcolumns are profiles\n');
%   fprintf(fid, '\trows are ROIs\n');
%   
%   fprintf(fid, '\nT1 ML:\n');
%   for r = 1:n.roi
%     fprintf(fid, '%s:\t', roi.label{r});
%     for p = 1:n.p
%       fprintf(fid, '%12.2f\t%c%8.2f\t',...
%         stat.t1.ml.mean(p,r), char(177), stat.t1.ml.std(p,r));
%     end
%     fprintf(fid, '\n');
%   end
%   
%   fprintf(fid, '\nT2 ML:\n');
%   for r = 1:n.roi
%     fprintf(fid, '%s:\t', roi.label{r});
%     for p = 1:n.p
%       fprintf(fid, '%12.2f\t%c%8.2f\t',...
%         stat.t2.ml.mean(p,r), char(177), stat.t2.ml.std(p,r));
%     end
%     fprintf(fid, '\n');
%   end
%   
%   fprintf(fid, '\n\n2016-05-31 parameter estimate run times (s)\n');
%   fprintf(fid, '\tcolumns are profiles\n');
%   fprintf(fid, '\trows are estimators\n');
%   fprintf(fid, '\nML: ');
%   for p = 1:n.p+1
%     fprintf(fid, '\t%12.2f', t(p).ml);
%   end
%   fprintf(fid, '\nRLS:');
%   for p = 1:n.p+1
%     fprintf(fid, '\t%12.2f', t(p).rls);
%   end
%   fclose(fid);
% end