% script HPDphant_2016_05_20_v2.m
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
%   2016-07-03      retry with mri_m0t1t2inveff_map.m for consistency
%   2016-07-12      added code to examine pooled std error of sample standard deviation
%   2016-08-02      new b1 rescaling from double-angle experiment

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/HPD_Phantom_05,20,16');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../crb');
addpath('../../etc');
addpath('../../others/mat-file-exchg');

% raw data path
% datdir = '/Volumes/a2/gnataraj/Data/2016,05,20_hpd_8ch';
datdir = '/Volumes/General Storage/Documents/school/michigan/research/data/2016,05,20_hpd_8ch';

% constant declarations
n.x = 256;
n.y = 256;
n.z = 8;
n.c = 8;
n.sl = 4;
n.p = 3;
n.r = 10;
n.k = sqrt((n.r-1)/2) * exp((gammaln((n.r-1)/2) - gammaln(n.r/2)));
n.v = 2 * ((n.r-1)/2 - exp(2 * (gammaln(n.r/2) - gammaln((n.r-1)/2))));

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 40;                             % mm

% acquisition parameters
P.sp.aex  = [15 5 15]' * pi/180;        % rad
P.sp.tr   = [12.2 12.2 13.9]';          % ms
C.sp      = length(P.sp.aex);          
E.sp      = 1;                          
P.sp.te   = ones(C.sp,E.sp) * 4.67;     % ms

P.de.aex  = [35 30 10 10]' * pi/180;    % rad
P.de.tr   = [24.4 17.5 28 17.5]';       % ms
C.de      = length(P.de.aex);           
E.de      = 2;
P.de.te   = ones(C.de,E.de) * 4.67;     % ms

% header options
bool.sc   = true;                       % uniformly scale kap from prior calibration
bool.sv   = true;                       % save coil-combined image data
bool.im   = true;                       % show images
bool.pr   = true;                       % print images and statistics

% coil combination options
opt.coil = {...
  'log2b', 10,...
  'fwhm', 1,...
  'nouter', 10,...
  'relPh', 1};

% b1 mapping options
opt.b1 = {...
  'coilOpt', opt.coil,...
  'reg.beta', 2^2,...
  'scale', 1,...
  'stop.iter', 2000,...
  'bool.chat', true,...
  'bool.disp', true};

% set data weights based on profile
wght = struct([]);
for p = 1:n.p
  switch p
    case 1                              % (2 sp, 1 de)
      wght(p).sp = [1 1 0]'     * ones(1,E.sp);
      wght(p).de = [0 1 0 0]'   * ones(1,E.de);
    case 2                              % (1 sp, 1 de) 
      wght(p).sp = [0 0 1]'     * ones(1,E.sp);
      wght(p).de = [0 0 1 0]'   * ones(1,E.de);
    case 3                              % (0 sp, 2 de)  
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [1 0 0 1]'   * ones(1,E.de);
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
end

% recon options
stop.iter = 300;                        % iteration stop criterion
bool.mag.ir = false;                    % use magnitude se-ir data/model
bool.mag.se = false;                    % use magnitude se data/model
bool.mag.sp = true;                     % use magnitude spgr data/model
bool.mag.de = true;                     % use magnitude dess data/model
bool.chat = true;                       % verbosity
bool.norm = false;                      % normalize coil-combined data before recon
bool.disp = false;                      % show image updates
bool.precon = true;                     % use preconditioner for rls estimation

% bloch-siegert b1 mapping
try 
  fprintf('\nTrying to load b1/kap maps...');
  load(sprintf('im,nu,sl-%u.mat', n.sl));
  fprintf('success!\n');
catch
  fprintf('failed: will reconstruct.\n'); 
  addpath([datdir '/pfile']);
  addpath([datdir '/mat']);
  addpath([datdir '/wav']);
  addpath('../../others/jfnielse/data');
  addpath('../../others/jfnielse/data/official/recon');
  addpath('../../others/jfnielse/data/official/img');
  
  % estimate separate b1 maps by repetition
  f.nframe = 9;
  f.necho = 2;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (672-127:672+128);
  f.yres = n.y;
  
  [~,P.bs.mag,~] = readwav('bsp.wav');    % G   
  P.bs.wrf  = 8;                          % kHz
  P.bs.dt   = 4e-3;                       % ms
  
  y = struct('im', [], 'k', []);
  nu = struct('b1', [], 'kap', []);
  for r = 1:n.r
    fprintf('\nProcessing repetition %u of %u...\n', r, n.r);
    
    % pfile details
    f.name = sprintf('P,bs-spgr-dess,rep%u.7', r);
    f.id = fopen(f.name, 'r', 'l'); 
      fseek(f.id,0,'bof');
    f.ver = str2double(num2str(fread(f.id,1,'float32')));
      fseek(f.id,0,'bof');
    f.hdr = read_rdb_hdr(f.id, f.ver);

    % bs coil data (frames 8-9)
    [y(r).im.bs.p, y(r).k.bs.p] = recon_gn(...
      f.name, 1+C.de+C.sp, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
    [y(r).im.bs.m, y(r).k.bs.m] = recon_gn(...
      f.name, 2+C.de+C.sp, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);

    % crop to single image slice and readout portion of interest
    % reorient to view superior -> inferior, right -> left            
    y(r).im.bs.p = ...                    % [n.x n.y n.c]
      permute(y(r).im.bs.p(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);                 
    y(r).im.bs.m = ...                    % [n.x n.y n.c]
      permute(y(r).im.bs.m(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);

    % bs b1 mapping
    nu(r) = mri_bs_kapb1_map(y(r).im.bs, P.bs, opt.b1{:});
  end
  
  % save b1/kap maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_05,20,16');
    save(sprintf('im,nu,sl-%u.mat', n.sl), 'nu');
    cd(curdir);
  end
end

% optional: uniform b1/kap scaling
if bool.sc
  scale = 1.050;
  for r = 1:n.r
    nu(r).b1.mom  = nu(r).b1.mom * scale;
    nu(r).b1.rls  = nu(r).b1.rls * scale;
    nu(r).kap.mom = nu(r).kap.mom * scale;
    nu(r).kap.rls = nu(r).kap.rls * scale;
  end
end

% spgr/dess coil combination
try 
  fprintf('\nTrying to load coil-combined spgr/dess data...');
  load(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl));
  fprintf('success!\n');
catch
  fprintf('failed.\nProceeding to extract and coil-combine raw data.\n');
  addpath([datdir '/pfile']);
  addpath([datdir '/mat']);
  addpath([datdir '/wav']);
  addpath('../../others/jfnielse/data');
  addpath('../../others/jfnielse/data/official/recon');
  addpath('../../others/jfnielse/data/official/img');
  
  % coil-combine maps by repetition
  f.nframe = 9;
  f.necho = 2;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (672-127:672+128);
  f.yres = n.y;
  
  % clear y to save memory
  if exist('y', 'var')
    clear('y');
  end
  
  ycc = struct('sp', [], 'de', []);
  smap = cell(n.r,n.p);
  for r = 1:n.r
    fprintf('\nGathering data for repetition %u of %u...\n', r, n.r);
    
    % pfile details
    f.name = sprintf('P,bs-spgr-dess,rep%u.7', r);
    f.id = fopen(f.name, 'r', 'l'); 
      fseek(f.id,0,'bof');
    f.ver = str2double(num2str(fread(f.id,1,'float32')));
      fseek(f.id,0,'bof');
    f.hdr = read_rdb_hdr(f.id, f.ver);

    % spgr coil data (frames 5-7)
    y.k.sp = NaN(f.rl, n.y, n.z, n.c, C.sp, E.sp);
    y.im.sp = NaN(f.rl, n.y, n.z, n.c, C.sp, E.sp);
    for i = 1:C.sp
      for j = 1:E.sp
        [y.im.sp(:,:,:,:,i,j), y.k.sp(:,:,:,:,i,j)] = recon_gn(...
          f.name, i+C.de, j, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
      end
    end

    % dess coil data (frames 1-4)
    y.k.de = NaN(f.rl, n.y, n.z, n.c, C.de, E.de);
    y.im.de = NaN(f.rl, n.y, n.z, n.c, C.de, E.de);
    for i = 1:C.de
      for j = 1:E.de
        [y.im.de(:,:,:,:,i,j), y.k.de(:,:,:,:,i,j)] = recon_gn(...
          f.name, i, j, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
      end
    end

    % crop to single image slice and readout portion of interest
    % reorient to view anterior -> posterior, right -> left           
    y.im.sp = ...                         % [n.x n.y C.sp E.sp n.c]
      permute(y.im.sp(f.cropx,:,n.sl,:,:,:), [2 1 5 6 4 3]); 
    y.im.de = ...                         % [n.x n.y C.de E.de n.c]
      permute(y.im.de(f.cropx,:,n.sl,:,:,:), [2 1 5 6 4 3]);  

    % separately coil-combine spgr and dess data together by profile
    ycc(r).sp = NaN(n.x, n.y, C.sp, E.sp);
    ycc(r).de = NaN(n.x, n.y, C.de, E.de);
    for p = 1:n.p
      fprintf('\nCoil-combining data for repetition %u, profile %u...\n', r, p);
      idx.sp = wght(p).sp(:,1)>0;
      idx.de = wght(p).de(:,1)>0;
      tmp = cat(3,...                     % [n.x n.y C(p).sp*E.sp+C(p).de*E.de n.c]
        reshape(y.im.sp(:,:,idx.sp,:,:), [n.x n.y sum(idx.sp)*E.sp n.c]),...
        reshape(y.im.de(:,:,idx.de,:,:), [n.x n.y sum(idx.de)*E.de n.c]));
      [tmp, smap{r,p}] = mri_multidata_coil_combine(tmp, opt.coil{:});
      ycc(r).sp(:,:,idx.sp,:) = ...
        reshape(tmp(:,:,1:(sum(idx.sp)*E.sp)), [n.x n.y sum(idx.sp) E.sp]);
      ycc(r).de(:,:,idx.de,:) = ...
        reshape(tmp(:,:,(sum(idx.sp)*E.sp+1):end), [n.x n.y sum(idx.de) E.de]);
      clear('idx');
    end
  end
  
  % save coil-combined images
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_05,20,16');
    save(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl), 'ycc', 'smap');
    cd(curdir);
  end
end

% m0,t1,t2 mapping by repetition
try
  fprintf('\nTrying to load parameter maps...');
  load(sprintf('../../data/HPD_Phantom_05,20,16/im,x,sl-%u.mat', n.sl));
  fprintf('success!\n');
catch
  fprintf('failed.\n');
  x = struct('m0', [], 't1', [], 't2', []);
  t = struct('ml', [], 'rls', []);
  mask = struct('t', [], 'b', [], 'thresh', []);
  
  for r = 1:n.r
    % fix mask across profiles 
    tmp = squeeze(abs(ycc(r).de(:,:,1,1)));
    mask(r).thresh = 0.05;
    mask(r).t = imfill(tmp > mask(r).thresh * max(col(tmp)), 'holes');
    mask(r).t = imdilate(~imdilate(~mask(r).t, strel('disk', 5)), strel('disk', 5));
    mask(r).b = imdilate(mask(r).t, strel('disk', 10));

    % use magnitude data as appropriate
    if bool.mag.sp
      ycc(r).sp = abs(ycc(r).sp);
    end
    if bool.mag.de
      ycc(r).de = abs(ycc(r).de);
    end

    % test different scan profiles
    for p = 1:n.p
      fprintf('\nTesting repetition %u of %u, profile %d of %d: (%u,%u) sp/de scans.\n',...
        r, n.r, p, n.p, sum(wght(p).sp(:,1)), sum(wght(p).de(:,1)));
      
      % parameter estimation
      opt.map = {...
        'mask.disp', mask(r).t,...
        'mask.est', mask(r).b,...
        'nu.kap', double(nu(r).kap.rls),...
        'wght.sp', wght(p).sp,...
        'wght.de', wght(p).de,...
        'dict.t1', dict(p).t1,...
        'dict.t2', dict(p).t2,...
        'stop.iter', stop.iter,...
        'bool.mag.sp', bool.mag.sp,...
        'bool.mag.de', bool.mag.de,...
        'bool.chat', bool.chat,...
        'bool.norm', bool.norm,...
        'bool.disp', bool.disp,...
        'bool.precon', bool.precon};
      [x(r,p), t(r,p)] = mri_m0t1t2inveff_map(ycc(r), P, opt.map{:});
    end
  end
  
  % save x maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_05,20,16');
    save(sprintf('im,x,sl-%u.mat', n.sl), 'x', 't');
    cd(curdir);
  end
end

% images
rng.m0 = [0 n.c];
rng.t1 = [0 2000];
rng.t2 = [0 500];
if bool.im
  % m0 images
  tmp = [x(:).m0];
  figure, im('row', n.p, cat(3, tmp(:).ml), rng.m0, 'cbar', ' '); colormap('jet');
  figure, im('row', n.p, cat(3, tmp(:).rls), rng.m0, 'cbar', ' '); colormap('jet');
  
  % t1 image
  tmp = [x(:).t1];
  figure, im('row', n.p, cat(3, tmp(:).ml), rng.t1, 'cbar', ' '); colormap('jet');
  figure, im('row', n.p, cat(3, tmp(:).rls), rng.t1, 'cbar', ' '); colormap('jet');
  
  % t2 image
  tmp = [x(:).t2];
  figure, im('row', n.p, cat(3, tmp(:).ml), rng.t2, 'cbar', ' '); colormap('jet');
  figure, im('row', n.p, cat(3, tmp(:).rls), rng.t2, 'cbar', ' '); colormap('jet');
end

% convert to matrix form for easier indexing
xmat.t1.ml  = NaN(n.x, n.y, n.r, n.p);
xmat.t1.rls = NaN(n.x, n.y, n.r, n.p);
xmat.t2.ml  = NaN(n.x, n.y, n.r, n.p);
xmat.t2.rls = NaN(n.x, n.y, n.r, n.p);
for r = 1:n.r
  for p = 1:n.p
    xmat.t1.ml(:,:,r,p)   = x(r,p).t1.ml;
    xmat.t1.rls(:,:,r,p)  = x(r,p).t1.rls;
    xmat.t2.ml(:,:,r,p)   = x(r,p).t2.ml;
    xmat.t2.rls(:,:,r,p)  = x(r,p).t2.rls;
  end
end

% define rois
roi.x = [130 162 182 182 163 132 100  80  80 100 110 152 153 110];
roi.y = [ 74  83 110 143 170 180 170 144 111  85 106 106 148 148];
n.roi = length(roi.x);
roi.r = 3*ones(1, n.roi);
roi.label = cell(n.roi,1);

% roi_masks
roi.mask = false(n.x, n.y, n.roi);
[roi.grid.x, roi.grid.y] = ndgrid(1:n.x, 1:n.y);
for r = 1:n.roi
  tmp = (roi.grid.x-roi.x(r)).^2 + (roi.grid.y-roi.y(r)).^2 <= roi.r(r)^2;
  roi.mask(:,:,r) = tmp;
  roi.label{r} = sprintf('Vial%2u', r);
end

% UMVU sample mean and sample standard deviation maps
t1.ml.smean  = squeeze(mean(xmat.t1.ml, 3));
t1.ml.sstd   = squeeze(std(xmat.t1.ml, 0, 3)) * n.k;
t1.rls.smean = squeeze(mean(xmat.t1.rls, 3));
t1.rls.sstd  = squeeze(std(xmat.t1.rls, 0, 3)) * n.k;
t2.ml.smean  = squeeze(mean(xmat.t2.ml, 3));
t2.ml.sstd   = squeeze(std(xmat.t2.ml, 0, 3)) * n.k;
t2.rls.smean = squeeze(mean(xmat.t2.rls, 3));
t2.rls.sstd  = squeeze(std(xmat.t2.rls, 0, 3)) * n.k;
 
% % summary statistics
% stat.t1.ml.mean     = NaN(n.p, n.roi);
% stat.t1.ml.std      = NaN(n.p, n.roi);
% stat.t1.rls.mean    = NaN(n.p, n.roi);
% stat.t1.rls.std     = NaN(n.p, n.roi);
% stat.t2.ml.mean     = NaN(n.p, n.roi);
% stat.t2.ml.std      = NaN(n.p, n.roi);
% stat.t2.rls.mean    = NaN(n.p, n.roi);
% stat.t2.rls.std     = NaN(n.p, n.roi);
% 
% for p = 1:n.p
%   for r = 1:n.roi
%     stat.t1.ml.mean(p,r)  = mean(masker(t1.ml.mean(:,:,p), roi.mask(:,:,r)));
%     stat.t1.ml.std(p,r)   = mean(masker(t1.ml.std(:,:,p), roi.mask(:,:,r)));
%     stat.t1.rls.mean(p,r) = mean(masker(t1.rls.mean(:,:,p), roi.mask(:,:,r)));
%     stat.t1.rls.std(p,r)  = mean(masker(t1.rls.std(:,:,p), roi.mask(:,:,r)));
%     stat.t2.ml.mean(p,r)  = mean(masker(t2.ml.mean(:,:,p), roi.mask(:,:,r)));
%     stat.t2.ml.std(p,r)   = mean(masker(t2.ml.std(:,:,p), roi.mask(:,:,r)));
%     stat.t2.rls.mean(p,r) = mean(masker(t2.rls.mean(:,:,p), roi.mask(:,:,r)));
%     stat.t2.rls.std(p,r)  = mean(masker(t2.rls.std(:,:,p), roi.mask(:,:,r)));
%   end
% end

% summary statistics
%   smean.pmean is pooled sample mean over reps
%   smean.perr  is pooled std err of sample mean over reps
%   sstd.pmean  is pooled sample std dev over reps
%   sstd.perr   is pooled std err of sample std dev over reps
stat.t1.ml.smean.pmean  = NaN(n.p, n.roi);
stat.t1.ml.smean.perr   = NaN(n.p, n.roi);
stat.t1.ml.sstd.pmean   = NaN(n.p, n.roi);
stat.t1.ml.sstd.perr    = NaN(n.p, n.roi);
stat.t1.rls.smean.pmean = NaN(n.p, n.roi);
stat.t1.rls.smean.perr  = NaN(n.p, n.roi);
stat.t1.rls.sstd.pmean  = NaN(n.p, n.roi);
stat.t1.rls.sstd.perr   = NaN(n.p, n.roi);
stat.t2.ml.smean.pmean  = NaN(n.p, n.roi);
stat.t2.ml.smean.perr   = NaN(n.p, n.roi);
stat.t2.ml.sstd.pmean   = NaN(n.p, n.roi);
stat.t2.ml.sstd.perr    = NaN(n.p, n.roi);
stat.t2.rls.smean.pmean = NaN(n.p, n.roi);
stat.t2.rls.smean.perr  = NaN(n.p, n.roi);
stat.t2.rls.sstd.pmean  = NaN(n.p, n.roi);
stat.t2.rls.sstd.perr   = NaN(n.p, n.roi);

for p = 1:n.p
  for r = 1:n.roi
    stat.t1.ml.smean.pmean(p,r)   = mean(masker(t1.ml.smean(:,:,p), roi.mask(:,:,r)));
    stat.t1.ml.sstd.pmean(p,r)    = mean(masker(t1.ml.sstd(:,:,p), roi.mask(:,:,r)));
    stat.t1.ml.smean.perr(p,r)    = stat.t1.ml.sstd.pmean(p,r) / sqrt(n.r);
    stat.t1.ml.sstd.perr(p,r)     = stat.t1.ml.sstd.pmean(p,r) * (n.k*sqrt(n.v/(n.r-1)));
    
    stat.t1.rls.smean.pmean(p,r)  = mean(masker(t1.rls.smean(:,:,p), roi.mask(:,:,r)));
    stat.t1.rls.sstd.pmean(p,r)   = mean(masker(t1.rls.sstd(:,:,p), roi.mask(:,:,r)));
    stat.t1.rls.smean.perr(p,r)   = stat.t1.rls.sstd.pmean(p,r) / sqrt(n.r);
    stat.t1.rls.sstd.perr(p,r)    = stat.t1.rls.sstd.pmean(p,r) * (n.k*sqrt(n.v/(n.r-1)));
    
    stat.t2.ml.smean.pmean(p,r)   = mean(masker(t2.ml.smean(:,:,p), roi.mask(:,:,r)));
    stat.t2.ml.sstd.pmean(p,r)    = mean(masker(t2.ml.sstd(:,:,p), roi.mask(:,:,r)));
    stat.t2.ml.smean.perr(p,r)    = stat.t2.ml.sstd.pmean(p,r) / sqrt(n.r);
    stat.t2.ml.sstd.perr(p,r)     = stat.t2.ml.sstd.pmean(p,r) * (n.k*sqrt(n.v/(n.r-1)));
    
    stat.t2.rls.smean.pmean(p,r)  = mean(masker(t2.rls.smean(:,:,p), roi.mask(:,:,r)));
    stat.t2.rls.sstd.pmean(p,r)   = mean(masker(t2.rls.sstd(:,:,p), roi.mask(:,:,r)));
    stat.t2.rls.smean.perr(p,r)   = stat.t2.rls.sstd.pmean(p,r) / sqrt(n.r);
    stat.t2.rls.sstd.perr(p,r)    = stat.t2.rls.sstd.pmean(p,r) * (n.k*sqrt(n.v/(n.r-1)));
  end
end
    
% empirical standard deviations
roi.t1.t = 5:7;
roi.t2.t = 6:7;
[roi.psi.t, iroi.t1.t, iroi.t2.t] = intersect(roi.t1.t, roi.t2.t);
roi.t1.b = 3:9;
roi.t2.b = 4:8;
[roi.psi.b, iroi.t1.b, iroi.t2.b] = intersect(roi.t1.b, roi.t2.b);

% relative scaling of parameters
roi.scale.t1 = 0.1;
roi.scale.t2 = 1;
stat.psi.ml.sstd.pmean = ...
  roi.scale.t1 * stat.t1.ml.sstd.pmean + ...
  roi.scale.t2 * stat.t2.ml.sstd.pmean;
stat.psi.ml.sstd.perr = ...
  roi.scale.t1 * stat.t1.ml.sstd.perr + ...
  roi.scale.t2 * stat.t2.ml.sstd.perr;

% worst-case ML pooled means of sample standard deviations
[worst.sigt1.pmean.t, tmp1] = max(stat.t1.ml.sstd.pmean(:,roi.t1.t), [], 2);
[worst.sigt1.pmean.b, tmp2] = max(stat.t1.ml.sstd.pmean(:,roi.t1.b), [], 2);
[worst.sigt2.pmean.t, tmp3] = max(stat.t2.ml.sstd.pmean(:,roi.t2.t), [], 2);
[worst.sigt2.pmean.b, tmp4] = max(stat.t2.ml.sstd.pmean(:,roi.t2.b), [], 2);
[worst.psi.pmean.t, tmp5] = max(stat.psi.ml.sstd.pmean(:,roi.psi.t), [], 2);
[worst.psi.pmean.b, tmp6] = max(stat.psi.ml.sstd.pmean(:,roi.psi.b), [], 2);

% worst-case ML pooled standard errors of sample standard deviations
worst.sigt1.pstd.t = NaN(n.p,1);
worst.sigt1.pstd.b = NaN(n.p,1);
worst.sigt2.pstd.t = NaN(n.p,1);
worst.sigt2.pstd.b = NaN(n.p,1);
worst.psi.pstd.t = NaN(n.p,1);
worst.psi.pstd.b = NaN(n.p,1);
for p = 1:n.p
  worst.sigt1.pstd.t(p) = stat.t1.ml.sstd.perr(p, tmp1(p)+min(roi.t1.t)-1);
  worst.sigt1.pstd.b(p) = stat.t1.ml.sstd.perr(p, tmp2(p)+min(roi.t1.b)-1);
  worst.sigt2.pstd.t(p) = stat.t2.ml.sstd.perr(p, tmp3(p)+min(roi.t1.t)-1);
  worst.sigt2.pstd.b(p) = stat.t2.ml.sstd.perr(p, tmp4(p)+min(roi.t1.b)-1);
  worst.psi.pstd.t(p) = stat.psi.ml.sstd.perr(p, tmp5(p)+min(roi.t1.t)-1);
  worst.psi.pstd.b(p) = stat.psi.ml.sstd.perr(p, tmp6(p)+min(roi.t1.b)-1);
end

% nist nmr statistics (aggregate measurements, no pooling needed)
stat.t1.nist.smean = ...
  [2480   2173   1907   1604    1332    1044    801.7  608.6  458.4  336.5  244.2   176.6   126.9   90.9]';
stat.t1.nist.sstd = ...
  [  10.8   14.7   10.3    7.2     0.8     3.2    1.70   1.03   0.33   0.18   0.09    0.09    0.03    0.05]'; 
stat.t2.nist.smean = ...
  [ 581.3  403.5  278.1  190.94  133.27   96.89  64.07  46.42  31.97  22.56  15.813  11.237   7.911  5.592]';
stat.t2.nist.sstd = ...
  [   0.39   0.55   0.28   0.011   0.073   0.049  0.034  0.014  0.083  0.012  0.0061  0.0057  0.0037 0.0055]';

% print
if bool.pr
  % sample mean +/- std err of sample mean
  fid = fopen('../../data/HPD_Phantom_05,20,16/2016-05-20,sample-mean', 'w');
  fprintf(fid, '2016-05-20 parameter estimate pooled statistics\n');
  fprintf(fid, '\tcolumns are profiles\n');
  fprintf(fid, '\trows are ROIs\n');
  
  fprintf(fid, '\nT1 ML sample mean +/- std err of sample mean:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.ml.smean.pmean(p,r), char(177), stat.t1.ml.smean.perr(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t1.nist.smean(r), char(177), stat.t1.nist.sstd(r));
  end
  fprintf(fid, '\nT1 RLS sample mean +/- std err of sample mean:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.rls.smean.pmean(p,r), char(177), stat.t1.rls.smean.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\nT2 ML sample mean +/- std err of sample mean:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.ml.smean.pmean(p,r), char(177), stat.t2.ml.smean.perr(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t2.nist.smean(r), char(177), stat.t2.nist.sstd(r));
  end
  fprintf(fid, '\nT2 RLS sample mean +/- std err of sample mean:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.rls.smean.pmean(p,r), char(177), stat.t2.rls.smean.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  fclose(fid);
  
  
  % sample std dev +/- std err of sample std dev
  fid = fopen('../../data/HPD_Phantom_05,20,16/2016-05-20,sample-std-dev', 'w');
  fprintf(fid, '2016-05-20 parameter estimate pooled statistics\n');
  fprintf(fid, '\tcolumns are profiles\n');
  fprintf(fid, '\trows are ROIs\n');
  
  fprintf(fid, '\nT1 ML sample std dev +/- std err of sample std dev:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.ml.sstd.pmean(p,r), char(177), stat.t1.ml.sstd.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\nT1 RLS sample std dev +/- std err of sample std dev:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.rls.sstd.pmean(p,r), char(177), stat.t1.rls.sstd.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\nT2 ML sample std dev +/- std err of sample std dev:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.ml.sstd.pmean(p,r), char(177), stat.t2.ml.sstd.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\nT2 RLS sample std dev +/- std err of sample std dev:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.rls.sstd.pmean(p,r), char(177), stat.t2.rls.sstd.perr(p,r));
    end
    fprintf(fid, '\n');
  end
  
  % print profile comparison
  fprintf(fid, '\n\nEmpirical performance summary over tight t1/t2 ranges\n');
  fprintf(fid, '\tcolumns are profiles\n');
  fprintf(fid, '\nsigwt1\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.sigt1.pmean.t(p), char(177), worst.sigt1.pstd.t(p));
  end
  fprintf(fid, '\nsigwt2\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.sigt2.pmean.t(p), char(177), worst.sigt2.pstd.t(p));
  end
  fprintf(fid, '\nPsiwt\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.psi.pmean.t(p), char(177), worst.psi.pstd.t(p));
  end
  
  fprintf(fid, '\n\nEmpirical performance summary over broad t1/t2 ranges\n');
  fprintf(fid, '\tcolumns are profiles\n\n');
  fprintf(fid, '\nsigwb1\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.sigt1.pmean.b(p), char(177), worst.sigt1.pstd.b(p));
  end
  fprintf(fid, '\nsigwb2\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.sigt2.pmean.b(p), char(177), worst.sigt2.pstd.b(p));
  end
  fprintf(fid, '\nPsiwb\t');
  for p = 1:n.p
    fprintf(fid, '%12.2f\t%c%8.2f\t',...
      worst.psi.pmean.b(p), char(177), worst.psi.pstd.b(p));
  end
  fclose(fid);
  
  
  fid = fopen('../../data/HPD_Phantom_05,20,16/2016-05-20,time', 'w');
  fprintf(fid, '2016-05-20 parameter estimate run times\n');
  fprintf(fid, '\tcolumns are profiles\n');
  fprintf(fid, '\trows are repetitions\n');
  
  fprintf(fid, '\nML:\n');
  for r = 1:n.r
    fprintf(fid, 'Rep%2u:\t', r);
    for p = 1:n.p
      fprintf(fid, '%12.2f\t', t(r,p).ml);
    end
    fprintf(fid, '\n');
  end
  
  fprintf(fid, '\nRLS:\n');
  for r = 1:n.r
    fprintf(fid, 'Rep%2u:\t', r);
    for p = 1:n.p
      fprintf(fid, '%12.2f\t', t(r,p).rls);
    end
    fprintf(fid, '\n');
  end
end
