% script HPDphant_2016_06_20.m
% joint recon of bs followed by all ir-se/se/spgr/dess scan profiles
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2016-05-26    protocol prepared
%   v0.2    2016-06-20    data acquired
%   v1.1    2016-06-20    modified from 2016-06-17; see log for details
%   v1.2    2016-06-21    roi analysis; studying ir-se/se t1 bias
%   v1.3    2016-06-27    using mri_m0t1t2inveff_map.m to avoid t1 bias due to imperfect inversion
%   v1.4    2016-06-29    now using sequential t1, then t2 estimation from ir, then se
%   v1.5    2016-07-03    coil-combine separately by profile
%   v1.6    2016-08-02    new b1 rescaling from double-angle experiment

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/HPD_Phantom_06,20,16');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../etc');
addpath('../../others/mat-file-exchg');

% constant declarations
n.x = 256;
n.y = 256;
n.z = 8;
n.c = 8;
n.sl = 6;
n.p = 5;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

% acquisition parameters
P.ir.ti   = [50 150 450 1350]';         % ms
C.ir      = length(P.ir.ti);
E.ir      = 1;
P.ir.tr   = ones(C.ir,1) * 1400;        % ms
P.ir.te   = ones(C.ir,E.ir) * 14;       % ms
P.ir.ainv = ones(C.ir,1) * pi;          % rad
P.ir.aex  = ones(C.ir,1) * pi/2;        % rad
P.ir.aref = ones(C.ir,1) * pi;          % rad

P.se.te   = [10 30 60 150]';            % ms
C.se      = length(P.se.te);
E.se      = 1;         
P.se.tr   = ones(C.se,1) * 1000;        % ms
P.se.aex  = ones(C.se,1) * pi/2;        % rad
P.se.aref = ones(C.se,1) * pi;          % rad

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
bool.cmp  = true;                       % make t1/t2 comparison plots
bool.sig  = false;                       % examine signal models

% mask threshold
mask.thresh = 0.05;

% coil combination options
opt.coil = {...
  'log2b', 10,...
  'fwhm', 1,...
  'nouter', 10,...
  'relPh', 1,...
  'disp', 0};

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
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [1 1 0]'     * ones(1,E.sp);
      wght(p).de = [0 1 0 0]'   * ones(1,E.de);
    case 2                              % (1 sp, 1 de)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [0 0 1]'     * ones(1,E.sp);
      wght(p).de = [0 0 1 0]'   * ones(1,E.de);
    case 3                              % (0 sp, 2 de)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);                   
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);   
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [1 0 0 1]'   * ones(1,E.de);
    case 4                              % (4 ir, 0 se)
      wght(p).ir = [1 1 1 1]'   * ones(1,E.ir);
      wght(p).se = [0 0 0 0]'   * ones(1,E.se);
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [0 0 0 0]'   * ones(1,E.de);
    case 5                              % (0 ir, 4 se)
      wght(p).ir = [0 0 0 0]'   * ones(1,E.ir);
      wght(p).se = [1 1 1 1]'   * ones(1,E.se);   
      wght(p).sp = [0 0 0]'     * ones(1,E.sp);
      wght(p).de = [0 0 0 0]'   * ones(1,E.de);
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
    case 4
      dict(p).t2.n  = 1;
      dict(p).t2.ub = 10^6;
      dict(p).t2.lb = 10^6;
  end
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
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/mat');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/wav');
  addpath('../../others/jfnielse/data');
  addpath('../../others/jfnielse/data/official/recon');
  addpath('../../others/jfnielse/data/official/img');
  
  % pfile details
  f.name = 'P,bs-spgr-dess.7';
  f.id = fopen(f.name, 'r', 'l'); 
    fseek(f.id,0,'bof');
  f.ver = str2double(num2str(fread(f.id,1,'float32')));
    fseek(f.id,0,'bof');
  f.hdr = read_rdb_hdr(f.id, f.ver);
  f.nframe = 9;
  f.necho = 2;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (512-127:512+128);
  f.yres = n.y;
  
  % bs coil data (frames 8-9)
  [y.im.bs.p, y.k.bs.p] = recon_gn(...
    f.name, 1+C.de+C.sp, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
  [y.im.bs.m, y.k.bs.m] = recon_gn(...
    f.name, 2+C.de+C.sp, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
  
  % crop to single image slice and readout portion of interest
  % reorient to view superior -> inferior, right -> left            
  y.im.bs.p = ...                       % [n.x n.y n.c]
    permute(y.im.bs.p(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);                 
  y.im.bs.m = ...                       % [n.x n.y n.c]
    permute(y.im.bs.m(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);

  % bs b1 mapping
  [~,P.bs.mag,~] = readwav('bsp.wav');  % G   
  P.bs.wrf  = 8;                        % kHz
  P.bs.dt   = 4e-3;                     % ms
  nu = mri_bs_kapb1_map(y.im.bs, P.bs, opt.b1{:});
  
  % save b1/kap maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,20,16');
    save(sprintf('im,nu,sl-%u.mat', n.sl), 'nu');
    cd(curdir);
  end
end

% optional: uniform b1/kap scaling
if bool.sc
  scale = 1.050;
  nu.b1.mom = nu.b1.mom * scale;
  nu.b1.rls = nu.b1.rls * scale;
  nu.kap.mom = nu.kap.mom * scale;
  nu.kap.rls = nu.kap.rls * scale;
end
% nu.kap.rls = ones(n.x,n.y);

% ir-se/se coil combination
try
  fprintf('\nTrying to load coil-combined ir/se data...');
  load('im,irse-se,coil-comb.mat');
  fprintf('success!\n');
catch
  fprintf('failed.\nProceeding to extract coil-combine raw data.\n');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/mat');
  
  % pfile details
  p.fr = n.x;                           % readout length
  p.ep = 0;                             % extended precision 
  p.hosp = 0;                           % hospital data
  p.basetr = 1;                         % baseline tr used
  p.coil = 1:n.c;                       % which coils to extract
  p.sl = 1;                             % which slices to extract
  p.echo = 1;                           % which echoes to extract
  
  % inversion-recovery coil data
  y.k.ir  = NaN(n.x, n.y, n.c, C.ir, E.ir);
  y.im.ir = NaN(n.x, n.y, n.c, C.ir, E.ir);
  for i = 1:C.ir
    for j = 1:E.ir
      p.fn = sprintf('P,ir,ti-%dms.7', P.ir.ti(i));
      [y.k.ir(:,:,:,i,j),~,~] = rawload_jfn(...
        p.fn, p.fr, p.ep, p.hosp, p.basetr, p.coil, p.sl, p.echo);
      for c = 1:n.c
        % shift field-of-view in phase-encoding direction
        y.k.ir(:,2:2:end,c,i,j) = -y.k.ir(:,2:2:end,c,i,j);

        % images are 2d ifft of data
        y.im.ir(:,:,c,i,j) = fftshift(ifft2(ifftshift(y.k.ir(:,:,c,i,j))));
      end
    end
  end
  
  % spin echo coil data
  y.k.se  = NaN(n.x, n.y, n.c, C.se, E.se);
  y.im.se = NaN(n.x, n.y, n.c, C.se, E.se);
  for i = 1:C.se
    for j = 1:E.se
      p.fn = sprintf('P,se,te-%dms.7', P.se.te(i));
      [y.k.se(:,:,:,i,j),~,~] = rawload_jfn(...
        p.fn, p.fr, p.ep, p.hosp, p.basetr, p.coil, p.sl, p.echo);
      for c = 1:n.c
        % shift field-of-view in phase-encoding direction
        y.k.se(:,2:2:end,c,i,j) = -y.k.se(:,2:2:end,c,i,j);

        % images are 2d ifft of data
        y.im.se(:,:,c,i,j) = fftshift(ifft2(ifftshift(y.k.se(:,:,c,i,j))));
      end
    end
  end
  
  % reorient to view superior -> inferior, right -> left
  y.im.ir = ...                         % [n.x n.y C.ir E.ir n.c]
    permute(flip(y.im.ir, 1), [2 1 4 5 3]);
  y.im.se = ...                         % [n.x n.y C.se E.se n.c]
    permute(flip(y.im.se, 1), [2 1 4 5 3]);
  
  % coil-combine ir and se data together
  tmp = cat(3,...                       % [n.x n.y C.ir*E.ir+C.se*E.se n.c]
    reshape(y.im.ir, [n.x n.y C.ir*E.ir n.c]),...
    reshape(y.im.se, [n.x n.y C.se*E.se n.c]));
  [tmp, smap] = mri_multidata_coil_combine(tmp, opt.coil{:});
  ycc.ir = reshape(tmp(:,:,1:(C.ir*E.ir)), [n.x n.y C.ir E.ir]);
  ycc.se = reshape(tmp(:,:,(C.ir*E.ir+1):end), [n.x n.y C.se E.se]);
  
  % save coil-combined images 
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,20,16');
    save('im,irse-se,coil-comb.mat', 'ycc', 'smap');
    cd(curdir);
  end
end

% spgr/dess coil combination
try 
  fprintf('\nTrying to load coil-combined spgr/dess data...');
  tmp = ycc;
  clear('ycc');
  load(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl));
  ycc = catstruct(tmp, ycc);
  fprintf('success!\n');
catch
  fprintf('failed.\nProceeding to extract and coil-combine raw data.\n');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/mat');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,20_hpd_8ch/wav');
  addpath('../../others/jfnielse/data');
  addpath('../../others/jfnielse/data/official/recon');
  addpath('../../others/jfnielse/data/official/img');
  
  % pfile details
  f.name = 'P,bs-spgr-dess.7';
  f.id = fopen(f.name, 'r', 'l'); 
    fseek(f.id,0,'bof');
  f.ver = str2double(num2str(fread(f.id,1,'float32')));
    fseek(f.id,0,'bof');
  f.hdr = read_rdb_hdr(f.id, f.ver);
  f.nframe = 9;
  f.necho = 2;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (512-127:512+128);
  f.yres = n.y;
  
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
  ycc.sp = NaN(n.x, n.y, C.sp, E.sp);
  ycc.de = NaN(n.x, n.y, C.de, E.de);
  smap = cell(3,1);
  for p = 1:3
    fprintf('\nCoil-combining data for profile %u...\n', p);
    idx.sp = wght(p).sp(:,1)>0;
    idx.de = wght(p).de(:,1)>0;
    tmp = cat(3,...                     % [n.x n.y C(p).sp*E.sp+C(p).de*E.de n.c]
      reshape(y.im.sp(:,:,idx.sp,:,:), [n.x n.y sum(idx.sp)*E.sp n.c]),...
      reshape(y.im.de(:,:,idx.de,:,:), [n.x n.y sum(idx.de)*E.de n.c]));
    [tmp, smap{p}] = mri_multidata_coil_combine(tmp, opt.coil{:});
    ycc.sp(:,:,idx.sp,:) = ...
      reshape(tmp(:,:,1:(sum(idx.sp)*E.sp)), [n.x n.y sum(idx.sp) E.sp]);
    ycc.de(:,:,idx.de,:) = ...
      reshape(tmp(:,:,(sum(idx.sp)*E.sp+1):end), [n.x n.y sum(idx.de) E.de]);
    clear('idx');
  end
  
  % save coil-combined images
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,20,16');
    save(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl), 'ycc', 'smap');
    cd(curdir);
  end
  
  % append se/ir coil-combined images
  tmp = ycc;
  clear('ycc');
  load('im,irse-se,coil-comb.mat');
  ycc = catstruct(tmp, ycc);
end

% create masks from short te image
tmp = squeeze(abs(ycc.se(:,:,1)));
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

% use magnitude data as appropriate
if bool.mag.ir
  ycc.ir = abs(ycc.ir);
end
if bool.mag.se
  ycc.se = abs(ycc.se);
end
if bool.mag.sp
  ycc.sp = abs(ycc.sp);
end
if bool.mag.de
  ycc.de = abs(ycc.de);
end

% test different scan profiles
try 
  load(sprintf('../../data/HPD_Phantom_06,20,16/im,x,sl-%u.mat', n.sl));
catch
  x = struct('m0', [], 't1', [], 't2', [], 'inveff', []);
  t = struct('ml', [], 'rls', []);
  for p = 1:n.p
    fprintf('\nTesting profile %d: (%u,%u,%u,%u) ir/se/sp/de scans.\n',...
      p, sum(wght(p).ir(:,1)), sum(wght(p).se(:,1)), sum(wght(p).sp(:,1)), sum(wght(p).de(:,1)));

    % set x0
    if p==5
      % for 4-se, set t1 from 4-ir ml est
      x0.t1 = x(4).t1.ml;
    else
      x0.t1 = [];
    end

    % parameter estimation
    opt.map = {...
      'mask.disp', mask.t,...
      'mask.est', mask.b,...
      'nu.kap', double(nu.kap.rls),...
      'wght', wght(p),...
      'dict.t1', dict(p).t1,...
      'dict.t2', dict(p).t2,...
      'dict.inveff.n', dict(p).inveff.n,...
      'x0.t1', x0.t1,...
      'stop.iter', stop.iter,...
      'bool.mag.ir', bool.mag.ir,...
      'bool.mag.se', bool.mag.se,...
      'bool.mag.sp', bool.mag.sp,...
      'bool.mag.de', bool.mag.de,...
      'bool.chat', bool.chat,...
      'bool.norm', bool.norm,...
      'bool.disp', bool.disp,...
      'bool.precon', bool.precon};
    [tmp, t(p)] = mri_m0t1t2inveff_map(ycc, P, opt.map{:});

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
  
  % save x maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,20,16');
    save(sprintf('im,x,sl-%u.mat', n.sl), 'x', 't');
    cd(curdir);
  end
end

% define rois
roi.x = [129 160 180 179 160 128  97  77  77  98 107 150 150 107];
roi.y = [ 75  85 112 145 171 181 171 143 111  84 106 107 149 149];
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
roi.all = sum(roi.mask, 3) > 0;

% omit profile 4 for analysis
if n.p == 5
  x(5).inveff = x(4).inveff;
  x(4) = [];
  n.p = 4;
end

% images
rng.m0 = [0 n.c];
rng.t1 = [0 2000];
rng.t2 = [0 500];
rng.inveff = [0.5 1];

opt.im.cmap = {'jet', 'gray'};
opt.im.prof = {'(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)'};
opt.im.t1 = {'T1 ML', 'T1 RLS'};
opt.im.t2 = {'T2 ML', 'T2 RLS'};

% colors
color.o = [1 0.6 0.2];                  % orange
color.f = [34 139 34]/256;              % forest green
color.m = [102 0 102]/256;              % maroon
color.l = [229 204 255]/256;            % lavendar
color.s = [255 0 127]/256;              % salmon
color.p = [255 153 255]/256;            % pink

% tight and loose rois
roi.t1.t = 5:7;
roi.t2.t = 6:7;
roi.t1.b = 3:9;
roi.t2.b = 4:8;

if bool.im
  for m = 1:length(opt.im.cmap)
    % t1 image
    tmp = [x(:).t1];
    figure; im('notick', 'row', 2, cat(3, tmp(:).ml, tmp(:).rls), rng.t1, 'cbar', ' ');...
      colormap(opt.im.cmap{m});...
      tmp = colorbar; set(tmp, 'YTick', linspace(rng.t1(1), rng.t1(2), 6));...
      hold on;...
      text(col(n.x/2-1:n.x:7*n.x/2), zeros(n.p,1), col(opt.im.prof),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k');
      text(zeros(2,1), col(n.y/2-1:n.y:3*n.y/2), col(opt.im.t1),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
      if strcmp(opt.im.cmap{m}, 'gray')
        for r = 1:n.roi
          if min(roi.t1.t) <= r && r <= max(roi.t1.t)
            tmp = color.o;
          elseif min(roi.t1.b) <= r && r <= max(roi.t1.b)
            tmp = 'y';
          elseif r < min(roi.t1.b)
            tmp = 'k';
          else
            tmp = 'w';
          end
          text(3*n.x+roi.x(r), roi.y(r), num2str(r),...
            'VerticalAlignment', 'middle',...
            'HorizontalAlignment', 'center',...
            'FontSize', 6,...
            'Color', tmp,... 
            'FontWeight', 'bold');
        end
      end
      hold off;
      if bool.pr
        tmp = sprintf('../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t1,%s.eps', opt.im.cmap{m});
        print('-depsc', tmp);
      end
      
    % t2 image
    tmp = [x(:).t2];
    figure; im('notick', 'row', 2, cat(3, tmp(:).ml, tmp(:).rls), rng.t2, 'cbar', ' ');...
      colormap(opt.im.cmap{m});...
      tmp = colorbar; set(tmp, 'YTick', linspace(rng.t2(1), rng.t2(2), 6));
      hold on;
      text(zeros(2,1), col(n.y/2-1:n.y:3*n.y/2), col(opt.im.t2),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
      if strcmp(opt.im.cmap{m}, 'gray')
        for r = 1:n.roi
          if min(roi.t2.t) <= r && r <= max(roi.t2.t)
            tmp = color.o;
          elseif min(roi.t2.b) <= r && r <= max(roi.t2.b)
            tmp = 'y';
          elseif r < min(roi.t2.b)
            tmp = 'k';
          else
            tmp = 'w';
          end
          text(3*n.x+roi.x(r), roi.y(r), num2str(r),...
            'VerticalAlignment', 'middle',...
            'HorizontalAlignment', 'center',...
            'FontSize', 6,...
            'Color', tmp,... 
            'FontWeight', 'bold');
        end
      end
      hold off;
      if bool.pr
        tmp = sprintf('../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t2,%s.eps', opt.im.cmap{m});
        print('-depsc', tmp);
      end
  end
else
  tmp = [x(:).m0];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.m0, 'cbar'); 
  tmp = [x(:).t1];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.t1, 'cbar');
  tmp = [x(:).t2];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.t2, 'cbar');
  tmp = [x(:).inveff];
  figure, im('row', 2, cat(4, cat(3, tmp(:).ml), cat(3, tmp(:).rls)), rng.inveff, 'cbar');
end

% summary statistics
stat.t1.ml.mean     = NaN(n.p, n.roi);
stat.t1.ml.std      = NaN(n.p, n.roi);
stat.t1.rls.mean    = NaN(n.p, n.roi);
stat.t1.rls.std     = NaN(n.p, n.roi);
stat.t2.ml.mean     = NaN(n.p, n.roi);
stat.t2.ml.std      = NaN(n.p, n.roi);
stat.t2.rls.mean    = NaN(n.p, n.roi);
stat.t2.rls.std     = NaN(n.p, n.roi);
stat.logt1.ml.mean  = NaN(n.p, n.roi);
stat.logt1.ml.std   = NaN(n.p, n.roi);
stat.logt1.rls.mean = NaN(n.p, n.roi);
stat.logt1.rls.std  = NaN(n.p, n.roi);
stat.logt2.ml.mean  = NaN(n.p, n.roi);
stat.logt2.ml.std   = NaN(n.p, n.roi);
stat.logt2.rls.mean = NaN(n.p, n.roi);
stat.logt2.rls.std  = NaN(n.p, n.roi);
for p = 1:n.p
  for r = 1:n.roi
    stat.m0.ml.mean(p,r)      = mean(masker(x(p).m0.ml, roi.mask(:,:,r)));
    
    stat.t1.ml.mean(p,r)      = mean(masker(x(p).t1.ml, roi.mask(:,:,r)));
    stat.t1.ml.std(p,r)       = std(masker(x(p).t1.ml, roi.mask(:,:,r)));
    stat.t1.rls.mean(p,r)     = mean(masker(x(p).t1.rls, roi.mask(:,:,r)));
    stat.t1.rls.std(p,r)      = std(masker(x(p).t1.rls, roi.mask(:,:,r)));
    stat.t2.ml.mean(p,r)      = mean(masker(x(p).t2.ml, roi.mask(:,:,r)));
    stat.t2.ml.std(p,r)       = std(masker(x(p).t2.ml, roi.mask(:,:,r)));
    stat.t2.rls.mean(p,r)     = mean(masker(x(p).t2.rls, roi.mask(:,:,r)));
    stat.t2.rls.std(p,r)      = std(masker(x(p).t2.rls, roi.mask(:,:,r)));
    
    stat.logt1.ml.mean(p,r)   = mean(masker(log10(x(p).t1.ml), roi.mask(:,:,r)));
    stat.logt1.ml.std(p,r)    = std(masker(log10(x(p).t1.ml), roi.mask(:,:,r)));
    stat.logt1.rls.mean(p,r)  = mean(masker(log10(x(p).t1.rls), roi.mask(:,:,r)));
    stat.logt1.rls.std(p,r)   = std(masker(log10(x(p).t1.rls), roi.mask(:,:,r)));
    stat.logt2.ml.mean(p,r)   = mean(masker(log10(x(p).t2.ml), roi.mask(:,:,r)));
    stat.logt2.ml.std(p,r)    = std(masker(log10(x(p).t2.ml), roi.mask(:,:,r)));
    stat.logt2.rls.mean(p,r)  = mean(masker(log10(x(p).t2.rls), roi.mask(:,:,r)));
    stat.logt2.rls.std(p,r)   = std(masker(log10(x(p).t2.rls), roi.mask(:,:,r)));
  end
end

% nist statistics
stat.t1.nist.mean = ...
  [2480   2173   1907   1604    1332    1044    801.7  608.6  458.4  336.5  244.2   176.6   126.9   90.9]';
stat.t1.nist.std = ...
  [  10.8   14.7   10.3    7.2     0.8     3.2    1.70   1.03   0.33   0.18   0.09    0.09    0.03    0.05]'; 
stat.t2.nist.mean = ...
  [ 581.3  403.5  278.1  190.94  133.27   96.89  64.07  46.42  31.97  22.56  15.813  11.237   7.911  5.592]';
stat.t2.nist.std = ...
  [   0.39   0.55   0.28   0.011   0.073   0.049  0.034  0.014  0.083  0.012  0.0061  0.0057  0.0037 0.0055]';

% log-transformed nist-statistics
% not quite right... should rescale data but not possible for nist
stat.t1.nist.cfvar    = stat.t1.nist.mean ./ stat.t1.nist.std;
stat.logt1.nist.mean  = log10(stat.t1.nist.mean);
stat.logt1.nist.std   = stat.logt1.nist.mean ./ stat.t1.nist.cfvar;
stat.t2.nist.cfvar    = stat.t2.nist.mean ./ stat.t2.nist.std;
stat.logt2.nist.mean  = log10(stat.t2.nist.mean);
stat.logt2.nist.std   = stat.logt2.nist.mean ./ stat.t2.nist.cfvar;

% print
if bool.pr
  fid = fopen('../../data/HPD_Phantom_06,20,16/2016-06-20,stat', 'w');
  fprintf(fid, '2016-06-20 parameter estimate statistics\n');
  fprintf(fid, '\tcolumns are profiles, (last is NIST) \n');
  fprintf(fid, '\trows are ROIs\n');
  
  fprintf(fid, '\nT1 ML:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.ml.mean(p,r), char(177), stat.t1.ml.std(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t1.nist.mean(r), char(177), stat.t1.nist.std(r));
  end
  
  fprintf(fid, '\nT1 RLS:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t1.rls.mean(p,r), char(177), stat.t1.rls.std(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t1.nist.mean(r), char(177), stat.t1.nist.std(r));
  end
  
  fprintf(fid, '\nT2 ML:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.ml.mean(p,r), char(177), stat.t2.ml.std(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t2.nist.mean(r), char(177), stat.t2.nist.std(r));
  end
  
  fprintf(fid, '\nT2 RLS:\n');
  for r = 1:n.roi
    fprintf(fid, '%s:\t', roi.label{r});
    for p = 1:n.p
      fprintf(fid, '%12.2f\t%c%8.2f\t',...
        stat.t2.rls.mean(p,r), char(177), stat.t2.rls.std(p,r));
    end
    fprintf(fid, '%12.2f\t%c%8.2f\n',...
      stat.t2.nist.mean(r), char(177), stat.t2.nist.std(r));
  end
  
  fprintf(fid, '\n\n2016-06-20 parameter estimate run times\n');
  fprintf(fid, '\tcolumns are profiles: last two are (4IR), then (4SE)\n');
  fprintf(fid, '\trows are estimators\n');
  fprintf(fid, '\nML: ');
  for p = 1:n.p+1
    fprintf(fid, '\t%12.2f', t(p).ml);
  end
  fprintf(fid, '\nRLS:');
  for p = 1:n.p+1
    fprintf(fid, '\t%12.2f', t(p).rls);
  end
  fclose(fid);
end

% t1/t2 roi fill areas
roi.t1.rng.t  = [800 1400];
roi.t2.rng.t  = [50 120];
roi.t1.rng.b  = [400 2000];
roi.t2.rng.b  = [40 200];

roi.t1.x.t    = log10(kron(roi.t1.rng.t, [1 1]));
roi.t1.y.t    = log10([roi.t1.rng.t fliplr(roi.t1.rng.t)]);
roi.t2.x.t    = log10(kron(roi.t2.rng.t, [1 1]));
roi.t2.y.t    = log10([roi.t2.rng.t fliplr(roi.t2.rng.t)]);
roi.t1.x.b    = log10(kron(roi.t1.rng.b, [1 1]));
roi.t1.y.b    = log10([roi.t1.rng.b fliplr(roi.t1.rng.b)]);
roi.t2.x.b    = log10(kron(roi.t2.rng.b, [1 1]));
roi.t2.y.b    = log10([roi.t2.rng.b fliplr(roi.t2.rng.b)]);

% fig options
fig.c = {'b', 'r', color.f, color.s, color.m};
fig.m = {'o', 'v', 'd', 's', '^'};
fig.f = {'c', 'm', 'g', color.p, color.l};

% fig.tick.t1 = [1.9:0.2:3.5];
% fig.ticklab.t1 = cell(1,length(fig.tick.t1));
% for i=1:length(fig.tick.t1)
%   fig.ticklab.t1{i} = sprintf('%.4g', 10^fig.tick.t1(i));
% end
% 
% fig.tick.t2 = [0.7:0.3:2.8];
% fig.ticklab.t2 = cell(1,length(fig.tick.t2));
% for i=1:length(fig.tick.t2)
%   fig.ticklab.t2{i} = sprintf('%.4g', 10^fig.tick.t2(i));
% end

fig.tick.t1 = log10([80:10:100 200:100:1000 2000:1000:3000]);
fig.ticklab.t1 = cell(1,length(fig.tick.t1));
for i=1:length(fig.tick.t1)
  if rem(fig.tick.t1(i),1)==0
    fig.ticklab.t1{i} = sprintf('%u', 10^fig.tick.t1(i));
  else
    fig.ticklab.t1{i} = ' ';
  end
end

fig.tick.t2 = log10([5:1:10 20:10:100 200:100:700]);
fig.ticklab.t2 = cell(1,length(fig.tick.t2));
for i=1:length(fig.tick.t2)
  if rem(fig.tick.t2(i),1)==0
    fig.ticklab.t2{i} = sprintf('%u', 10^fig.tick.t2(i));
  else 
    fig.ticklab.t2{i} = ' ';
  end
end

if bool.cmp
  % fig 1: t1,ml
  figure; hold on;
  fill(roi.t1.x.b, roi.t1.y.b, 'y');
  fill(roi.t1.x.t, roi.t1.y.t, color.o);
  t1{1} = plot(log10(1:10000), log10(1:10000), 'k--', 'LineWidth', 1);
  for p = 1:n.p
    t1{p+1} = errorbarxy(...
      stat.logt1.nist.mean, col(stat.logt1.ml.mean(p,:)),...
      stat.logt1.nist.std, col(stat.logt1.ml.std(p,:)),...
      'Color', fig.c{p},...
      'LineStyle', 'none',...
      'Marker', fig.m{p},...
      'LineWidth', 1.0,...
      'MarkerSize', 8,...
      'MarkerFaceColor', fig.f{p});
    t1{p+1} = col(t1{p+1});
  end
  roi.txton = 1:n.roi;
  text(stat.logt1.nist.mean(roi.txton), col(stat.logt1.ml.mean(1,roi.txton)),...
    num2str(roi.txton'),...
    'VerticalAlignment', 'top',...
    'HorizontalAlignment', 'left',...
    'FontSize', 16);
  hold off; axis tight; axis square; grid on;
  set(gca, 'XTick', fig.tick.t1);
  set(gca, 'YTick', fig.tick.t1);
  set(gca, 'XLim', minmax(fig.tick.t1));
  set(gca, 'YLim', minmax(fig.tick.t1));
  set(gca, 'XTickLabel', fig.ticklab.t1);
  set(gca, 'YTickLabel', fig.ticklab.t1);
  set(gca, 'Layer', 'top');
  xlabel('NIST T1 Estimates (ms)', 'FontSize', 16);
  xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
  ylabel('Candidate Profile T1 ML Estimates (ms)', 'FontSize', 16);
  tmp = [t1{2:end}];
  tmp = legend([t1{1} tmp(1,:)],...
    'Ideal', '(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)', 'Location', 'NW');
  set(tmp, 'FontSize', 18);
  if bool.pr
    print('-depsc', '../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t1-ml-compare.eps');
  end
  
  % fig 2: t2,ml
  figure; hold on;
  fill(roi.t2.x.b, roi.t2.y.b, 'y');
  fill(roi.t2.x.t, roi.t2.y.t, color.o);
  t2{1} = plot(log10(1:1000), log10(1:1000), 'k--', 'LineWidth', 1);
  for p = 1:n.p
    t2{p+1} = errorbarxy(...
      stat.logt2.nist.mean, col(stat.logt2.ml.mean(p,:)),...
      stat.logt2.nist.std, col(stat.logt2.ml.std(p,:)),...
      'Color', fig.c{p},...
      'LineStyle', 'none',...
      'Marker', fig.m{p},...
      'LineWidth', 1.0,...
      'MarkerSize', 8,...
      'MarkerFaceColor', fig.f{p});
    t2{p+1} = col(t2{p+1});
  end
  roi.txton = 1:n.roi;
  text(stat.logt2.nist.mean(roi.txton), col(stat.logt2.ml.mean(1,roi.txton)),...
    num2str(roi.txton'),...
    'VerticalAlignment', 'top',...
    'HorizontalAlignment', 'left',...
    'FontSize', 16);
  hold off; axis tight; axis square; grid on;
  set(gca, 'XTick', fig.tick.t2);
  set(gca, 'YTick', fig.tick.t2);
  set(gca, 'XLim', minmax(fig.tick.t2));
  set(gca, 'YLim', minmax(fig.tick.t2));
  set(gca, 'XTickLabel', fig.ticklab.t2);
  set(gca, 'YTickLabel', fig.ticklab.t2);
  set(gca, 'Layer', 'top');
  xlabel('NIST T2 Estimates (ms)', 'FontSize', 16);
  xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
  ylabel('Candidate Profile T2 ML Estimates (ms)', 'FontSize', 16);
  tmp = [t2{2:end}];
  tmp = legend([t2{1} tmp(1,:)],...
    'Ideal', '(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)', 'Location', 'NW');
  set(tmp, 'FontSize', 18);
  if bool.pr
    print('-depsc', '../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t2-ml-compare.eps');
  end
  
  % fig 1: t1,rls
  figure; hold on;
  fill(roi.t1.x.b, roi.t1.y.b, 'y');
  fill(roi.t1.x.t, roi.t1.y.t, color.o);
  t1{1} = plot(log10(1:10000), log10(1:10000), 'k--', 'LineWidth', 1);
  for p = 1:n.p
    t1{p+1} = errorbarxy(...
      stat.logt1.nist.mean, col(stat.logt1.rls.mean(p,:)),...
      stat.logt1.nist.std, col(stat.logt1.rls.std(p,:)),...
      'Color', fig.c{p},...
      'LineStyle', 'none',...
      'Marker', fig.m{p},...
      'LineWidth', 1.0,...
      'MarkerSize', 8,...
      'MarkerFaceColor', fig.f{p});
    t1{p+1} = col(t1{p+1});
  end
  roi.txton = 1:n.roi;
  text(stat.logt1.nist.mean(roi.txton), col(stat.logt1.rls.mean(1,roi.txton)),...
    num2str(roi.txton'),...
    'VerticalAlignment', 'top',...
    'HorizontalAlignment', 'left',...
    'FontSize', 16);
  hold off; axis tight; axis square; grid on;
  set(gca, 'XTick', fig.tick.t1);
  set(gca, 'YTick', fig.tick.t1);
  set(gca, 'XLim', minmax(fig.tick.t1));
  set(gca, 'YLim', minmax(fig.tick.t1));
  set(gca, 'XTickLabel', fig.ticklab.t1);
  set(gca, 'YTickLabel', fig.ticklab.t1);
  set(gca, 'Layer', 'top');
  xlabel('NIST T1 Estimates (ms)', 'FontSize', 16);
  xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
  ylabel('Candidate Profile T1 RLS Estimates (ms)', 'FontSize', 16);
  tmp = [t1{2:end}];
  tmp = legend([t1{1} tmp(1,:)],...
    'Ideal', '(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)', 'Location', 'NW');
  set(tmp, 'FontSize', 18);
  if bool.pr
    print('-depsc', '../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t1-rls-compare.eps');
  end
  
  % fig 2: t2,ml
  figure; hold on;
  fill(roi.t2.x.b, roi.t2.y.b, 'y');
  fill(roi.t2.x.t, roi.t2.y.t, color.o);
  t2{1} = plot(log10(1:1000), log10(1:1000), 'k--', 'LineWidth', 1);
  for p = 1:n.p
    t2{p+1} = errorbarxy(...
      stat.logt2.nist.mean, col(stat.logt2.rls.mean(p,:)),...
      stat.logt2.nist.std, col(stat.logt2.rls.std(p,:)),...
      'Color', fig.c{p},...
      'LineStyle', 'none',...
      'Marker', fig.m{p},...
      'LineWidth', 1.0,...
      'MarkerSize', 8,...
      'MarkerFaceColor', fig.f{p});
    t2{p+1} = col(t2{p+1});
  end
  roi.txton = 1:n.roi;
  text(stat.logt2.nist.mean(roi.txton), col(stat.logt2.rls.mean(1,roi.txton)),...
    num2str(roi.txton'),...
    'VerticalAlignment', 'top',...
    'HorizontalAlignment', 'left',...
    'FontSize', 16);
  hold off; axis tight; axis square; grid on;
  set(gca, 'XTick', fig.tick.t2);
  set(gca, 'YTick', fig.tick.t2);
  set(gca, 'XLim', minmax(fig.tick.t2));
  set(gca, 'YLim', minmax(fig.tick.t2));
  set(gca, 'XTickLabel', fig.ticklab.t2);
  set(gca, 'YTickLabel', fig.ticklab.t2);
  set(gca, 'Layer', 'top');
  xlabel('NIST T2 Estimates (ms)', 'FontSize', 16);
  xlhand = get(gca, 'xlabel'); set(gca, 'FontSize', 16);
  ylabel('Candidate Profile T2 RLS Estimates (ms)', 'FontSize', 16);
  tmp = [t2{2:end}];
  tmp = legend([t2{1} tmp(1,:)],...
    'Ideal', '(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)', 'Location', 'NW');
  set(tmp, 'FontSize', 18);
  if bool.pr
    print('-depsc', '../../data/HPD_Phantom_06,20,16/2016-06-20,hpd,t2-rls-compare.eps');
  end
end

% signal models for ml estimates
if bool.sig
  stat.m0.ml.mean = NaN(n.p, n.roi);
  stat.m0.nist.mean = NaN(n.p, n.roi);
  
  sig.ir.n      = 300;
  sig.ir.tr     = ones(sig.ir.n,1) * 1400;
  sig.ir.ti     = col(linspace(1, max(P.ir.tr(:,1)), sig.ir.n));
  sig.ir.te     = ones(sig.ir.n,1) * 14;
  sig.ir.model  = NaN(sig.ir.n,n.p,n.roi);
  sig.ir.nist   = NaN(sig.ir.n,n.p,n.roi);
  
  sig.se.n      = 300;
  sig.se.tr     = ones(sig.se.n,1) * 1000;
  sig.se.te     = col(linspace(1, 200, sig.se.n));
  sig.se.model  = NaN(sig.se.n,n.p,n.roi);
  sig.se.nist   = NaN(sig.se.n,n.p,n.roi);
        
  for p = 4:n.p
    for r = 6
      stat.m0.ml.mean(p,r) = mean(masker(x(p).m0.ml, roi.mask(:,:,r)));
      
      % m0 ml estimate using nist t1/t2
      tmp = NaN(C.ir+C.se,1);
      for c = 1:C.ir
        tmp2 = IR_gradx_v2(...
          1,...
          stat.t1.nist.mean(r),...
          stat.t2.nist.mean(r),...
          P.ir.tr(c,1), P.ir.ti(c,1), P.ir.te(c,1),...
          'inveff', mean(masker(x(p).inveff.ml, roi.mask(:,:,r))),...
          'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
        tmp(c) = tmp2(1);
      end
      for c = 1:C.se
        tmp2 = SE_gradx(...
          1,...
          stat.t1.nist.mean(r),...
          stat.t2.nist.mean(r),...
          P.se.tr(c,1), P.se.te(c,1),...
          'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
        tmp(c+C.ir) = tmp2(1);
      end
      tmp2 = col(mean(masker(cat(3, ycc.ir, ycc.se), roi.mask(:,:,r)), 1));
      tmp3 = diag([wght(p).ir(:,1); wght(p).se(:,1)]);
      stat.m0.nist.mean(p,r) = (tmp'*tmp3*tmp2) / (tmp'*tmp3*tmp);
        
      % se-ir signals
      if sum(wght(p).ir(:,1))>0
        for i = 1:sig.ir.n
          sig.ir.model(i,p,r) = IR_fun_v4(...
            stat.m0.ml.mean(p,r),...
            stat.t1.ml.mean(p,r),...
            stat.t2.ml.mean(p,r),...
            sig.ir.tr(i), sig.ir.ti(i), sig.ir.te(i),...
            'inveff', mean(masker(x(p).inveff.ml, roi.mask(:,:,r))),...
            'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
          sig.ir.nist(i,p,r) = IR_fun_v4(...
            stat.m0.nist.mean(p,r),...
            stat.t1.nist.mean(r),...
            stat.t2.nist.mean(r),...
            sig.ir.tr(i), sig.ir.ti(i), sig.ir.te(i),...
            'inveff', mean(masker(x(p).inveff.ml, roi.mask(:,:,r))),...
            'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
        end
        
        figure; hold on;...
          plot(sig.ir.ti, real(sig.ir.model(:,p,r)), 'r--');...
          plot(sig.ir.ti, imag(sig.ir.model(:,p,r)), 'b--');...
          plot(sig.ir.ti, abs(sig.ir.model(:,p,r)), 'k--');...
          plot(sig.ir.ti, real(sig.ir.nist(:,p,r)), 'r');...
          plot(sig.ir.ti, imag(sig.ir.nist(:,p,r)), 'b');...
          plot(sig.ir.ti, abs(sig.ir.nist(:,p,r)), 'k');
        for c = 1:C.ir
          if wght(p).ir(c,1)>0
            tmp = mean(masker(ycc.ir(:,:,c), roi.mask(:,:,r)));
            scatter(P.ir.ti(c), real(tmp), 'r*');
            scatter(P.ir.ti(c), imag(tmp), 'b*');
            scatter(P.ir.ti(c), abs(tmp), 'k*');
          end
        end
        hold off;
        title(sprintf('IR: Vial #%u', r));
        xlabel('TI (ms)');
        legend(...
          'model-ml-re', 'model-ml-im', 'model-ml-abs',...
          'model-nist-re', 'model-nist-im', 'model-nist-abs',...
          'data-re', 'data-im', 'data-abs', 'Location', 'Best');
      end
      
      % se signals
      if sum(wght(p).se(:,1))>0
        for i = 1:sig.se.n
          sig.se.model(i,p,r) = SE_fun_v4(...
            stat.m0.ml.mean(p,r),...
            stat.t1.ml.mean(p,r),...
            stat.t2.ml.mean(p,r),...
            sig.se.tr(i), sig.se.te(i),...
            'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
          sig.se.nist(i,p,r) = SE_fun_v4(...
            stat.m0.ml.mean(p,r),...
            stat.t1.nist.mean(r),...
            stat.t2.nist.mean(r),...
            sig.se.tr(i), sig.se.te(i),...
            'kap', mean(masker(double(nu.kap.rls), roi.mask(:,:,r))));
        end
        
        figure; hold on;...
          plot(sig.se.te, real(sig.se.model(:,p,r)), 'r--');...
          plot(sig.se.te, imag(sig.se.model(:,p,r)), 'b--');...
          plot(sig.se.te, abs(sig.se.model(:,p,r)), 'k--');...
          plot(sig.se.te, real(sig.se.nist(:,p,r)), 'r');...
          plot(sig.se.te, imag(sig.se.nist(:,p,r)), 'b');...
          plot(sig.se.te, abs(sig.se.nist(:,p,r)), 'k');
        for c = 1:C.se
          if wght(p).se(c,1)>0
            tmp = mean(masker(ycc.se(:,:,c), roi.mask(:,:,r)));
            scatter(P.se.te(c), real(tmp), 'r*');
            scatter(P.se.te(c), imag(tmp), 'b*');
            scatter(P.se.te(c), abs(tmp), 'k*');
          end
        end
        hold off;
        title(sprintf('SE: Vial %u', r));
        xlabel('TE (ms)');
        legend(...
          'model-ml-re', 'model-ml-im', 'model-ml-abs',...
          'model-nist-re', 'model-nist-im', 'model-nist-abs',...
          'data-re', 'data-im', 'data-abs', 'Location', 'Best');
      end
    end
  end
end