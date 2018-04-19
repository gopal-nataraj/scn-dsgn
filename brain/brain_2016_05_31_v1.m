% script brain_2016_05_31.m
% joint recon of bs followed by all ir-se/se/spgr/dess scan profiles
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2016-05-26    protocol prepared
%   v0.2    2016-05-31    data acquired
%   v1.1    2016-06-07    original (irse/se only)
%   v2.1    2016-06-13    now using mri_m0t1t2_map.m
%   v2.2    2016-06-16    test different profiles
%   v2.3    2016-06-30    now using mri_m0t1t2inveff_map.m to avoid se-ir t1 bias
%   v2.4    2016-07-03    coil-combine separately by profile
%   v2.5    2016-07-14    register coil-combined images prior to mapping
%   v2.6    2016-08-02    new b1 rescaling from double-angle experiment

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/Brain_05,31,16');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../model/spgr');
addpath('../../model/dess');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../etc');
addpath('../../../../../contrib/mat/web');

% raw data path
% datdir = '/Volumes/a2/gnataraj/Data/2016,05,31_brain_32ch';
datdir = '/Volumes/General Storage/Documents/school/michigan/research/data/2016,05,31_brain_32ch';

% constant declarations
n.x = 256;
n.y = 256;
n.z = 8;
n.c = 32;
n.sl = 5;
n.p = 5;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

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
bool.sv   = false;                      % save coil-comb, param reconstructions
bool.reg  = true;                       % register coil-comb images
bool.rs   = false;                      % roi manual selection
bool.im   = true;                       % show images
bool.roi  = false;                      % show rois in all images
bool.pr   = false;                      % print images and statistics

% mask threshold
mask.thresh = 0.05;

% coil combination options
opt.coil = {...
  'nouter', 10,...
  'thresh', 0.1,...
  'disp', 0};

opt.b1 = {...
  'coilOpt', opt.coil,...
  'reg.beta', 2^2,...
  'scale', 1,...
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
  addpath([datdir '/pfile']);
  addpath([datdir '/matlab']);
  addpath([datdir '/wav']);
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
  % reorient to view anterior -> posterior, right -> left            
  y.im.bs.p = ...                       % [n.x n.y n.c]
    permute(flip(flip(y.im.bs.p(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 4 3]);                 
  y.im.bs.m = ...                       % [n.x n.y n.c]
    permute(flip(flip(y.im.bs.m(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 4 3]);

  % bs b1 mapping
  [~,P.bs.mag,~] = readwav('bsp.wav');  % G   
  P.bs.wrf  = 8;                        % kHz
  P.bs.dt   = 4e-3;                     % ms
  nu = mri_bs_kapb1_map(y.im.bs, P.bs, opt.b1{:});
  
  % save b1/kap maps
  if bool.sv
    curdir = pwd;
    cd('../../data/Brain_05,31,16');
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
  fprintf('failed.\nProceeding to extract coil-combined raw data.\n');
  addpath([datdir '/pfile']);
  addpath([datdir '/matlab']);
  
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
    
  % reorient to view anterior -> posterior, right -> left
  y.im.ir = ...                         % [n.x n.y C.ir E.ir n.c]
    permute(flip(y.im.ir, 2), [2 1 4 5 3]);
  y.im.se = ...                         % [n.x n.y C.se E.se n.c]
    permute(flip(y.im.se, 2), [2 1 4 5 3]);
  
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
    cd('../../data/Brain_05,31,16');
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
  fprintf('failed.\nProceeding to extract and coil-combined raw data.\n');
  addpath([datdir '/pfile']);
  addpath([datdir '/matlab']);
  addpath([datdir '/wav']);
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
    permute(flip(flip(y.im.sp(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 5 6 4 3]); 
  y.im.de = ...                         % [n.x n.y C.de E.de n.c]
    permute(flip(flip(y.im.de(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 5 6 4 3]);  
  
%   % coil-combine spgr and dess data together
%   tmp = cat(3,...                       % [n.x n.y C.sp*E.sp+C.de*E.de n.c]
%     reshape(y.im.sp, [n.x n.y C.sp*E.sp n.c]),...
%     reshape(y.im.de, [n.x n.y C.de*E.de n.c]));
%   [tmp, smap] = mri_multidata_coil_combine(tmp, opt.coil{:});
%   ycc.sp = reshape(tmp(:,:,1:(C.sp*E.sp)), [n.x n.y C.sp E.sp]);
%   ycc.de = reshape(tmp(:,:,(C.sp*E.sp+1):end), [n.x n.y C.de E.de]);

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
    cd('../../data/Brain_05,31,16');
    save(sprintf('im,spgr-dess,coil-comb,sl-%u.mat', n.sl), 'ycc', 'smap');
    cd(curdir);
  end
  
  % append se/ir coil-combined images
  tmp = ycc;
  clear('ycc');
  load('im,irse-se,coil-comb.mat');
  ycc = catstruct(tmp, ycc);
end

% register coil-combined images and nu.kap.rls to first ycc.ir image
if bool.reg
  fprintf('\nRegistering coil-combined images...\n');
  reg.fixed = abs(ycc.ir(:,:,1));
  reg.ttype = 'rigid';
  [reg.optim, reg.metric] = imregconfig('multimodal');
  
  for i = 1:C.ir
    fprintf('\tRegistering IR dataset %u to reg.fixed...', i);
    tmp = imregtform(abs(ycc.ir(:,:,i,1)), reg.fixed, reg.ttype, reg.optim, reg.metric);
    for j = 1:E.ir
      ycc.ir(:,:,i,j) = imwarp(ycc.ir(:,:,i,j), tmp, 'OutputView', imref2d(size(reg.fixed)));
    end
    fprintf('done.\n');
  end
  for i = 1:C.se
    fprintf('\tRegistering SE dataset %u to reg.fixed...', i);
    tmp = imregtform(abs(ycc.se(:,:,i,1)), reg.fixed, reg.ttype, reg.optim, reg.metric);
    for j = 1:E.se
      ycc.se(:,:,i,j) = imwarp(ycc.se(:,:,i,j), tmp, 'OutputView', imref2d(size(reg.fixed)));
    end
    fprintf('done.\n');
  end
  for i = 1:C.sp
    fprintf('\tRegistering SPGR dataset %u to reg.fixed...', i);
    tmp = imregtform(abs(ycc.sp(:,:,i,1)), reg.fixed, reg.ttype, reg.optim, reg.metric);
    for j = 1:E.sp
      ycc.sp(:,:,i,j) = imwarp(ycc.sp(:,:,i,j), tmp, 'OutputView', imref2d(size(reg.fixed)));
    end
    fprintf('done.\n');
  end
  for i = 1:C.de
    fprintf('\tRegistering DESS dataset %u to reg.fixed...', i);
    tmp = imregtform(abs(ycc.de(:,:,i,1)), reg.fixed, reg.ttype, reg.optim, reg.metric);
    for j = 1:E.de
      ycc.de(:,:,i,j) = imwarp(ycc.de(:,:,i,j), tmp, 'OutputView', imref2d(size(reg.fixed)));
    end
    fprintf('done.\n');
  end
  
  fprintf('\tRegistering kap map to reg.fixed...');
  tmp = imregtform(abs(nu.kap.rls), reg.fixed, reg.ttype, reg.optim, reg.metric);
  nu.kap.rls = imwarp(nu.kap.rls, tmp, 'OutputView', imref2d(size(reg.fixed)));
  fprintf('done.\n');
end

% debug: try using method-of-moments b0map
if 0
  num = (angle(ycc.de(:,:,2,1) .* conj(ycc.de(:,:,2,2))) - pi);
  den = 4*pi*P.de.te(1);
  nu.b0 = div0(num, den);
  figure;
    dyrng = 1/(2*P.de.te(1));
    clim = [-dyrng 0];
    im(nu.b0, clim, 'off-resonance (kHz)', 'cbar'); 
    drawnow;
else
  nu.b0 = zeros(n.x,n.y);
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
  load(sprintf('../../data/Brain_05,31,16/im,x,sl-%u.mat', n.sl));
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
      'nu.kap', ones(size(mask.t)),...double(nu.kap.rls),...
      'nu.b0', nu.b0,...
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
    cd('../../data/Brain_05,31,16');
    save(sprintf('im,x,sl-%u.mat', n.sl), 'x', 't');
    cd(curdir);
  end
end

% roi selection
if bool.rs
  tmp = 'ROI selection: [m]anual or [s]aved? ';
  roi.sel = input(tmp, 's');
else
  roi.sel = 's';
end
roi.label = {'wm'; 'gm'; 'sf'};
n.roi = length(roi.label);

switch roi.sel
  case 'm'
    roi.mask = false(n.x, n.y, n.roi);
    for r = 1:n.roi
      while true
        fprintf('Select %s polygon.\n', roi.label{r});
        tmp = roipoly(x(n.p).t1.rls/2000);
        roi.mask(:,:,r) = roi.mask(:,:,r) | tmp;
        tmp = input('Add another polygon to this ROI [*/n]? ', 's');
        if strcmp(tmp,'n'), break; end;
      end
    end

    % save mask
    dir = pwd; 
    cd('../../data/Brain_05,31,16');
    if exist('roi-mask.mat', 'file')
      tmp = input('Overwrite previously-defined ROI mask [y/n]? ', 's');
      switch tmp
        case 'y'
          save('roi-mask.mat', 'roi'); 
          fprintf('roi-mask.mat overwritten.\n');
        case 'n'
          fprintf('New ROI mask not saved.\n');
        otherwise
          error('Unrecognized input.');
      end
    else
      save('roi-mask.mat', 'roi');
    end
    cd(dir);

  case 's'
    try
      addpath('../../data/Brain_05,31,16');
      load('roi-mask.mat');
    catch
      error('ROI masks not found!');
    end

  otherwise
    error('Unrecognized input!');
end

% boundaries
roi.bound.wm = bwboundaries(roi.mask(:,:,1), 'noholes');
roi.bound.gm = bwboundaries(roi.mask(:,:,2), 'noholes');

% omit profile 4 for analysis
if n.p == 5
  x(5).inveff = x(4).inveff;
  x(4) = [];
  n.p = 4;
end

% images
rng.m0 = [0 n.c];
rng.t1 = [600 1600];
rng.t2 = [20 120];
rng.inveff = [0.5 1];

opt.im.cmap = {'jet', 'gray'};
opt.im.prof = {'(2SP,1DE)', '(1SP,1DE)', '(0SP,2DE)', '(4IR,4SE)'};
opt.im.t1 = {'T1 ML', 'T1 RLS'};
opt.im.t2 = {'T2 ML', 'T2 RLS'};

if bool.im
  for m = 1:length(opt.im.cmap)
    % t1 image
    tmp = [x(:).t1];
    figure; im('notick', 'row', 2, cat(3, tmp(:).ml, tmp(:).rls), rng.t1, 'cbar', ' ');...
      colormap(opt.im.cmap{m});...
      tmp = colorbar; set(tmp, 'YTick', linspace(rng.t1(1), rng.t1(2), 6));...
      hold on;...
      if m==2
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1) + (n.p-1)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
          if bool.roi
            plot(tmp(:,1) + (n.p-2)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
            plot(tmp(:,1) + (n.p-3)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
            plot(tmp(:,1) + (n.p-4)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
          end
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1) + (n.p-1)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
          if bool.roi
            plot(tmp(:,1) + (n.p-2)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
            plot(tmp(:,1) + (n.p-3)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
            plot(tmp(:,1) + (n.p-4)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
          end
        end
      end
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
      hold off;
      if bool.pr
        tmp = sprintf('../../data/Brain_05,31,16/2016-05-31,brain,t1,%s.eps', opt.im.cmap{m});
        print('-depsc', tmp);
      end
      
    % t2 image
    tmp = [x(:).t2];
    figure; im('notick', 'row', 2, cat(3, tmp(:).ml, tmp(:).rls), rng.t2, 'cbar', ' ');...
      colormap(opt.im.cmap{m});...
      tmp = colorbar; set(tmp, 'YTick', linspace(rng.t2(1), rng.t2(2), 6));
      hold on;
      if m==2 && bool.roi
        for b = 1:length(roi.bound.wm)
          tmp = roi.bound.wm{b};
          plot(tmp(:,1) + (n.p-1)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-2)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-3)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-4)*n.x, tmp(:,2), 'm', 'LineWidth', 0.5);
        end
        for b = 1:length(roi.bound.gm)
          tmp = roi.bound.gm{b};
          plot(tmp(:,1) + (n.p-1)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-2)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-3)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
          plot(tmp(:,1) + (n.p-4)*n.x, tmp(:,2), 'c', 'LineWidth', 0.5);
        end
      end
      text(zeros(2,1), col(n.y/2-1:n.y:3*n.y/2), col(opt.im.t2),...
        'VerticalAlignment', 'bottom',...
        'HorizontalAlignment', 'center',...
        'FontSize', 16,...
        'Color', 'k',...
        'Rotation', 90);
      hold off;
      if bool.pr
        tmp = sprintf('../../data/Brain_05,31,16/2016-05-31,brain,t2,%s.eps', opt.im.cmap{m});
        print('-depsc', tmp);
      end
  end
elseif bool.rs
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
  end
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
  fid = fopen('../../data/Brain_05,31,16/2016-05-31,brain,stat', 'w');
  fprintf(fid, '2016-05-31 parameter estimate statistics\n');
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
  
  fprintf(fid, '\n\n2016-05-31 parameter estimate run times (s)\n');
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
