% script HPDphant_2016_06_17.m
% joint recon of bs followed by all ir-se/se/spgr/dess scan profiles
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2016-05-26    protocol prepared
%   v0.2    2016-06-17    data acquired
%   v1.1    2016-06-18    original: bad ir/se data; spgr/dess coil-combine failing

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/HPD_Phantom_06,17,16');
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
n.sl = 5;
n.p = 4;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

% acquisition parameters
P.ir.ti   = [50 150 450 1350]';         % ms
C.ir      = length(P.ir.ti);
E.ir      = 1;
P.ir.tr   = ones(C.ir,1) * 1400;        % ms
P.ir.te   = ones(C.ir,E.ir) * 10;       % ms
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
bool.rs   = false;                       % roi manual selection
bool.im   = false;                       % show images
bool.roi  = true;                       % show rois in images
bool.pr   = false;                       % print images and statistics

% mask threshold
mask.thresh = 0.05;

% recon options
dict.t1.n = 200;                        % num t1 dict elements
dict.t2.n = 200;                        % num t2 dict elements
stop.iter = 10;                        % iteration stop criterion
bool.mag.sp = true;                     % use magnitude spgr data/model
bool.mag.de = true;                     % use magnitude dess data/model
bool.chat = true;                       % verbosity
bool.disp = false;                      % show image updates
bool.precon = true;                     % use preconditioner for rls estimation

% bloch-siegert b1 mapping
try 
  fprintf('\nTrying to load b1/kap maps...');
  load(sprintf('im,nu,sl-%u.mat', n.sl));
  fprintf('success!\n');
catch
  fprintf('failed: will reconstruct.\n'); 
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/matlab');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/wav');
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
    permute(flip(flip(y.im.bs.p(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 4 3]);                 
  y.im.bs.m = ...                       % [n.x n.y n.c]
    permute(flip(flip(y.im.bs.m(f.cropx,:,n.sl,:,:,:), 1), 2), [2 1 4 3]);

  % bs b1 mapping
  [~,P.bs.mag,~] = readwav('bsp.wav');  % G   
  P.bs.wrf  = 8;                        % kHz
  P.bs.dt   = 4e-3;                     % ms
  b1mapArg = {...
    'coilOpt', {'nouter', 10, 'thresh', 0.1},...
    'reg.beta', 2^2,...
    'scale', 1,...
    'bool.chat', true,...
    'bool.disp', true};
  nu = mri_bs_kapb1_map(y.im.bs, P.bs, b1mapArg{:});
  
  % save b1/kap maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,17,16');
    save(sprintf('im,nu,sl-%u.mat', n.sl), 'nu');
    cd(curdir);
  end
end

% optional: uniform b1/kap scaling
if bool.sc
  scale = 1.016/0.913;
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
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/matlab');
  
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
    permute(flip(y.im.ir, 2), [2 1 4 5 3]);
  y.im.se = ...                         % [n.x n.y C.se E.se n.c]
    permute(flip(y.im.se, 2), [2 1 4 5 3]);
  
  % coil-combine ir and se data together
  tmp = cat(3,...                       % [n.x n.y C.ir*E.ir+C.se*E.se n.c]
    reshape(y.im.ir, [n.x n.y C.ir*E.ir n.c]),...
    reshape(y.im.se, [n.x n.y C.se*E.se n.c]));
  [tmp, smap] = mri_multidata_coil_combine(tmp);
  ycc.ir = reshape(tmp(:,:,1:(C.ir*E.ir)), [n.x n.y C.ir E.ir]);
  ycc.se = reshape(tmp(:,:,(C.ir*E.ir+1):end), [n.x n.y C.se E.se]);
  
  % save coil-combined images 
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,17,16');
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
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/matlab');
  addpath('/Volumes/a2/gnataraj/Data/2016,06,17_hpd_8ch/wav');
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
  
  % coil-combine spgr and dess data together
  tmp = cat(3,...                       % [n.x n.y C.sp*E.sp+C.de*E.de n.c]
    reshape(y.im.sp, [n.x n.y C.sp*E.sp n.c]),...
    reshape(y.im.de, [n.x n.y C.de*E.de n.c]));
%   [tmp, smap] = mri_multidata_coil_combine(tmp);
  [tmp, smap] = mri_multidata_coil_combine(tmp, 'nouter', 0);
  ycc.sp = reshape(tmp(:,:,1:(C.sp*E.sp)), [n.x n.y C.sp E.sp]);
  ycc.de = reshape(tmp(:,:,(C.sp*E.sp+1):end), [n.x n.y C.de E.de]);
  
  % save coil-combined images
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_06,17,16');
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
if bool.mag.sp
  ycc.sp = abs(ycc.sp);
end
if bool.mag.de
  ycc.de = abs(ycc.de);
end

% test different scan profiles
for p = 1:n.p
  % set data weights based on profile
  switch p
    case 1                              % (2 sp, 1 de)
      wght.ir = [0 0 0 0]'  * ones(1,E.ir);                   
      wght.se = [0 0 0 0]'  * ones(1,E.se);   
      wght.sp = [1 1 0]'    * ones(1,E.sp);
      wght.de = [0 1 0 0]'  * ones(1,E.de);
    case 2                              % (1 sp, 1 de)
      wght.ir = [0 0 0 0]'  * ones(1,E.ir);                   
      wght.se = [0 0 0 0]'  * ones(1,E.se);   
      wght.sp = [0 0 1]'    * ones(1,E.sp);
      wght.de = [0 0 1 0]'  * ones(1,E.de);
    case 3                              % (0 sp, 2 de)
      wght.ir = [0 0 0 0]'  * ones(1,E.ir);                   
      wght.se = [0 0 0 0]'  * ones(1,E.se);   
      wght.sp = [0 0 0]'    * ones(1,E.sp);
      wght.de = [1 0 0 1]'  * ones(1,E.de);
    case 4                              % (4 ir, 4 se)
      wght.ir = [1 1 1 1]'  * ones(1,E.ir);     
      wght.se = [1 1 1 1]'  * ones(1,E.se);   
      wght.sp = [0 0 0]'    * ones(1,E.sp);
      wght.de = [0 0 0 0]'  * ones(1,E.de);
  end
  fprintf('\nTesting profile %d: (%u,%u,%u,%u) ir/se/sp/de scans.\n',...
    p, sum(wght.ir(:,1)), sum(wght.se(:,1)), sum(wght.sp(:,1)), sum(wght.de(:,1)));
  
  % parameter estimation
  mapArg = {...
    'mask.disp', mask.t,...
    'mask.est', mask.b,...
    'nu.kap', double(nu.kap.rls),...
    'wght.ir', wght.ir,...
    'wght.se', wght.se,...
    'wght.sp', wght.sp,...
    'wght.de', wght.de,...
    'dict.t1.n', dict.t1.n,...
    'dict.t2.n', dict.t2.n,...
    'stop.iter', stop.iter,...
    'bool.mag.sp', bool.mag.sp,...
    'bool.mag.de', bool.mag.de,...
    'bool.chat', bool.chat,...
    'bool.disp', bool.disp,...
    'bool.precon', bool.precon};
  [x(p), t(p)] = mri_m0t1t2_map(ycc, P, mapArg{:});
end