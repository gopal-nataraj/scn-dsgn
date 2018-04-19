% script HPDphant_2016_07_25.m
% comparison of b1 maps from bloch-siegert vs double-angle
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v1.0    2016-07-22    protocol prepared; data acquisition
%   v1.1    2016-07-23    recon original; discovered aliasing artifacts
%   v2.0    2016-07-25    protocol revised; new data acquisition; data clipped
%   v2.1    2016-07-26    protocol revised; new data acquisition; data clipped
%   v2.2    2016-07-27    protocol repeated with lower gains, bs separate

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../../irt; 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/HPD_Phantom_07,27,16');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../etc');

% raw data path
% datdir = '/Volumes/a2/gnataraj/Data/2016,07,27_hpd_8ch';
datdir = '/Volumes/General Storage/Documents/school/michigan/research/data/2016,07,27_hpd_8ch';

% constant declarations 
n.x = 256;
n.y = 256;
n.z = 6;
n.c = 8;
n.sl = 3;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

% acquisition parameters
P.sp.aex  = [90 45]' * pi/180;          % rad
P.sp.tr   = [4000 4000]';               % ms
C.sp      = length(P.sp.aex);           
E.sp      = 1;
P.sp.te   = ones(C.sp,E.sp) * 4.67;     

% mask threshold
mask.thresh = 0.05;

% coil combination options
opt.coil = {...
  'log2b', 10,...
  'fwhm', 1,...
  'nouter', 10,...
  'relPh', 1,...
  'disp', 0};

% da b1 mapping options
opt.b1.da = {...
  'tips', P.sp.aex,...
  'l2b_bmap', 40,...
  'chat', 1,...
  'niter', 5000,...
  'isave', [0; 5000]};

% bs b1 mapping options
opt.b1.bs = {...
  'coilOpt', opt.coil,...
  'reg.beta', 2^2,...
  'scale', 1,...
  'stop.iter', 2000,...
  'bool.chat', true,...
  'bool.disp', true};

% output options
bool.sv = 1;
bool.im = 1;

% b1 mapping
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
  
  % da pfile details
  f.name = 'P,spgr-bs.7';
  f.id = fopen(f.name, 'r', 'l'); 
    fseek(f.id,0,'bof');
  f.ver = str2double(num2str(fread(f.id,1,'float32')));
    fseek(f.id,0,'bof');
  f.hdr = read_rdb_hdr(f.id, f.ver);
  f.nframe = 4;
  f.necho = 1;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (512-127:512+128);
  f.yres = n.y;
  
  % double-angle coil data (frames 1-2)
  y.k.sp = NaN(f.rl, n.y, n.z, n.c, C.sp, E.sp);
  y.im.sp = NaN(f.rl, n.y, n.z, n.c, C.sp, E.sp);
  for i = 1:C.sp
    for j = 1:E.sp
      [y.im.sp(:,:,:,:,i,j), y.k.sp(:,:,:,:,i,j)] = recon_gn(...
        f.name, i, j, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
    end
  end
  
  % bs pfile details
  f.name = 'P,bs.7';
  f.id = fopen(f.name, 'r', 'l'); 
    fseek(f.id,0,'bof');
  f.ver = str2double(num2str(fread(f.id,1,'float32')));
    fseek(f.id,0,'bof');
  f.hdr = read_rdb_hdr(f.id, f.ver);
  f.nframe = 2;
  f.necho = 1;
  f.disp = false;
  f.fft3 = true;
  f.flip = 1;
  f.rl = 1024;
  f.cropx = (512-127:512+128);
  f.yres = n.y;
  
  % bloch-siegert coil data (frames 1-2)
  [y.im.bs.p, y.k.bs.p] = recon_gn(...
    f.name, 1, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
  [y.im.bs.m, y.k.bs.m] = recon_gn(...
    f.name, 2, 1, f.necho, f.nframe, f.disp, f.fft3, f.flip, f.yres);
  
  % crop to single image slice and readout portion of interest
  % reorient to view superior -> inferior, right -> left    
  y.im.sp = ...                         % [n.x n.y C.sp E.sp n.c]
    permute(y.im.sp(f.cropx,:,n.sl,:,:,:), [2 1 5 6 4 3]);
  y.im.bs.p = ...                       % [n.x n.y n.c]
    permute(y.im.bs.p(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);                 
  y.im.bs.m = ...                       % [n.x n.y n.c]
    permute(y.im.bs.m(f.cropx,:,n.sl,:,:,:), [2 1 4 3]);
  
  % da b1 mapping
  tmp = reshape(y.im.sp, [n.x n.y C.sp*E.sp n.c]);
  [ycc.sp, smap] = ...                  % [n.x n.y C.sp E.sp]
    mri_multidata_coil_combine(tmp, opt.coil{:});
  tmp = permute(ycc.sp, [1 2 5 3 4]);   % [n.x n.y 1 C.sp E.sp]
  tmp = mri_b1map(tmp, opt.b1.da{:});
  nu.da.kap.mom = tmp(:,:,1);
  nu.da.kap.rls = tmp(:,:,2);
  
  % bs b1 mapping
  [~,P.bs.mag,~] = readwav('bsp.wav');  % G   
  P.bs.wrf  = 8;                        % kHz
  P.bs.dt   = 4e-3;                     % ms
  nu.bs = mri_bs_kapb1_map(y.im.bs, P.bs, opt.b1.bs{:});
  
  % save b1/kap maps
  if bool.sv
    curdir = pwd;
    cd('../../data/HPD_Phantom_07,27,16');
    save(sprintf('im,nu,sl-%u.mat', n.sl), 'nu');
    cd(curdir);
  end
end

% image
if bool.im
  figure, im(cat(3, nu.da.kap.mom, nu.da.kap.rls, nu.bs.kap.mom, nu.bs.kap.rls), [0 1.2], 'cbar');
end