% script brain_2016_05_31_irse_se.m
% joint recon of irse and se data collected 2016-05-31
%
% written by: gopal nataraj
% copyright 2016, university of michigan
%
% version control
%   v0.1    2016-05-26    protocol prepared
%   v0.2    2016-05-31    data acquired
%   v1.1    2016-06-07    original

% irt setup
if ~exist('irtdir', 'var')
  curdir = pwd; 
  cd ../../../irt; 
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../../data/Brain_05,31,16');
addpath('../../map/b1');
addpath('../../map/sense');
addpath('../../map/t1-t2');
addpath('../../model/ir');
addpath('../../model/se');
addpath('../../etc');

% constant declarations
n.x = 256;
n.y = 256;
n.z = 1;
n.c = 32;

fov.x = 240;                            % mm
fov.y = 240;                            % mm
fov.z = 5;                              % mm

% variable acquisition parameters
P.ir.ti   = [50 150 450 1350]';         % ms
P.se.te   = [10 30 60 150]';            % ms

% fixed acquisition parameters
D.ir      = length(P.ir.ti);
P.ir.tr   = ones(D.ir,1) * 1400;        % ms
P.ir.te   = ones(D.ir,1) * 10;          % ms
P.ir.aex  = ones(D.ir,1) * pi/2;        % rad
P.ir.aref = ones(D.ir,1) * pi;          % rad

D.se      = length(P.se.te);
P.se.tr   = ones(D.se,1) * 1000;        % ms
P.se.aex  = ones(D.se,1) * pi/2;        % rad
P.se.aref = ones(D.se,1) * pi;          % rad

% mask threshold
mask.thresh = 0.05;

% recon options
wght.ir = [1 1 1 1]';                   % inversion recovery data weights
wght.se = [1 1 1 1]';                   % spin echo data weights
stop.iter = 100;                        % iteration stop criterion
bool.sv = true;                         % save data files
bool.pr = false;                        % print images
bool.chat = true;                       % verbosity
bool.disp = false;                      % show image updates
bool.precon = true;                     % use preconditioner for rls estimation

% pfile extraction
try
  load('im,irse-se,coil-comb.mat');
catch
  addpath('/Volumes/a2/gnataraj/Data/2016,05,31_brain_32ch/pfile');
  addpath('/Volumes/a2/gnataraj/Data/2016,05,31_brain_32ch/matlab');
  
  % pfile details
  p.fr = n.x;                           % readout length
  p.ep = 0;                             % extended precision 
  p.hosp = 0;                           % hospital data
  p.basetr = 1;                         % baseline tr used
  p.coil = 1:n.c;                       % which coils to extract
  p.sl = 1;                             % which slices to extract
  p.echo = 1;                           % which echoes to extract
  
  % inversion-recovery coil data
  y.k.ir  = NaN(n.x, n.y, n.c, D.ir);
  y.im.ir = NaN(n.x, n.y, n.c, D.ir);
  for i = 1:D.ir
    p.fn = sprintf('P,ir,ti-%dms.7', P.ir.ti(i));
    [y.k.ir(:,:,:,i),~,~] = rawload_jfn(...
      p.fn, p.fr, p.ep, p.hosp, p.basetr, p.coil, p.sl, p.echo);
    for c = 1:n.c
      % shift field-of-view in phase-encoding direction
      y.k.ir(:,2:2:end,c,i) = -y.k.ir(:,2:2:end,c,i);
      
      % images are 2d ifft of data
      y.im.ir(:,:,c,i) = fftshift(ifft2(ifftshift(y.k.ir(:,:,c,i))));
    end
  end
  
  % spin echo coil data
  y.k.se  = NaN(n.x, n.y, n.c, D.se);
  y.im.se = NaN(n.x, n.y, n.c, D.se);
  for i = 1:D.se
    p.fn = sprintf('P,se,te-%dms.7', P.se.te(i));
    [y.k.se(:,:,:,i),~,~] = rawload_jfn(...
      p.fn, p.fr, p.ep, p.hosp, p.basetr, p.coil, p.sl, p.echo);
    for c = 1:n.c
      % shift field-of-view in phase-encoding direction
      y.k.se(:,2:2:end,c,i) = -y.k.se(:,2:2:end,c,i);
      
      % images are 2d ifft of data
      y.im.se(:,:,c,i) = fftshift(ifft2(ifftshift(y.k.se(:,:,c,i))));
    end
  end
    
  % view anterior -> posterior
  y.im.ir = permute(flip(y.im.ir, 2), [2 1 4 3]);
  y.im.se = permute(flip(y.im.se, 2), [2 1 4 3]);
  
  % coil-combine ir and se data together
  [tmp, smap] = mri_multidata_coil_combine(cat(3, y.im.ir, y.im.se));
  ycc.ir = tmp(:,:,1:D.ir);               % [nx ny D.ir]
  ycc.se = tmp(:,:,D.ir+1:end);           % [nx ny D.se]
  
  % save coil-combined images 
  if bool.sv
    curdir = pwd;
    cd('../../data/Brain_05,31,16');
    save('im,irse-se,coil-comb.mat', 'ycc', 'smap'); 
    cd(curdir);
  end
end

% create masks from short te image
tmp = squeeze(abs(ycc.se(:,:,1)));
mask.t = imfill(tmp > mask.thresh * max(col(tmp)), 'holes');
mask.t = imdilate(~imdilate(~mask.t, strel('disk', 5)), strel('disk', 5));
mask.b = imdilate(mask.t, strel('disk', 10));

% parameter estimation
mapArg = {...
  'mask.disp', mask.t,...
  'mask.est', mask.b,...
  'wght.ir', wght.ir,...
  'wght.se', wght.se,...
  'stop.iter', stop.iter,...
  'bool.chat', bool.chat,...
  'bool.disp', bool.disp,...
  'bool.precon', bool.precon};
[x, t] = mri_irse_se_m0t1t2_map(ycc, P, mapArg{:});

% images
figure; im(cat(3, x.m0.ml, x.m0.rls), [0 8], 'cbar');
figure; im(cat(3, x.t1.ml, x.t1.rls), [500 1500], 'cbar');
figure; im(cat(3, x.t2.ml, x.t2.rls), [20 120], 'cbar');
