function V = datasynth(nvoxels,scheme,model,varargin)
% V = datasynth(nv,scheme)
% V = datasynth(nv,scheme,model,...,name,value,...)
% INPUT:
%   nv:      number of voxels
%   scheme: cf camino scheme file format, specifies
%           diffusion acq. params: time, directions
%           contains N diffusion directions
%            --can be a string (filename)
%            --a java camino scheme object
%            --a matlab struct
%   model:  'brownian' (only one implemented currently)
%           'mexp'     (multiexp)
%           'tensor', 'dti', 'exp', 'standard'
%           (others?...)
% 
%   options: datasynth without input lists. See EXAMPLE.
% 
% OUTPUT:
%   V:      a N x nv array containing samples of the data produced by the simulations.
% 
% EXAMPLE:
% V = datasynth(1,'../test/bmx7.scheme','brownian','N_walkers',10,'SNR',16);
%
% philip.batchelor@kcl.ac.uk

%% TODO:
% 
% --multiple models:
%    NB: imaging parameters are common for all models
%        simulation parameters depend on model.
% --improve options processing 
%
%% V = mexp(...);
% either b is set, or compute b from deltas?
% 
% Multi-exp data: cell array of n1 x n2 x n3 x 3 x 3 tensors? (ncomps)
%            or   ND   array of n1 x n2 x n3 x 3 x 3 x ncomps?
%                 vfs (not necessary if ncomps = 1, make defaults?)


%-------------------------------------------------------
% INITIALISATIONS
%-------------------------------------------------------

import simulation.DiffusionSimulation;
import simulation.GeometryFactory;
import simulation.SimulationParams;
import simulation.StepGeneratorFactory;


%-------------------------------------------------------
% OPTIONAL INPUTS
%-------------------------------------------------------

defaults = {
    'N_walkers'      ,  10000;  % number of walkers
    'tmax'           ,  100000; % number of timesteps 
    'p'              ,  1E-4;   % membrane transition probability
    'G'              ,  0.022   % mT/m?
    'delta'          ,  0.032   % s?
    'DELTA'          ,  0.04    % s?
    'DIFF_CONST'     ,  2.02E-9;% m^2 s^{-1}
    'initial'        ,  SimulationParams.UNIFORM;
    'geomType'       ,  GeometryFactory.CELL_STRIPED;
    'voxelSize'      ,  10.0;   % size of the central voxel
    'stepType'       ,  StepGeneratorFactory.FIXEDLENGTH;
    'L'              ,  20.0;   % the substrate size
    'l'              ,  1.0;
    'stripethickness',  3;      % stripe width for striped substrate
    'p_perc'         ,  0.5;    % percolation prob
    'fixedFrac'      ,  0.75;   % fill fraction 
    'modFill'        ,  4;      % 
    'modFree'        ,  1;      %
    'SNR'            ,  -1;     % 
    'modQ'           ,  -1;     % fixed |q|
    'qScale'         ,  1.0;    % 
    'tau'            ,  1.0;    %				% 
    'tauScale'       ,  1.0;    % scale for diffusion time
    'M'              ,  -1;     % number of b=0 images
    'N'              ,  -1;     % number of b>0 images
    'sequenceIndex'  ,  0; 
    'diffusionTime'  ,  -1;     %
    'seed'           , 	36558013;
	   };

if nargin == 0; defaults, return; end;
if nargin < 3 || isempty(model); model = 'brownian'; end;

simp = parseargs(cell2struct(defaults(:,2),defaults(:,1),1),varargin{:});
qScale        = simp.qScale;
modQ          = simp.modQ;
tauScale      = simp.tauScale;
diffusionTime = simp.diffusionTime;
sequenceIndex = simp.sequenceIndex;

%------------------------------------------
% Arguments passed to CL_Initializer
%------------------------------------------
S = javaArray('java.lang.String',8);
S(1) = java.lang.String('-snr');             S(2) = java.lang.String(num2str(simp.SNR));
S(3) = java.lang.String('-seed');            S(4) = java.lang.String(num2str(simp.seed));
S(5) = java.lang.String('-stripethickness'); S(6) = java.lang.String(num2str(simp.stripethickness));
S(7) = java.lang.String('-latticesize');     S(8) = java.lang.String(num2str(simp.L));
tools.CL_Initializer.CL_init(S);
% NB: could pass all as strings, but seems cleaner like this
% ideally, we would want to not use CL_Initializer for 
% Matlab?

%-------------------------------------------------------
% Imaging Parameters: Read Scheme
%-------------------------------------------------------

if ischar(scheme)
  if exist(scheme,'file') == 2 % use camino?
    [simp.DELTA,simp.N,simp.M,simp.g] = readscheme(scheme);
    schemeVersion = tools.CL_Initializer.getSchemeVersion(scheme);
    if(schemeVersion == 0); imPars = imaging.SchemeV0(scheme, qScale, tauScale);
    else;                   imPars = imaging.SchemeV1(scheme);end
  end
elseif isstruct(scheme)
  simp.DELTA = scheme.DELTA; simp.N = scheme.N; simp.M = scheme.M; simp.g = scheme.g;
  if (simp.M >= 0) && (simp.N > 0)
    imPars = imaging.SchemeV0(simp.M, simp.N, qScale*modQ, tauScale*diffusionTime);
  end
elseif isjava(scheme) && (isa('SchemeV0',scheme) || isa('SchemeV1',scheme))
  imPars = scheme; % check correct class?
else
  imPars = SchemeV0.getSchemeV0(sequenceIndex, qScale, tauScale);
end;

M = simp.M; % Number of b  = 0 images
N = simp.N; % Number of b ~= 0 images

%-------------------------------------------------------
% MODELS
%-------------------------------------------------------

switch model,
 case 'brownian', 
  
 V = browniansimulation(nvoxels,M,N,imPars,simp);
 case 'exponential', %TODO

 otherwise % standard single-exp
end



%=======================================================
function V = browniansimulation(nvoxels,M,N,imPars,simp)
% Run the 'proper' simulation to generate the data.

%-------------------------------------------------------
% Brownian Motion Parameters (Simulation)
%-------------------------------------------------------
geomParams = tools.CL_Initializer.getGeometryParams(simp.geomType);
simParams  = simulation.SimulationParams(simp.N_walkers,simp.tmax,simp.p,simp.initial, simp.geomType, geomParams,simp.stepType, simp.voxelSize, imPars);
stepParams = simulation.StepGeneratorFactory.getStepParamsArray(simp.stepType,simParams, imPars);
simParams.setStepParams(stepParams);

%-------------------------------------------------------
% Run the Simulation
%-------------------------------------------------------

data = simulation.DiffusionSimulation(simParams, imPars);     

%-------------------------------------------------------
% Convert Output to Matlab
%-------------------------------------------------------

V = zeros(N+M,nvoxels);
for n = 1:nvoxels; V(:,n) = data.nextVoxel(); end
return;

%=======================================================
function V = mexp(M,N,simp)
% Simulate multi-exponential data.
% 
% 
% Allow for D_xy parameters
% 
% Allow for eigenvector-eigenvalue parameters?
% (i.e. singular tensors too)
V = []; warning('TODO');
return;

%=======================================================
function V = host(M,N,simp)
% Simulate higher order symmetric tensor data.
% 
% Allow for 'directions' input?
% 
% Allow for 'spherical harmonics'?
V = []; warning('TODO');
return;
%=======================================================
function [imPars,simp] = schemeread(scheme,simp)
qScale        = simp.qScale;
modQ          = simp.modQ;
tauScale      = simp.tauScale;
diffusionTime = simp.diffusionTime;
sequenceIndex = simp.sequenceIndex;

if ischar(scheme)
  if exist(scheme,'file') == 2 % use camino?
    [simp.DELTA,simp.N,simp.M,simp.g] = readscheme(scheme);
    schemeVersion = tools.CL_Initializer.getSchemeVersion(scheme),
    if(schemeVersion == 0) 
      imPars = imaging.SchemeV0(scheme, qScale, tauScale);
    else
      imPars = imaging.SchemeV1(scheme);
    end
  end
elseif isstruct(scheme)
  simp.DELTA = scheme.DELTA; 
  simp.N = scheme.N; 
  simp.M = scheme.M; 
  simp.g = scheme.g;
  if (simp.M >= 0) & (simp.N > 0)
    imPars = imaging.SchemeV0(M, N, qScale * modQ, tauScale * diffusionTime);
  end
elseif isjava(scheme) & (isa('SchemeV0',scheme) | isa('SchemeV1',scheme))
  imPars = scheme; % check correct class?
else
  imPars = SchemeV0.getSchemeV0(sequenceIndex, qScale, tauScale);
end;
