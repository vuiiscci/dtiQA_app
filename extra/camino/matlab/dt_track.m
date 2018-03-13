function trs = dt_track(Ds,seedROI,varargin)
% trs = dt_track(Ds)
% trs = dt_track(Ds,seedROI)
% trs = dt_track(Ds,seedROI,...,[optionstr,value],...)
% 
% Does DT tractography using CAMINO. 
%
% INPUT:
% 
%  Ds: a tensor field, i.e. a 3D n1 x n2 x n3 volume of 3x3 tensors, stored as 
%      a n1 x n2 x n3 x 3 x 3 array.
% 
%  seedROI: the seed mask, a scalar volume of same dimensions as
%      the tensor field, with n at seed position, zero
%      (false) otherwise. The different values of nonzero denote
%      different seed regions. For the mo, every tract goes into
%      one cell array.
% 
%      
%      This is a java 'short', i.e. Matlab 'int16'. seedROI is thus a n1 x n2 x n3 array.
%      If not specified, or [], seeding is from everywhere.
%    
%  option,value pair: see track options in camino, 
%          Strings              Default values         Description
%          ---------------------------------------------------------
%	   'interpolated',      false                  linear or not
%	   'iterations',        5000        
%	   'minTractPoints'     0                      min number of pts
%	   'VoxelDims',         [1,1,1]               
%          'anisThresh',        0                      threshold on FA?
%	   'ipThresh',          cos(pi * 80.0 / 180.0) threshold on curv
%	   'stepSize',          0                      
%          'waypointVol',       []                     
%
%
% OUTPUT:
%  trs: a cell array of fiber tracks, 
%       length(trs) = number of nonzeros in seedROI,
%       each cell is a row array of 3D point coordinates, i.e.
%       trs{k} = [x1 y1 z1
%                 x2 y2 z2
%                 ........]
% 
% 
% % EXAMPLE 1: Synthetic Helicoid
% %-------------------------------
% 
% % Create a helicoid in a 64 x 64 x 64 volume:
% n1 = 64; n2 = 64; n3 = 64
% t = linspace(0,2*pi,1000)';
% c   = [16*cos(t(:)) + 32,16*sin(t(:)) + 32,64/(2*pi)*t(:)] + 1;
% c_t = [-16*sin(t), 16*cos(t), 64/(2*pi)*ones(size(t))];
% rc  = round(c); 
% rc(rc(:,1) > 64 | rc(:,1) < 1,:) = [];
% rc(rc(:,2) > 64 | rc(:,2) < 1,:) = [];
% rc(rc(:,3) > 64 | rc(:,3) < 1,:) = [];
% 
% % This puts a tensor with PD // to c in voxels through which c goes
% D = zeros(64,64,64,3,3);
%  for k = 1:length(rc)-1
%    D(rc(k,1),rc(k,2),rc(k,3),:,:) = c_t(k,:)'*c_t(k,:);
%  end;
% 
% % EXAMPLE 2: (real data)
% %-----------------------
% load ~/projects/data/chris_clarks_data/2006-03-15_123654/tensor;
% addpath('~/projects/DTITOOLS/')
% FA = invariants(D,'fa');
% ROI = zeros(size(FA),'int16');
% ROI(~isnan(FA)) = int16(1); D(isnan(D)) = 0; D = D + 1e-6*randn(size(D));
% for x = 1:96; for y = 1:96; for z = 1:50 
% D(x,y,z,:,:) = sqrtm(squeeze(D(x,y,z,:,:))'*squeeze(D(x,y,z,:,:)));
% end; end; end;
% trs = dt_track(D,ROI);
% 
% -----------------------------------------------------------------------
% philip.batchelor@ucl.ac.uk
% ========================================================================
% WARNING: Uses CAMINO from D.~Alexander et al in the backrgound, requires
% running Matlab with jvm (that's the default), and camino.jar. See SETUP.M
% Adapted from StreamlineTractrography.java from  @author Philip Cook
%
% SEE ALSO: http://www.cs.ucl.ac.uk/research/medic/camino/
% ========================================================================
%
% TODO: --find ROIs by non-rigid reg of known one?
%       --multiple ROIS?
%       --pass anisThresh using this DT_TractographyImage constructor
%       
%-------------------------------------------------------
% Convert Inputs to Java Arrays for passing to Camino
%-------------------------------------------------------
% Assume 1 tensor for the mo:

%-------------------------------------------------------
% CHECK FOR CAMINO.JAR
%-------------------------------------------------------

%%global CAMINO_PATH
%%global CAMINO_TXT
v = version('-java');
if strcmp(v,'Java is not enabled'); error('This code requires Matlab with jvm'); end;

%%if exist('tractography.DT_TractographyImage') == 8
%%  disp(CAMINO_TXT);
%%else; error('camino.jar doesnt appear to be loaded'); end

%------------------------------------------------------- 
if ndims(Ds)==5;mps=ones(size(Ds,1),size(Ds,2),size(Ds,3),1);end;

if nargin<2 || isempty(seedROI)
  seedROI=ones([size(Ds,1),size(Ds,2),size(Ds,3)],'int16');

end;
%
% Convert the Matlab data to a java array:
DTs = javaArray('misc.DT',size(Ds,1),size(Ds,2),size(Ds,3),1);

for x = 1:size(Ds,1)
  for y = 1:size(Ds,2)
    for z = 1:size(Ds,3)
      DTs(x,y,z,1) = misc.DT(Ds(x,y,z,1,1),...
			     Ds(x,y,z,1,2),...
			     Ds(x,y,z,1,3),...
			     Ds(x,y,z,2,2),...
			     Ds(x,y,z,2,3),...
			     Ds(x,y,z,3,3));
    end
  end
end

%-------------------------------------------------------
% Parse Optional Inputs
%-------------------------------------------------------

% Defaults: (NB not all used?)
%----------
defaults = {
	   'interpolated',      false
	   'voxelCoordinates',  false 
	   'iterations',        5000
	   'minTractPoints'     0
	   'VoxelDims',         [1,1,1]
           'anisThresh',        0
	   'ipThresh',          cos(pi * 80.0 / 180.0)
	   'stepSize',          0
           'waypointVol',       []
	   };

% Parse Options:
%---------------
opts = parseargs(cell2struct(defaults(:,2),defaults(:,1),1),varargin{:});

xVoxelDim   = opts.VoxelDims(1);
yVoxelDim   = opts.VoxelDims(2);
zVoxelDim   = opts.VoxelDims(3);
%"-anisthresh" 
% Initialisations:
%-----------------	   
 
 %-------------------------------------------------------
 % Get Seeds
 %-------------------------------------------------------

 rois = tractography.FreeFormROI.getRegions(seedROI, xVoxelDim, yVoxelDim, zVoxelDim);
 
 %-------------------------------------------------------
 % For the moment: only tensor data, standard tracking
 %-------------------------------------------------------
 
 % 1. Get a 'Tractography Image'
 %------------------------------
 
 
 %% For Danny: here is the call where things could be made much
 %  better: if Voxel..Source could take a Matlab array of n1 x n2 x
 %  n3 x 3 x 3 (which is seen as a java array with corresponding
 %  dimensions in java).
 image = tractography.DT_TractographyImage(DTs,ones([size(seedROI),1]),ones(size(seedROI),'int32'),int32(size(seedROI)),[xVoxelDim, yVoxelDim,zVoxelDim]);
 
 % 2. Fibre Tracker
 %-----------------
 
 if (opts.interpolated) 
   tracker = tractography.LinIntDTEulerFibreTracker(image, opts.stepSize, opts.ipThresh);
 else
   tracker = tractography.NonInterpolatedFibreTracker(image, opts.ipThresh);
 end
 
 
 %% OBSOLETED BY PHIL COOK
 % minWaypointIndex = 0; maxWaypointIndex = 0;
 % if (~isempty(opts.waypointVol))
 %   minMax = getMinMax(opts.waypointVol);
 %   minWaypointIndex = minMax(1); maxWaypointIndex = minMax(2);
 %   numWaypoints = maxWaypointIndex - minWaypointIndex + 1;
 % end
minTractPoints = opts.minTractPoints;
 
 [XDataDim,YDataDim,ZDataDim] = size(seedROI);
 filter = tractography.StreamlineROI_Filter(int32(XDataDim), int32(YDataDim),int32(ZDataDim),xVoxelDim, yVoxelDim, zVoxelDim);
 
 %% NEW? seedVoxels = getSeedVoxels(seedROI, xVoxelDim, yVoxelDim,
 %zVoxelDim);
 
 % The actual fibre tracking
 %--------------------------
 for region = 1:rois.length()
   seeds = rois(region).getVoxelCentres();
   fidx = 1;  

   
   for sp = 1:seeds.length
     %% NEW? x = seedVoxels(sp-1).x;
     %% NEW? y = seedVoxels(sp-1).y;
     %% NEW? z = seedVoxels(sp-1).z;
     %% NEW? numPDs = image.numberOfPDs(x,y,z);
     %% NEW? if numPDs == 0; numPDS = 1; end;
     % PAC - call trackFromSeed(Point3D) instead of trackFromSeed(Point3D, int) 
      seedPoint = seeds(sp);  
       tc = filter.processTracts(tracker.trackFromSeed(seedPoint));
       %% OLD: seedPoint,opts.waypointVol,minWaypointIndex,maxWaypointIndex);
       for fibre = 0:tc.numberOfTracts()-1
	 tract = tc.getTract(fibre);
	 numPoints = tract.numberOfPoints();
	 tr = zeros(numPoints,3);
	 for p = 0:(numPoints-1)
	   point = tract.getPoint(p);
	   tr(p+1,:) = [point.x,point.y,point.z];
	 end
	 %% NB: inefficient! could be sped up
	 if length(tr) > minTractPoints; trs{fidx} = tr; fidx = fidx + 1; end;
       end % fibre
     end % sp
 end
