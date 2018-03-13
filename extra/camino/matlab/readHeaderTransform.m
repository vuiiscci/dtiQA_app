function [trans] = readHeaderTransform(hdrFile, inv)
%
% [trans] = readHeaderTransform(hdrFile, inv)
%
% Gets the voxel to physical space transformation matrix (or its inverse) for the header.
% 
% The returned matrix is a 4x4 affine transformation. This has the form
%
%   voxToPhys = [A11, A12, A13, t_x; A21, A22, A23, t_y ; A31, A32, A33, t_z; 0, 0, 0, 1]
%
% where A_ij are elements of a 3x3 affine transform, and T = [t_x; t_y; t_z] is a translation. 
%
% For the inverse transform, physToVox is composed of inv(A) and a translation -inv(A) * T
%
% Given a voxel index [i;j;k;1], the physical space coordinate phys = [x;y;z;1] is
%
%   phys = voxToPhys * [i; j; k; 1] 
%
% Args:
%     hdrFile = file to read
%     inv = if 1, return the inverse transform physToVox, if 0 or not specified, return the forward transform voxToPhys
% 

if (nargin < 2) 
  inv = 0;
end

ih = javaMethod('readHeader', 'imaging.ImageHeader', hdrFile);

if (inv == 1) 
  trans = ih.getPhysicalToVoxelTransform().entries;
else
  trans = ih.getVoxelToPhysicalTransform().entries;
end

