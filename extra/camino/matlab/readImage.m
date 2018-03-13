function image = readImage(file)
%
% READIMAGE reads an Analyze, NIFTI-1 or ITK MetaIO image file. 
% Specify the full path to the file including the extension.
%
% 
% image = READIMAGE(file) returns all volumes in the image. For  
% a 4D image with C components and spatial dimensions X, Y, Z, the 
% returned image dimensions [X][Y][Z][C].
%
% 

ih = javaMethod('readHeader', 'imaging.ImageHeader', file);

image = ih.readVolumeData();
