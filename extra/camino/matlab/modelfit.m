function output = modelfit(options, inputdata, mask)
% MODELFIT is a wrapper for the Camino modelfit command.  It fits various
% models to diffusion MRI data and outputs the fitted parameters in each
% voxel of the input data array.  Input and output can either be matlab
% arrays of data files.
%
% MODELFIT(OPTIONS) fits the model specified by OPTIONS.inversion in every
% voxel of OPTIONS.inputfile and outputs the results to options.outputfile
% assuming the acquisition specified in OPTIONS.schemefile.  All those
% fields must exist in the options structure.
%
% MODELFIT(OPTIONS, INPUTDATA) as above, but uses the array INPUTDATA
% instead of reading data from OPTIONS.inputfile.
%
% MODELFIT(OPTIONS, INPUTDATA, MASK) as above, but uses the array MASK as
% a background mask and returns default values in all voxels in which
% mask contains a zero.  MODELFIT(OPTIONS, INPUTDATA) assumes all voxels
% are foreground, unless OPTIONS.bgmask specifies a file containing a
% background mask.  Note that we can only provide a MASK array as a matlab
% argument if INPUTDATA is passed as a matlab argument.
%
% OUTPUT = MODELFIT(OPTIONS, INPUTDATA) as above, but returns the array of
% fitted model parameters in the array OUTPUT rather than writing it to
% OPTIONS.outputfile.
%
% OUTPUT = MODELFIT(OPTIONS, INPUTDATA, MASK) as MODELFIT(OPTIONS,
% INPUTDATA, MASK), but returns the array of fitted models rather than
% outputting it to the file OPTIONS.outputfile.
%
% OUTPUT = MODELFIT(OPTIONS) as above, but reads the input data from
% OPTIONS.inputfile.  Since the program has no way to determine the
% dimensions of the voxel array, the shape of the return array is [NVOX, 1,
% 1, C] where NVOX is the total number of voxels in the image and C is the
% number of components of the output in each voxel.
%
% All options to the MODELFIT program are contained as fields in the
% OPTIONS structure.  The names of the fields reflect the command line
% arguments listed in the modelfit man page.  For example, -inputfile
% <filename> specifies the name of an input file for the modelfit unix
% command.  Equivalently OPTIONS.inputfile specifies the input file name
% for the matlab command.
%
% MODELFIT expects and returns image arrays in voxel order, ie, the first
% index specifies the component in one voxel and the remaining indices
% specify the position of the image voxel.  The following series of 
% commands
%
% setup;
% options.inversion = 1;
% options.schemefile = 'A.scheme';
% options.inputfile = 'SubjectA.Bfloat';
% options.outputfile = 'SubjectA_DT.Bdouble';
% modelfit(options);
%
% Is equivalent to
%
% setup;
% options.inversion = 1;
% options.schemefile = 'A.scheme';
% fid = fopen('SubjectA.Bfloat', 'r', 'b');
% d = fread(fid, 'float');
% fclose(fid);
% idata = reshape(d, [60 128 128 45]); % 60 measurements in each of 128x128x45 voxels.
% odata = modelfit(options, idata);
% fid = fopen('SubjectA_DT.Bdouble', 'w', 'b');
% fwrite(fid, odata, 'double');
% fclose(fid);

% $Id$


% Set up the input to the routine
usedatafile = 0;
if(nargin == 1)
    usedatafile = 1;
    if(~isfield(options, 'inputfile'))
        error('No input array or file specified.');
    end
end
args = opttoargs(options);

usemaskarray = 0;
if(nargin >= 3)
    usemaskarray = 1;
end

% Figure out the output options
outputtofile = 0;
if(nargout == 0)
    outputtofile = 1;
    if(~isfield(options, 'outputfile'))
        error('No return value or output file specified.');
    end
end

% Check essential options are specified
if(~isfield(options, 'schemefile'))
    error('No schemefile specified in options.');
end

% Import the required Camino packages
import tools.*;
import apps.*;
import data.*;

% Initialize the Camino variables
CL_Initializer.CL_init(args);
CL_Initializer.checkParsing(args);
CL_Initializer.initImagingScheme();

inv = ModelFit.getIndexedInversion(CL_Initializer.inversionIndices, CL_Initializer.imPars);

% Initialize the input
if(~usedatafile)
    % Note that the input data array is in voxel order, but we need to
    % permute it to scanner order to construct the MatlabArrayDataSource.
    tools.CL_Initializer.setDataSource(MatlabArrayDataSource(permute(inputdata, [2 3 4 1])));
    
    if(usemaskarray)
        tools.CL_Initializer.setMaskSource(MatlabArrayDataSource(mask));
    end
else
    CL_Initializer.initDataSynthesizer();
end

% Initialize the output options
om = data.OutputManager();
if(~outputtofile)
    if(~usedatafile)
        [comp x y z] = size(inputdata);
    else
        % Use a default array size computed from the size of the data file.
        f = dir(options.inputfile);
        comp = CL_Initializer.imPars.numMeasurements();
        bpv = bytesperval(options);
        x = f.bytes/(bpv*comp);
        y = 1;
        z = 1;
        [x y z]
    end
    om.setOutputArray(x, y, z);
end

% Run the procedure
ModelFit.processVoxels(inv, om);

% Collect the output for return
if(~outputtofile)
    output = om.getOutputArray();
    output = permute(output, [4 1 2 3]);
end

