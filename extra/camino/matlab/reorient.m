function output = reorient(options, inputdata)
% REORIENT is a wrapper for the Camino reorient command.  It reorients
% diffusion tensors for consistency with tissue orientation after an image
% warp.
%
% REORIENT(OPTIONS) reorients using the options specified in the
% OPTIONS stucture.  These are the same as the command line options
% listed in the man page.  The OPTIONS structure must contain arguments
% -inputfile and -outputfile to tell the program where to read the input
% data from and where to output the results.  It must also contain
% -datadims to specify the dimensions of the input file and -voxeldims to
% specify the voxel dimensions.
%
% REORIENT(OPTIONS, INPUTDATA) as above, but uses the array INPUTDATA
% instead of reading data from OPTIONS.inputfile.  The program now
% determines the image dimensions from the input array, so -datadims is not
% necessary.
%
% OUTPUT = REORIENT(OPTIONS, INPUTDATA) as above, but returns the array of
% fitted model parameters in the array OUTPUT rather than writing it to
% OPTIONS.outputfile.
%
% OUTPUT = REORIENT(OPTIONS) as above, but reads the input data from
% OPTIONS.inputfile. -datadims required in OPTIONS.
%

% $Id$


% Set up the input to the routine
usedatafile = 0;
if(nargin == 1)
    usedatafile = 1;
    if(~isfield(options, 'inputfile'))
        error('No input array or file specified.');
    end
    if(~isfield(options, 'datadims'))
        error('Need to specify dimensions of input data file.');
    end
end
args = opttoargs(options);

% Figure out the output options
outputtofile = 0;
if(nargout == 0)
    outputtofile = 1;
    if(~isfield(options, 'outputfile'))
        error('No return value or output file specified.');
    end
end

% Import the required Camino packages
import tools.*;
import apps.*;
import data.*;
import misc.*;
import numerics.*;

% Initialize the Camino options
Reorient.initOptions(args)

% Initialize the input
if(~usedatafile)
    % Note that the input data array is in voxel order, but we need to
    % permute it to scanner order to construct the MatlabArrayDataSource.
    tools.CL_Initializer.setDataSource(MatlabArrayDataSource(permute(inputdata, [2 3 4 1])));
end

% Initialize the output options
om = data.OutputManager();
if(~outputtofile)
    if(~usedatafile)
        [comp x y z] = size(inputdata);
        if(comp<8)
            warning('Expecting DT data with at least 8 components per voxel.');
        end
    else
        x = CL_Initializer.dataDims(1);
        y = CL_Initializer.dataDims(2);
        z = CL_Initializer.dataDims(3);
        [x y z]
    end
    om.setOutputArray(x, y, z);
end

% Run the procedure
Reorient.initVariables();
Reorient.execute(om);
om.close()

% Collect the output for return
if(~outputtofile)
    output = om.getOutputArray();
    output = permute(output, [4 1 2 3]);
end

