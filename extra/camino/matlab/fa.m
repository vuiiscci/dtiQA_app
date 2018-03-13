function output = fa(options, inputdata)

% Set up the input to the routine
usedatafile = 0;
if(nargin == 1)
    usedatafile = 1;
    if(~isfield(options, 'inputfile'))
        error('No input array or file specified.');
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
import data.*;
import imaging.DW_Scheme;
import tools.CL_Initializer;

import java.io.IOException;
import java.util.Random;
import java.util.logging.Logger;

import numerics.MTRandom;

% Initialize the Camino options
thefa = apps.FracAnis(args);

% Initialize the input
if(~usedatafile)
    % Note that the input data array is in voxel order, but we need to
    % permute it to scanner order to construct the MatlabArrayDataSource.
    tools.CL_Initializer.setDataSource(MatlabArrayDataSource(permute(inputdata, [2 3 4 1])));
end

% Initialize the output options
% assume the data is in the form 8*x*y*z
om = data.OutputManager();
if(~outputtofile)
    if(~usedatafile)
        [comp x y z] = size(inputdata); %so comp=8
    else
        error('cannot read from a file and output to a matlab array');
    end
    om.setOutputArray(x, y, z);
end

% Run the procedure
thefa.initDefaultVals();
thefa.initVariables();
thefa.execute(om);


% Collect the output for return
if(~outputtofile)
    output = om.getOutputArray();
    output = permute(output, [4 1 2 3]);
end

om.close();


