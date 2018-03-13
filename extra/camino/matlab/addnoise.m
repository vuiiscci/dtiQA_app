function output = addnoise(options, inputdata)

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
theaddnoise = apps.AddNoise(args);

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
%         if(comp<8)
%             warning('Expecting DT data with at least 8 components per voxel.');
%         end
    else
    error('cannot read from a file and output to a matlab array');
%         x = CL_Initializer.dataDims(1);
%         y = CL_Initializer.dataDims(2);
%         z = CL_Initializer.dataDims(3);
%         [x y z]
    end
    om.setOutputArray(x, y, z);
end

% Run the procedure
theaddnoise.initDefaultVals();
theaddnoise.initVariables();
theaddnoise.execute(om);
% theaddnoise.execute(om);
% om.close()

% Collect the output for return
if(~outputtofile)
    output = om.getOutputArray();
    output = permute(output, [4 1 2 3]);
end

om.close()

%$ bin/datasynth -seed 0  -testfunc 1 -schemefile test/bmx6.scheme -voxels 10 -snr 16 > test/matlab_in.in
%$ cat matlab_out.out | ../bin/float2txt
%$ cat matlab_out.out | ../bin/float2txt > matlab_out.text

%fid1 = fopen('matlab_in.in'); [A,n] = fread(fid1, 660, 'float', 'b'); fclose(fid1);


