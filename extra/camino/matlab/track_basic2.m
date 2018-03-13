function output = track_basic2(options, inputdata)
% TRACK is a wrapper for the Camino track command.  

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

import imaging.*;
import inverters.ModelIndex;
import misc.LoggedException;

import numerics.*;
import tools.*;
import tractography.*;

import java.util.Random;
import java.util.zip.*;
import java.util.logging.*;

import java.io.*;


theST = apps.StreamlineTractography(args);

% Initialize the input
if(~usedatafile) %ie. if the input is already in matlab
    % Note that the input data array is in voxel order, but we need to
    % permute it to scanner order to construct the MatlabArrayDataSource.
    tools.CL_Initializer.setDataSource(MatlabArrayDataSource(inputdata));
end

% Initialize the output options
om = data.OutputManager();

if(~outputtofile)
    if(~usedatafile)
%        [comp x y z] = size(inputdata)
        size(inputdata)
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
theST.initDefaultVals();
theST.initOptions(args);
theST.initVariables();
theST.execute(om);

% Collect the output for return
if(~outputtofile)
    output = om.getOutputArray();
%    output = permute(output, [4 1 2 3]);
end

om.close();


