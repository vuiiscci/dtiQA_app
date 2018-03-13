function output = track_basic(options, inputdata)
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

% Initialize the output options
om = data.OutputManager();

% Run the procedure
theST.initDefaultVals();
theST.initOptions(args);
theST.initVariables();
theST.execute(om);

om.close();


