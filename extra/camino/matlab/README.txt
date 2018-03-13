Using Camino in matlab
----------------------

Add the jar file:

  javaaddpath('/path/to/camino.jar')

And the commands themselves

  addpath('/path/to/camino/matlab')


You will need to increase the maximum memory space that matlab allows
java to use.  Create a text file called java.opts containing the
string "-Xmx512M" in directory $MATLABROOT/bin/$ARCH directory (or in
the MATLAB startup directory). $MATLABROOT is the MATLAB root
directory and $ARCH is your system architecture.

You can find the MATLAB root directory on your machine, by typing

% matlabroot

at the MATLAB Command Prompt:


You can find the system architecture by executing the following at the
MATLAB Command Prompt:

% computer('arch')
