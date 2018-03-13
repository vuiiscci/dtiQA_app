Introductory user guide for the Camino toolkit
______________________________________________

This document provides a quick introduction to the Camino toolkit to
help new users get started.  Each program in the toolkit has a man
page, which provides more detailed information on specific tools.

There are also tutorials and general documentation on the website:

http://camino.org.uk


Installation
============

Download the latest version from Sourceforge - you will get a file that looks like

camino-code-163f67cbf550560aa351b3d0a3bbbd7a22863cb4.zip

The string following code- is the commit hash. Make a note of this as it uniquely
identifies the exact version of the code you are running.

We recommend moving the unzipped code to a directory named camino:

 $ unzip camino-code-commithash.zip

 $ mv camino-code-commithash camino

Next, compile the code. You will need to have the Oracle Java SDK 1.7 or later and have 
javac and java in your system path.

 $ cd camino
 $ make


Cygwin
======

Although the individual components of Camino will run from a commandline under windows, 
in order to get the most out of Camino it is necessary to have a UNIX-like shell 
environment that allows data pipes and redirection. Without these facilities, it is 
extremely difficult to use Camino in the way it was designed.

Fortunately, it is extremely easy to install Camino under windows with Cygwin, this 
section explains the procedure step-by-step. Firstly, download the Cygwin installer 
from Sun is installed. Make sure you have the SDK as well as the usual JRE! 
This is not installed under windows as standard!

Once Cygwin and the Java SDK are installed on your system, check that the location of 
the Java SDK is added to your windows path. You can do this as follows:

1) From the desktop, click start and right-click on "My Computer"
2) In the "System Properties" window that appears, select the "advanced" tab
3) Click the "Environment variables" button.
4) Highlight the "path" variable and click the "edit" button
5) If the path to your Java SDK is not in the list, add the FULL path to the end 
   of the list, using a semi-colon to separate it from the previous entry

Now start Cygwin and follow the instructions for installing Camino under linux/unix. 
For instructions on how to install geomview under windows and Cygwin, click here (SaVi 
is not required)


Getting Started
===============

Read the camino man page for an overview of the toolkit:

(% cd camino)
% man -M man camino

Further help getting started is available from the Camino website

http://camino.org.uk



Testing (for developers)
=======================

There are two sets of tests distributed with the code. The unit tests are low 
level, object-orientated tests of the Java classes. To make and run the unit 
tests: 

% cd camino/test
% make
% ./runtest.sh

The other level of testing in Camino is at the application level. To run this test:

 $ cd camino
 $ test/ScriptTest > testoutput.mymachine

This is a regression test; you should see consistent results with different versions of 
Camino unless the developers have deliberately changed the results or added new tests. 
However, the results will not be identical on different platforms. 

If you plan to modify Camino code, then you should run ScriptTest on your local machine 
before making changes and save the results. You can then make your changes and compare the 
output of ScriptTest to your baseline results.


Contributions
=============

The developers welcome contributions to the Camino toolkit. Please contact the developers at 
camino@cs.ucl.ac.uk if you would like to contribute code.



Help, bug reports, feature requests
=======================

Please see the Camino website for documentation and tutorials. If you still need help, there 
is a Camino users mailing list that can be joined via the website.

Bugs and feature requests can be added on Sourceforge under Tickets. Bug reports should include
sufficient information to reproduce the problem.
