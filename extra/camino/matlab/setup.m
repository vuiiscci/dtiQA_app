% check if  exists and variable is set, otherwise try to find

CAMINO_PATH = '.';

CAMINO = [CAMINO_PATH,'/camino.jar'];
%%global CAMINO_message;
lin = '===============================================================';
C0 = ' __                          __';
C1 = '/  \   /\  |\    /| I |\  | /  \';
C2 = '|     /__\ | \  / | | | \ | |  |';
C3 = '\__/ /    \|  \/  | I |  \| \__/';

C = sprintf('\t%s\n\t%s\n\t%s\n\t%s',C0,C1,C2,C3);

CAMINO_TXT = sprintf('%s\n\t << This software uses >>\n\n%s\n\n\nSee http://www.cs.ucl.ac.uk/research/medic/camino/\n%s',lin,C,lin);
%-------------------------------------------------------

if exist(CAMINO) ~= 2
  w = [];
  if isunix; [s,w] = system('locate camino.jar'); end;
  % XP: get-childitem \\servername\c$ -include camino.jar -recurse -name 
  % separate different lines as a cell array? and try all lines
  disp(w);
  CAMINO = input('Input camino.jar path manually? (with quotes!) ');
  if    exist(CAMINO) ~= 2;  warning('no valid camino.jar');
  else; javaaddpath(CAMINO); end
else 
  javaaddpath(CAMINO); 
end;

if exist(CAMINO) == 2; disp(CAMINO_TXT); end;
%%if exist(CAMINO) == 2; warndlg('This software uses CAMINO;','CAMINO');end
%-------------------------------------------------------

import apps.*;
import imaging.*;
import misc.*;
import numerics.*;
import tools.*;
import data.*;
import tractography.*;
import Jama.*;
import Jama.util.*;

import java.util.Random;
import java.util.zip.*;
import java.io.*;
import java.lang.*;

%-------------------------------------------------------
warning off;
global CAMINO_PATH CAMINO_TXT CAMINO_message;
warning on;