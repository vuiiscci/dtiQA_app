function bpv = bytesperval(options)
% BYTESPERVAL returns the number of bytes per data item in an input data
% file according to the data type specified in the OPTIONS structure.
% If the OPTIONS structure contains an inputdatatype field, that specifies
% the number of bytes.  If the field is absent, a default setting of
% "double" is assumed.

% $Id$

if(~isfield(options, 'inputdatatype'))
    bpv = 4;
    return;
end

if(strcmp(options.inputdatatype, 'double'))
    bpv = 8;
    return;
elseif(strcmp(options.inputdatatype, 'double'))
    bpv = 8;
    return;
elseif(strcmp(options.inputdatatype, 'float'))
    bpv = 4;
    return;
elseif(strcmp(options.inputdatatype, 'long'))
    bpv = 8;
    return;
elseif(strcmp(options.inputdatatype, 'int'))
    bpv = 4;
    return;
elseif(strcmp(options.inputdatatype, 'short'))
    bpv = 2;
    return;
elseif(strcmp(options.inputdatatype, 'char'))
    bpv = 1;
    return;
elseif(strcmp(options.inputdatatype, 'byte'))
    bpv = 1;
    return;
end

error(['Unrecognized data type: ', options.inputdatatype]);
