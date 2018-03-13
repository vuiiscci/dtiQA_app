function args = opttoargs(options)
% OPTTOARGS converts an options structure to a list of arguments for Camino
% commands.  Each field becomes a flag by prepending a '-' and each value
% is placed in the args array directly after the flag.

% $Id$

argind = 1;
f = fieldnames(options);
for j=1:length(f)
    fldname = f{j};
    args(argind) = java.lang.String(['-' fldname]);
    argind = argind + 1;
    if(strcmp(class(getfield(options, fldname)), 'char'))
        % Split on white space.
        remain = getfield(options, fldname);
        while true
            [str, remain] = strtok(remain);
            if isempty(str),  break;  end
            args(argind) = java.lang.String(str);
            argind = argind + 1;
        end
    else
        for i=getfield(options, fldname)
            args(argind) = java.lang.String(sprintf('%d', i));
            argind = argind + 1;
        end
    end
end


    