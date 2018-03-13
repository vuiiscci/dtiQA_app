function output = twotenfit(options, inputdata)


options.inversion = 31;
if isfield(options.cyl, 'cylsym')
    options.inversion = 11;
end

output = modelfit(options, inputdata);