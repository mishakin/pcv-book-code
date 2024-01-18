function args = parseArguments(args, varargin)
% PARSEARGUMENTS parses a list of pairs {ArgName, ArgValue} 
% into structure of default arguments
%
% 02/13 Falko Schindler

%% convert argument structure to argument pairs
if ~isempty(varargin) && isstruct(varargin{1})
    arg1 = [fieldnames(varargin{1}), struct2cell(varargin{1})]';
    varargin = {arg1{:}, varargin{2 : end}};
end

%% process pairs of arguments
fields = fieldnames(args);
for i = 1 : 2 : numel(varargin)
    target = textscan(varargin{i}, '%s', 'Delimiter', '.');
    try
        getfield(args, target{1}{:});
    catch
        if numel(target{1}) == 1 && any(strcmpi(varargin{i}, fields))
            target = {fields(find(strcmpi(varargin{i}, fields), 1))};
        else
            error('Unknown parameter ''%s''', varargin{i});
        end
    end
    args = setfield(args, target{1}{:}, varargin{i + 1});
end
