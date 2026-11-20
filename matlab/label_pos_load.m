function M = label_pos_load(file)
%LABEL_POS_LOAD  Hand-set inline-label positions, keyed by figtag.
%   M = LABEL_POS_LOAD(FILE) reads FILE (default 'label_pos_g.json') into a containers.Map figtag -> spec, where
%   spec = {label, [x y], rot_deg, hide; ...} is what compare_limits_g takes as
%   batch_params.label_pos.  A missing file gives an empty map.
if nargin < 1, file = fullfile(fileparts(mfilename('fullpath')), 'label_pos_g.json'); end
M = containers.Map();
if ~isfile(file), return; end
J = jsondecode(fileread(file));
for i = 1:numel(J)
    L = J(i).labels;
    if isempty(L), continue; end
    if iscell(L), L = [L{:}]; end
    spec = cell(numel(L), 4);
    for k = 1:numel(L)
        spec(k, :) = {L(k).label, [L(k).x, L(k).y], L(k).rot, logical(L(k).hide)};
    end
    M(J(i).figtag) = spec;
end
end
