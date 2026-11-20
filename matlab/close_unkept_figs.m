function close_unkept_figs()
%CLOSE_UNKEPT_FIGS  close all force, except figures kept open on purpose
%   (those carry appdata 'label_edit', set by compare_limits_g's save step).
f = findall(0, 'Type', 'figure');
keep = arrayfun(@(h) isappdata(h, 'label_edit'), f);
close(f(~keep), 'force');
end
