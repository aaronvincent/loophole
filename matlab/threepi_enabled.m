function tf = threepi_enabled()
% threepi_enabled  Toggle for the 3pi channels of the SP/PP sector.
%
% Off by default.  Start MATLAB with LOOP_3PI=1 in its environment to turn on
% chi chi -> 3pi and phi -> 3pi together.  An environment variable, not a
% workspace flag, because the channels sit under parfor deep in the solvers:
% parpool workers inherit the environment of the MATLAB that launches them.
persistent cached
if isempty(cached)
    cached = strcmp(getenv('LOOP_3PI'), '1');
end
tf = cached;
end
