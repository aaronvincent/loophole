function tf = pipiphi_enabled()
% pipiphi_enabled  Toggle for the chi chi -> pi pi phi channels of SP/PP and
% chi chi -> pi+pi-A' of axialVector (pipiA_weight.m).
%
% Off by default.  LOOP_PIPIPHI=1 in the environment of the MATLAB process turns
% them on, including in parpool workers; see threepi_enabled.m.
persistent cached
if isempty(cached)
    cached = strcmp(getenv('LOOP_PIPIPHI'), '1');
end
tf = cached;
end
