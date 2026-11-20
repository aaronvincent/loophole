function name = full_sigma_name(interaction)
% Map a sigma_chiN interaction name to its FULL (angular-integrated) variant.
%
% DMnucxsec.m selects the sigma_chiN approximation from the interaction name: the
% heavy/light PREFIX picks the heavy-mediator limit, the light-mediator closed forms,
% or the full angular-integrated form (see DMnucxsec.m). The canonical observable
% pipeline (model_observables.m, compare_*.m) wants the accurate FULL sigma even when
% the relic model is a heavy- or light-mediator approximation (and to stay consistent
% with the Julia/scan side, which is full-only). It routes the interaction string
% through here so 'heavyScalar*'/'lightScalar*' -> their 'scalar*' full sibling.
% scalarPP and axialVector already return the full form, so they pass through unchanged.
    switch interaction
        case {'heavyScalar',   'lightScalar'},   name = 'scalar';
        case {'heavyScalarPS', 'lightScalarPS'}, name = 'scalarPS';
        case {'heavyScalarSP', 'lightScalarSP'}, name = 'scalarSP';
        case {'heavyScalarPP', 'lightScalarPP'}, name = 'scalarPP';
        otherwise,                               name = interaction;
    end
end
