function s = SetDefaultParams_CCS(s)
% s = SetDefaultsms(s);
% Sets default smeters
% s: user-specified set of parameters that are used instead of defaults

if (isfield(s, 'p') == 0),
    s.p = 0.2;
end

if (isfield(s, 'delta') == 0),
    s.delta = 0.2;
end
end

