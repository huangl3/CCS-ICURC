function s = SetDefaultParams_ICURC(s)
% s = SetDefaultsms(s);
% Sets default smeters
% s: user-specified set of parameters that are used instead of defaults

if (isfield(s, 'TOL') == 0)
    s.TOL = 1e-4;
end

if (isfield(s, 'max_ite') == 0)
    s.max_ite = 500;
end

if (isfield(s, 'eta') == 0)
    s.eta = [1, 1, 1];
    s.steps_are1 = true;
elseif (s.eta == [1, 1, 1])
    s.steps_are1 = true;
else
    s.steps_are1 = false;
end

end