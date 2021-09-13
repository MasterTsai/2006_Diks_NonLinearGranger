function [ur,uv] = var_resid(rp,vp);

% VAR regressions for the calendar-adjusted data
% grid search, max lag 20 for independent variable, 40 for the dependent
% first treat vp as independent variable

% generate delay vectors
%
% old method (depends on GARCH toolbox):
% rpLAG = lagmatrix(rp,[1:40]);
% vpLAG = lagmatrix(vp,[1:20]);

rpLAG = delayvectors(rp,40);
vpLAG = delayvectors(vp,20);

len = length(rp);

% initialize the best AIC found
hbest = 100000;

for k=1:20
    for m=1:40
        regressors = [ones(len-40,1) vpLAG(41:end,1:k) rpLAG(41:end,1:m)];
        % OLS regression
        %
        % Old implementation (depends on GARCH toolbox)
        % [b,bint,r,rint,stats] = regress(rp(41:end),regressors);
        %
        % New implementation
        b = regressors\rp(41:end);
        % Calculate residuals
        r = rp(41:end) - regressors*b;
        s = std(r);
        
        % Determine AIC
        logl = 0;
        for i=41:len,
            logl = logl - log(s) - 0.5*(r(i-40)*r(i-40)/(s*s));
        end
        AIC = -2*logl + 2*(k+m+2);
        [AIC]
        % Update best result so far
        if (AIC/len < hbest),
            hbest = AIC/len;
            kbest = k;
            mbest =m;
        end
        [k m kbest mbest]
    end
end

[kbest mbest hbest]

% Repeat regression for best AIC found and store the residuals
regressors = [ones(len-40,1) vpLAG(41:end,1:kbest) rpLAG(41:end,1:mbest)];
%[b,bint,ur,rint,stats] = regress(rp(41:end),regressors);
b = regressors\rp(41:end);
% Calculate residuals
ur = rp(41:end) - regressors*b;


% now the same again, but with rp as the independent variable

% Old implementation
%
% rpLAG = lagmatrix(rp,[1:20]);
% vpLAG = lagmatrix(vp,[1:40]);
rpLAG = delayvectors(rp,20);
vpLAG = delayvectors(vp,40);

hbest = 100000;

for k=1:20
    for m=1:40
        [k m]
        regressors = [ones(len-40,1) rpLAG(41:end,1:k) vpLAG(41:end,1:m)];
        %[b,bint,v,rint,stats] = regress(vp(41:end),regressors);
        b = regressors\vp(41:end);
        v = vp(41:end) - regressors*b;
        s = std(v);
        logl = 0;
        for i=41:len,
            logl = logl - log(s) - 0.5*(v(i-40)*v(i-40)/(s*s));
        end
        AIC = -2*logl + 2*(k+m+2);
        [AIC]
        if (AIC/len < hbest),
            hbest = AIC/len;
            kbest = k;
            mbest =m;
        end
        [k m kbest mbest]
    end
end

[kbest mbest hbest]

regressors = [ones(len-40,1) rpLAG(41:end,1:kbest) vpLAG(41:end,1:mbest)];
%[b,bint,uv,rint,stats] = regress(vp(41:end),regressors);
b = regressors\vp(41:end);
uv = vp(41:end) - regressors*b;

return
