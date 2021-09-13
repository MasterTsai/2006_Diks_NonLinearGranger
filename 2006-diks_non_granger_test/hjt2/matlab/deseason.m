function [rp,vp] = deseason(data,textdata);

days = flipud(textdata(:,1));
days = days(1:end-1);

volume = flipud(data(:,5));
price = flipud(data(:,6));

% convert dates to day numbers
daynums = datenum(days,1900);

s = datevec(daynums);

month = s(:,2,:,:,:,:);

subplot(3,2,1)
plot(daynums,log(price));
subplot(3,2,2)
plot(daynums,log(volume));

Ndata = length(price);

dayofweek = zeros(Ndata-1,5);
monthofyear = zeros(Ndata-1,11);
logret = zeros(Ndata-1,1);
logvolinc = zeros(Ndata-1,1);
jul = zeros(Ndata-1,1);

for i=2:Ndata,
    logret(i-1) = log(price(i)) - log(price(i-1));
    logvolinc(i-1) = 100*(log(volume(i)) - log(volume(i-1)));
    jul(i-1) = datenum(daynums(i));
    dayofweek(i-1,weekday(daynums(i))-1) = 1;
    if month(i) < 12,
       monthofyear(i-1,month(i)) = 1;
   end
end

dummies = [dayofweek, monthofyear];

subplot(3,2,3)
plot(jul,logret);
subplot(3,2,4)
plot(jul,logvolinc);

%
% OLS regressions (mean equation)
%
% Old implementation (depends on GARCH toolbox) 
% [b,bint,rets,rint,stats] = regress(logret, dummies);
%
% New implementation
b = dummies\logret;

% Calculate OLS residuals
rets = logret-dummies*b;

% OLS regression (variance equation)
%
%[b,bints,err,rint,stats] = regress(log(rets.*rets), dummies);
b = dummies\log(rets.*rets);

fac = exp(dummies*b/2);
rp = rets./fac;

%
% OLS regressions
%
%[b,bint,vols,rint,stats] = regress(logvolinc, dummies);
b = dummies\logvolinc;
vols = logvolinc-dummies*b;

%[b,bints,err,rint,stats] = regress(log(vols.*vols), dummies);
b = dummies\log(vols.*vols);

fac = exp(dummies*b/2);
vp = vols./fac;

subplot(3,2,5)
plot(jul,rp);

subplot(3,2,6)
plot(jul,vp);

return
