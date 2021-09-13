function [y] = delayvectors(x,m)

s = length(x);

y = zeros(s,m);

for i=1:m
    for j=1:s
        t = j-i;
        if (t >= 1)
            y(j,i) = x(t);
        else
            y(j,i) = NaN;
        end
    end
end

return;