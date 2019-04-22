function [xy ij] = modify( x,y)
global MM
global Lgrid
xy =[ fix(x/Lgrid)*Lgrid+Lgrid/2,fix(y/Lgrid)*Lgrid+Lgrid/2];
ij = [MM+0.5-xy(2)/Lgrid,xy(1)/Lgrid+0.5];
end

