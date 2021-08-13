clc
clear all
% clf
%%
% figure(9)
datfiles  = dir('*nTPetu*');
for k = 1 : length(datfiles)
    data = load(datfiles(k).name);
    GetMaximumP(k, :) = data(:, 4);
end

[max(GetMaximumP(:, :))]'
