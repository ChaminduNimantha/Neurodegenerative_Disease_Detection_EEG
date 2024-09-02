function [ax, X] = BarwithErr(Mean,Error,K)
if nargin<3
    K.k1 = [];
    K.k2 = [];
    K.k3 = [];
    K.k  = [];
end
figure
ax = subplot(2,1,1);
hold on
if min(size(Mean))==1
    Mean = [Mean; 0*Mean];
    Error = [Error; 0*Error];
end
bar(1:size(Mean,1),Mean,'grouped');
%%
% Add errorbars
ngroups = size(Mean,1);
nbars = size(Mean,2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
X{nbars}= [];
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Mean(:,i), Error(:,i), 'k', 'linestyle', 'none','linewidth',1);
    X{i} = x; % Need to draw the significance stars
end
% draw significance asterisks
astr = {}; % significance comparison x locations
p = []; % significance p-value
for ii=1:length(K.k1)
    i = K.k1(ii);
    astr = [astr [X{1}(i) X{3}(i)]];
    p = [p K.p1(ii)];    
end
for ii=1:length(K.k2)
    i=K.k2(ii);
    astr = [astr [X{2}(i) X{3}(i)]];
    p = [p K.p2(ii)];
end
for ii=1:length(K.k3)
    i=K.k3(ii);
    astr = [astr [X{2}(i) X{1}(i)]];
    p = [p K.p3(ii)];
end
sigstar(astr,p);
