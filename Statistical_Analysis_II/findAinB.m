function [I tt] = findAinB(A,B)
% finds the indices in B that corresponds to values in A
offsetvec =A; %HDR(4).offsetvec;
T = B;%HDR(1).modelT;
I = [];
tt = [];
for t = 1:length(offsetvec)
    i = find(T == offsetvec(t));
%     if ~isempty(i)
        tt = [tt offsetvec(t)];
        I = [I i];
%     end
    
end
