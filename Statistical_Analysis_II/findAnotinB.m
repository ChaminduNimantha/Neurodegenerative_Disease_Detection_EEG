function [I C] = findAnotinB(A,B)
% finds the indices I and elements C in A that corresponds to values not in B
offsetvec =B; %HDR(4).offsetvec;
T = A;%HDR(1).modelT;
I = [];
tt = [];
for t = 1:length(offsetvec)
    i = find(T == offsetvec(t));
%     if ~isempty(i)
        I = [I i];
%     end
    
end

tt = A(setdiff(1:length(A),I));
I = [];
for t = 1:length(tt)
    i = find(T == tt(t));
        I = [I i];   
end
C = tt;