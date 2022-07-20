function pattern = generate_pattern(m)

% pattern = generate_pattern(m)
%
% given m=[m1 m2 m3 m4 ... mN] generares a matrix that can be used to generate the cartesian product
% of {1,2,...,m1} x {1,2,...,m2} x {1,2,...m3} x ....
% 
% e.g.
%
% generate_pattern([3 2 2])
% 
% pattern =
% 
%       1      2      3      1      2      3      1      2      3      1      2      3
%       1      1      1      2      2      2      1      1      1      2      2      2
%       1      1      1      1      1      1      2      2      2      2      2      2


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N=length(m);

% it is convenient from a computational point of view to generate the pattern as unsigned int to create the pattern
if max(m)<=intmax('uint8')
    pattern=zeros(N,prod(m),'uint8');
elseif max(m) <=intmax('uint16')
    warning('SparseGKit:uint16','more than 255 points are asked in one direction, using uint16 to handle this')
    pattern=zeros(N,prod(m),'uint16');
else
    warning('SparseGKit:double','more than 65535 points are asked in one direction, using double precision to handle this')
    pattern=zeros(N,prod(m));
end


% the algorithm works one direction at a time. at the n-th iteration the n-th row of pattern is generated,
% by repeating q times the vector BASE, which containes itselt a sequence,
% obtained repeating p times each number from j=1 to j=m(n), e.g. 
% generate_pattern([3 2 2])
% 
% pattern =
% 
%       1      2      3      1      2      3      1      2      3      1      2      3
%       1      1      1      2      2      2      1      1      1      2      2      2
%       1      1      1      1      1      1      2      2      2      2      2      2

for k=1:N
    p = prod([1 m(1:k-1)]);
    q = prod([m(k+1:end) 1]);
    
    % length of base vector
    lb=p*m(k);
    base = zeros(1,lb);
    
    % generate base vector
    bb=1;
    for j=1:m(k)
        base(bb:bb+p-1)=j;
        bb=bb+p;
    end
    
    % repeat base vector
    pp=1;
    for j=1:q
        pattern(k,pp:pp+lb-1)=base;
        pp=pp+lb;
    end
    
end