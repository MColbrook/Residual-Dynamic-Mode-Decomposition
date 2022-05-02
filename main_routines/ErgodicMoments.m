function MU = ErgodicMoments(X,N)
% Compute moments -N to N of vector X using ergodic formula
X=X(:);
M=length(X);
MU=zeros(1,N+1);
MU(1)=X'*X/(M*2*pi);
% MU2=MU;

w=conv(X,conj(flipud(X)));
MU(2:N+1)=transpose(w(length(X)-(1:N)))./(2*pi*((M-1):(-1):M-N));

% for j=1:N
%     MU(j+1)=X(j+1:end)'*X(1:end-j)/(2*pi*(M-j));
% end
% norm(MU-MU2)

MU=[conj(fliplr(MU(2:end))),MU];

end

