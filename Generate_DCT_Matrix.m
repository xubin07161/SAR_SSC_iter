function [DCT]=Generate_DCT_Matrix(n,Pn)

DCT=zeros(n,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:n-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);

end