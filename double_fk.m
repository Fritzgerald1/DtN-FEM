function [ D ] = double_fk(kh,kth1,kth2,mu1,mu2 )
alpha1 = sqrt((kth1/kh)^2-1);
alpha2 = -sqrt((kth1/kh)^2-1);
alpha3 = sqrt((kth2/kh)^2-1);
alpha4 = -sqrt((kth2/kh)^2-1);
matrix(1,1) = exp(-1i*kh*alpha1);
matrix(1,2) = -exp(-1i*kh*alpha2);
matrix(2,3) = exp(1i*kh*alpha3);
matrix(2,4) = -exp(1i*kh*alpha4);
matrix(3,1) = mu1*alpha1/(mu2*alpha3);
matrix(3,2) = mu1*alpha2/(mu2*alpha3);
matrix(3,3) = -1;
matrix(3,4) = 1;
matrix(4,1:2) = 1;
matrix(4,3:4) = -1;
[L,U] = lu(matrix);
D = log(abs(U(1,1))) + log(abs(U(2,2)))+log(abs(U(3,3)))+log(abs(U(4,4)));
end

