function Kg = Stiffness(Cf)
K1 = ImportKM("K.txt");
M1 = ImportKM("M.txt");
mu1 = 26.1*10^9;
Kin = sparse(K1.x,K1.y,K1.v);
Min = sparse(M1.x,M1.y,M1.v);
Kg=(Kin-Cf^2*Min)/mu1;
% full(Kg)
% Kg1 = Kg
end