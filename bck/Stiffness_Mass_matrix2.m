function [Kg,K,M] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,density,mu,Cf)
E = 70e9;
v = 0.33;
K = zeros(2*Nnode,2*Nnode);
M = zeros(2*Nnode,2*Nnode);
K1 = zeros(16,16,Nelement);
M1 = zeros(16,16,Nelement);

for i1=1:Nelement
	x1 = Coordinate(Ielement(i1,1),1); y1 = Coordinate(Ielement(i1,1),2);
	x3 = Coordinate(Ielement(i1,2),1); y3 = Coordinate(Ielement(i1,2),2);
	x5 = Coordinate(Ielement(i1,3),1); y5 = Coordinate(Ielement(i1,3),2);
	x7 = Coordinate(Ielement(i1,4),1); y7 = Coordinate(Ielement(i1,4),2);
	[Ke,Me] = Ksolve(E,v,x1,y1,x3,y3,x5,y5,x7,y7);
	Kr = Ke;
	Mr = density*Me;
	for j=1:16
		for k=1:16
			K1(j,k,i1) = Kr(j,k);
			M1(j,k,i1) = Mr(j,k);
		end
	end
end

for i1 = 1:Nelement
	a=Ielement(i1,:);
	for j = 1:8
		for k = 1:8
			K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Kr(j*2-1:j*2,k*2-1:k*2);
			M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Mr(j*2-1:j*2,k*2-1:k*2);
		end
	end
end

Kg = (K-Cf^2*M)/mu;
end
