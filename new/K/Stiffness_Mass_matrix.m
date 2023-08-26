function [K,M] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,lambda,mu,density)
E = mu*(2*mu+3*lambda)/(mu+lambda);
v = lambda/(2*(mu+lambda));
E = 70e9;
v = 0.33;

K = zeros(2*Nnode,2*Nnode);
M = zeros(2*Nnode,2*Nnode);
K1 = zeros(16,16,Nelement);
M1 = zeros(16,16,Nelement);

parfor i1=1:Nelement
	x1 = Coordinate(Ielement(i1,1),1); y1 = Coordinate(Ielement(i1,1),2);
	x2 = Coordinate(Ielement(i1,2),1); y2 = Coordinate(Ielement(i1,2),2);
	x3 = Coordinate(Ielement(i1,3),1); y3 = Coordinate(Ielement(i1,3),2);
	x4 = Coordinate(Ielement(i1,4),1); y4 = Coordinate(Ielement(i1,4),2);
	[Ke,Me] = Ksolve(E,v,x1,y1,x2,y2,x3,y3,x4,y4);
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
			K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + K1(j*2-1:j*2,k*2-1:k*2,i1);
			M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + M1(j*2-1:j*2,k*2-1:k*2,i1);
		end
	end
end
end
