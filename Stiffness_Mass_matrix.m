function [ Kg,K,M ] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,mu1,h,density1,he,Cf )
density = density1;
K=spalloc(2*Nnode,2*Nnode,600000);
M=spalloc(2*Nnode,2*Nnode,600000);
% K = zeros(2*Nnode,2*Nnode);
% M = zeros(2*Nnode,2*Nnode);
E = 70e9;
v = 0.33;

% for i1=1:Nelement
% 	x1 = Coordinate(Ielement(i1,1),1); y1 = Coordinate(Ielement(i1,1),2);
% 	x3 = Coordinate(Ielement(i1,2),1); y3 = Coordinate(Ielement(i1,2),2);
% 	x5 = Coordinate(Ielement(i1,3),1); y5 = Coordinate(Ielement(i1,3),2);
% 	x7 = Coordinate(Ielement(i1,4),1); y7 = Coordinate(Ielement(i1,4),2);
% 	[Ke,Me] = Ksolve(E,v,h,x1,y1,x3,y3,x5,y5,x7,y7);
% 
% 	Kr = Ke;   %%%%%%
% 	Mr = density*Me*he^2;
% 	a=Ielement(i1,:);
% 	for j=1:8
% 		for k=1:8
% 			K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Kr(j*2-1:j*2,k*2-1:k*2);
% 			M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Mr(j*2-1:j*2,k*2-1:k*2);
% 		end
% 	end
% end
% 
% Kg=(K-Cf^2*M)/mu1;
% end

	x1 = Coordinate(Ielement(1:Nelement,1),1); y1 = Coordinate(Ielement(1:Nelement,1),2);
	x3 = Coordinate(Ielement(1:Nelement,2),1); y3 = Coordinate(Ielement(1:Nelement,2),2);
	x5 = Coordinate(Ielement(1:Nelement,3),1); y5 = Coordinate(Ielement(1:Nelement,3),2);
	x7 = Coordinate(Ielement(1:Nelement,4),1); y7 = Coordinate(Ielement(1:Nelement,4),2);
	
	
	[Ke,Me] = Ksolve(E,v,h,x1,y1,x3,y3,x5,y5,x7,y7);

	Kr = Ke;   %%%%%%
	Mr = density*Me*he^2;
	a=Ielement(i1,:);
	for j=1:8
		for k=1:8
			K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=K((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Kr(j*2-1:j*2,k*2-1:k*2);
			M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=M((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2) + Mr(j*2-1:j*2,k*2-1:k*2);
		end
	end

Kg=(K-Cf^2*M)/mu1;
end
