% function [ Kg,K,M ] = Stiffness_Mass_matrix( Nnode,Nelement,Ielement,Coordinate,mu1,mu2,h,density1,density2,he,Cf )
ngq=8;
[xgq,wgq]=gauss_integration(ngq);
%K=sparse(Nnode,Nnode);
%M=sparse(Nnode,Nnode);
K2=spalloc(Nnode,Nnode,600000);
M2=spalloc(Nnode,Nnode,600000);
for i1=1:Nelement
    Ke2=zeros(8,8);
    Me2=zeros(8,8);
    for i2=1:ngq
        for i3=1:ngq
            J =zeros(2,2);
            s=xgq(i2);
            t=xgq(i3);
            Ndt(1)=0.25*(1-s)*(s+2*t);
            Ndt(2)=0.25*(1+s)*(2*t-s);
            Ndt(3)=0.25*(1+s)*(s+2*t);
            Ndt(4)=0.25*(1-s)*(2*t-s);
            Ndt(5)=0.5*(s^2-1);
            Ndt(6)=-t*(1+s);
            Ndt(7)=0.5*(1-s^2);
            Ndt(8)=t*(s-1);
            %J(2,1)=Ndt(1)*Coordinate(Ielement(i1,1),1)+Ndt(2)*Coordinate(Ielement(i1,2),1)+Ndt(3)*Coordinate(Ielement(i1,3),1)+Ndt(4)*Coordinate(Ielement(i1,4),1);
            %J(2,2)=Ndt(1)*Coordinate(Ielement(i1,1),2)+Ndt(2)*Coordinate(Ielement(i1,2),2)+Ndt(3)*Coordinate(Ielement(i1,3),2)+Ndt(4)*Coordinate(Ielement(i1,4),2);    
            Nds(1)=0.25*(1-t)*(2*s+t);
            Nds(2)=0.25*(1-t)*(2*s-t);
            Nds(3)=0.25*(1+t)*(2*s+t);
            Nds(4)=0.25*(1+t)*(2*s-t);
            Nds(5)=s*(t-1);
            Nds(6)=0.5*(1-t^2);
            Nds(7)=-s*(1+t);
            Nds(8)=0.5*(t^2-1);
            %J(1,1)=Nds(1)*Coordinate(Ielement(i1,1),1)+Nds(2)*Coordinate(Ielement(i1,2),1)+Nds(3)*Coordinate(Ielement(i1,3),1)+Nds(4)*Coordinate(Ielement(i1,4),1);
            %J(1,2)=Nds(1)*Coordinate(Ielement(i1,1),2)+Nds(2)*Coordinate(Ielement(i1,2),2)+Nds(3)*Coordinate(Ielement(i1,3),2)+Nds(4)*Coordinate(Ielement(i1,4),2);
            for i33 = 1:ngq
                J(2,1) = J(2,1)+Ndt(i33)*Coordinate(Ielement(i1,i33),1);
                J(2,2) = J(2,2)+Ndt(i33)*Coordinate(Ielement(i1,i33),2);
                J(1,1) = J(1,1)+Nds(i33)*Coordinate(Ielement(i1,i33),1);
                J(1,2) = J(1,2)+Nds(i33)*Coordinate(Ielement(i1,i33),2);
            end
            % NT*N %
            N1=-0.25*(1-s)*(1-t)*(1+s+t);
            N2=0.25*(1+s)*(1-t)*(s-t-1);
            N3=0.25*(1+s)*(1+t)*(s+t-1);
            N4=0.25*(1-s)*(1+t)*(t-s-1);
            N5=0.5*(1-s^2)*(1-t);
            N6=0.5*(1+s)*(1-t^2);
            N7=0.5*(1-s^2)*(1+t);
            N8=0.5*(1-s)*(1-t^2);
            N=[N1 N2 N3 N4 N5 N6 N7 N8];
            NTN=N.'*N;
            % BTCB %
            B=J\[Nds;Ndt];
                mu = mu1;
            C=[mu 0;0 mu];
            BTCB=B.'*C*B;
            Ke2=Ke2+h*wgq(i2)*wgq(i3)*BTCB*det(J);
            Me2=Me2+h*wgq(i2)*wgq(i3)*NTN*det(J);
        end
    end
    % combine the local one into the global one %
    Kr2=Ke2;   %%%%%%
    for i4=1:size(Kr2,1)
        for i5=1:size(Kr2,2)
            K2(Ielement(i1,i4),Ielement(i1,i5))=K2(Ielement(i1,i4),Ielement(i1,i5))+Kr2(i4,i5);
        end
    end
        density = density1;  %%
    Mr2=density*Me2*he^2;  %%%%%%%
    for i44=1:size(Mr2,1)
        for i55=1:size(Mr2,2)
            M2(Ielement(i1,i44),Ielement(i1,i55))=M2(Ielement(i1,i44),Ielement(i1,i55))+Mr2(i44,i55);
        end
    end
end
Kg2=(K2-Cf^2*M2)/mu1;

% end