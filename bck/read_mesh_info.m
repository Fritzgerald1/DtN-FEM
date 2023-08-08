function [Nnode, Nelement, Coordinate, Ielement, dx, dy] = read_mesh_info()
meshInfo = importfile("mesh.dat");

x = meshInfo.Var1(meshInfo.Type == "GRID");
y = meshInfo.Var2(meshInfo.Type == "GRID");
Coordinate = [x, y];
Nnode = size(Coordinate,1);

[~,~,dx] = find(diff(x),1);
[~,~,dy] = find(diff(y),1);

Ielement = [meshInfo.Var1(meshInfo.Type == "CQUAD8"),	meshInfo.Var2(meshInfo.Type == "CQUAD8"),meshInfo.Var3(meshInfo.Type == "CQUAD8"),meshInfo.Var4(meshInfo.Type == "CQUAD8"),...
	meshInfo.Var5(meshInfo.Type == "CQUAD8"),meshInfo.Var6(meshInfo.Type == "CQUAD8"),meshInfo.Var7(meshInfo.Type == "CQUAD8"),meshInfo.Var8(meshInfo.Type == "CQUAD8"),];
Nelement = size(Ielement,1);

end