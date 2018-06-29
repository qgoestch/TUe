function [Chi_el,Lelv,Chi_elV,Lelb,Chi_elB] = Orth_idx(mesh,Velements,Cell_center)
Chi_el = zeros(mesh.NEV,1);
Chi_elV = zeros(0,1);
Chi_elB = zeros(0,1);
Lelv = zeros(0,1);
Lelb = zeros(0,1);
for i = 1:mesh.NEV
    Nc0 = mesh.POS(Velements(i).nodes(1),:);
    Nc1 = mesh.POS(Velements(i).nodes(2),:);
    Nc2 = mesh.POS(Velements(i).nodes(3),:);
    Nc3 = mesh.POS(Velements(i).nodes(4),:); 
    face_center = zeros(mesh.Eltype,3);
    boundf = 0;
    for j = 1:3
        face_center(1,j) = (Nc0(j)+Nc1(j)+Nc2(j))/3;
        face_center(2,j) = (Nc1(j)+Nc2(j)+Nc3(j))/3;
        face_center(3,j) = (Nc0(j)+Nc1(j)+Nc3(j))/3;
        face_center(4,j) = (Nc0(j)+Nc2(j)+Nc3(j))/3;
    end
    A_vec = zeros(3,mesh.Eltype);
    A_vec(:,1) = Area_vec(Nc0,Nc1,Nc2);
    A_vec(:,2) = -Area_vec(Nc1,Nc2,Nc3);
    A_vec(:,3) = -Area_vec(Nc0,Nc1,Nc3);
    A_vec(:,4) = Area_vec(Nc0,Nc2,Nc3);
    f_vec = zeros(3,mesh.Eltype);
    c_vec = zeros(3,mesh.Eltype);
    for k = 1:mesh.Eltype
        f_vec(:,k) = (face_center(k,:) - Cell_center(i,:))';
        if Velements(i).neighbours(k) == 0
            c_vec(:,k) = f_vec(:,k);
            boundf = 1;
        else
            c_vec(:,k) = (Cell_center(Velements(i).neighbours(k),:) - Cell_center(i,:))';
        end
    end
    chi_f = zeros(mesh.Eltype,1);
    for l = 1:mesh.Eltype
        chi_f(l) = min(dot(A_vec(:,l),f_vec(:,l))/(norm(A_vec(:,l))*norm(f_vec(:,l))),...
            dot(A_vec(:,l),c_vec(:,l))/(norm(A_vec(:,l))*norm(c_vec(:,l))));
    end
    if boundf == 1
        Chi_elB = [Chi_elB ; min(chi_f)];
        Lelb = [Lelb , i];
    else
        Chi_elV = [Chi_elV ; min(chi_f)];
        Lelv = [Lelv , i];
    end
    Chi_el(i) = min(chi_f);
end
end