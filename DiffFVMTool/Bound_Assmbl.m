function[Belements] = Bound_Assmbl(mesh,Belements,Acoef)
% Boundary assembly function with already some data stored
for i = 1:mesh.NEB
    if mesh.ENT{Belements(i).elset,2}{1,1} == "Facades"
        alpha = 0.01;
        Belements(i).Abn = Acoef(alpha);
%% Syntax of the boundary assembling :
% Below must be added all the different materials present inside the model
% with the following syntax (some examples are already coded here suppress them if needed):
%
%  elseif mesh.ENT{Belements(i).elset,2}{1,1} == "<MaterialNameInSketchup>"
%     alpha = <AbsorptionCoeffValue(ex: 0.1)>;
%     Belements(i).Abn = Acoef(alpha);    
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Airbound"
        alpha = 1;
        Belements(i).Abn = Acoef(alpha); 
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Bichette"
        alpha = 0.2;
        Belements(i).Abn = Acoef(alpha);
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Ground"
        alpha = 0.03;
        Belements(i).Abn = Acoef(alpha);
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Isolant"
        alpha = 0.39;
        Belements(i).Abn = Acoef(alpha);
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Moquette"
        alpha = 0.01;
        Belements(i).Abn = Acoef(alpha);
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Plafond"
        alpha = 0.35;
        Belements(i).Abn = Acoef(alpha);
    elseif mesh.ENT{Belements(i).elset,2}{1,1} == "Hard wall"
        alpha = 0.017;
        Belements(i).Abn = Acoef(alpha);
    elseif isempty(Belements(i).Abn) 
        text(0,0.5,'Err: No Corresponding BC material in database','fontsize',16,'color','red')
        return
    end
end