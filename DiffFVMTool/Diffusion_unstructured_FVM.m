% Code developped by Quentin Goestchel for the analysis of the FVM method
% on diffusion equation.
%--- CAUTION !---% The information below are important:
% To run the code below, the file newid.m must be added in the directory
% MATLAB_<version>/toolbox/matlab/uitools and the gptoolbox folder must be
% in the current folder.
% For Mac (Unix) user, Matlab must be opened from the
% Terminal in order to enable the command line control of Gmsh
% For Windows user, Gmsh.exe must be in the work folder
% To download Gmsh: http://gmsh.info/
%% Beginning 
clear variables
close all
close all hidden
clc
% Import the gptoolbox, download it and change the path 
gp_subdirs = split(genpath(strcat(pwd,'/gptoolbox/')),':');
addpath(strjoin(gp_subdirs(~contains(gp_subdirs,'.git')),':'));
savepath
% Configuration of the plot interface 
set(gcf, 'Position', get(0, 'Screensize'));
ax = subplot(3, 2, 1);
set (ax,'visible','off')
figure(1)
subplot(3,2,1)
text(0,1,'Initialising modelisation...','fontsize',14);
%% Mesh importation

% Manual input of the filename, caracteristic length of the mesh,...
% Boundary conditions, importations, recording time and source emition duration 
prompt = {'Enter the name of the .geo file [Without extension]: ',...
    'Enter the caracteristic length of the mesh [in meter]: ',...
    'Enter the boundary conditions type desired [1:Sabine,2:Eyring,3:Modified] ',...
    'Is the mesh already imported? [Y/N](default value Y): ',...
    'Do you want a 3D plot of the mesh orthogonality? [Y/N](default value N): ',...
    'Enter recording time [in seconds]: ',...
    'Enter source signal duration [in seconds]: '};
defaultanswer = {'sphere','1','3','N','N','4','2'};
numlines = 1;
answer = newid(prompt,'Code inputs',numlines,defaultanswer);
options.Resize='on';
filename = num2str(answer{1});
lc = num2str(answer{2});
th = str2double(answer{3});
already_imported = num2str(answer{4});
plotidx = num2str(answer{5});
recording_time =  str2double(answer{6});
sourceon_time =  str2double(answer{7}); 
% To save time, the mesh importation can be stored in mesh_data.mat
if isempty(already_imported)
    str = 'Y';
end
if already_imported == 'N' % Condition if nothing has been stored before: 
    % The lines below command Gmsh to do the 3D meshing
    if lc == 'H'
        strimp = strcat({'gmsh -3 '},{pwd},{'/GeoModels/'},{filename},{'.geo'});
        system(strimp{1});
        lc = 2;
    else
        strimp = strcat({'gmsh -3 -clscale '},{lc},{' '},{pwd},{'/GeoModels/'},{filename},{'.geo'});
        system(strimp{1});
        lc = str2double(lc);
    end
    % The lines below open the .msh file, call the connectivity function
    % and save all the data in a .mat file
    msh_file = strcat(pwd,'/GeoModels/',filename,'.msh');
    [mesh,Velements,Belements,Cell_center,Vcell] = Input_mesh(msh_file,lc);
    % Saving of the mesh connectivity data:
    save('mesh_data.mat','mesh','Velements','Belements','Cell_center','Vcell')
else
    load('mesh_data.mat');
    lc = str2double(lc);
end
figure(1)
subplot(3,2,1)
text(0,0.9,'Volume meshed by Gmsh and connectivity process completed','fontsize',14');
%% Mesh orthogonality evaluation %%
% Call of the orthogonnality index computation function
[Chi_el,Lelv,Chi_elV,Lelb,Chi_elB] = Orth_idx(mesh,Velements,Cell_center);
% Sorting of the axis and data for volume and boundary related elements
X = [Lelv , Lelb];
Y = [Chi_elV ; Chi_elB];
figure(1)
subplot(3,2,2)
stem(Lelv,Chi_elV)
hold on
stem(Lelb,Chi_elB)
legend({'Volume elements','Boundary related elements'},'FontSize',22)
title('Mesh orthogonal quality','fontsize',18)
hline = refline([0 mean(Y)]);
hline.Color = 'c';
hline2 = refline( [0 median(Y)]);
hline2.Color = 'g';
set([hline hline2],'LineWidth',2)
legend('Location','southwest');
legend({'In volume','Boundary related','Mean','Median'},'FontSize',12);
xlabel('Elements index','fontsize',18)
ylabel('Orthogonal quality','fontsize',18)
axis([0 mesh.NEV 0 1])

% figure(2)
% stem(Lelv,Chi_elV)
% hold on
% stem(Lelb,Chi_elB)
% title('Mesh orthogonal quality','fontsize',18')
% hline = refline([0 mean(Y)]);
% hline.Color = 'c';
% hline2 = refline( [0 median(Y)]);
% hline2.Color = 'g';
% set([hline hline2],'LineWidth',2)
% legend('Location','southwest');
% legend([hline hline2],{'Mean','Median'},'FontSize',12);
% xlabel('Elements index','fontsize',18)
% ylabel('Orthogonal quality','fontsize',18)
% axis([0 mesh.NEV 0 1])


if isempty(plotidx)
    str = 'N';
end
if plotidx == 'Y'
    T = zeros(mesh.NEV,4);
    for i = 1:mesh.NEV
        T(i,:) = Velements(i).nodes;
    end
    figure(1)
    subplot(3,2,4)
    tetramesh(T,mesh.POS,Chi_el,'FaceAlpha',0.4)
    title('3D Mesh orthogonal quality','fontsize',18)
    colormap hot
    light;
    lighting phong;
    camlight('left');
    axis equal
    caxis([ 0 1 ]) %Set the scale in order to have an interpretable graphic result
    scale = colorbar; % Create a colorbar
    scale.Label.String = 'Orthogonality index';
    scale.Label.FontSize = 14;
    drawnow
    
%     figure(3)
%     tetramesh(T,mesh.POS,Chi_el,'FaceAlpha',0.4)
%     title('3D Mesh orthogonal quality','fontsize',18')
%     colormap hot
%     light;
%     lighting phong;
%     camlight('left');
%     axis equal
%     caxis([ 0 1 ]) %Set the scale in order to have an interpretable graphic result
%     scale = colorbar; % Create a colorbar
end
figure(1)
subplot(3, 2, 1)
text(0,0.8,['(Mean ; median) orthogonality index = ( ',num2str(mean(Y)),' ; ',num2str(median(Y)),' )'],'fontsize',14);
chig = mean(Y);
%% Data input %%%
c = 343; %Sound particle velocity [m.s^-1]
rho = 1.21; %Air density [Kg.m^-3] at 20°C
% Below, temporal discretization (time step) [s] which is set to fulfill
% Adapted empirical Navarro criterion
dt = sqrt(lc/2*(1*10^-8));
fm = 1/dt;

% Extreme limit of the consistance condition:
if lc/(2*dt) < c  
    figure(1)
    subplot(3, 2, 1)
    text(0,0.6,'Error: lc/(2dt) < to c.','fontsize',20,'color','red')
    text(0,0.4,'The propagation of the numerical scheme is slower than sound speed','fontsize',14);
    return
end 

figure(1)
subplot(3, 2, 1)
text(0,0.7,['"Adapted Navarro" empirical criterion = ',num2str(2*dt^2/lc)],'fontsize',14');

%% Boundary condition assembly %% 

% th type of constant A on boundary conditions. options:
%  1. sabine A
%  2. eyring A
%  3. modified Xiang A
if th==1
    Acoef = @(alpha) (c*alpha)/4; % Sabine
elseif th==2
    Acoef = @(alpha) (c*(-log(1-alpha)))/4; % Eyring
elseif th==3
    Acoef = @(alpha) (c* alpha)/(2*(2-alpha)); % Xiang
else
    close all hidden
    disp('Error : Wrong boundary condition input')
    return
end

if mesh.ENT{1,1} == 0 % Condition for meshes without any BC data
    figure(1)
    subplot(3, 2, 1)
    text(0,0.6,{'No mesh BC absorption data: BC are set to alpha = 0.1'},'fontsize',14);
    alpha = 0.1;
    Abn = Acoef(alpha);
else
    % All the boundary assembly process is made in the following function
    % To add a material, open the source code of the Bound_Assmbl function.
    Belements = Bound_Assmbl(mesh,Belements,Acoef);  
end
%% Geometrical parameters computation %%
% Total boundary surface
S = 0;
Serr = 0;
% Geometrical parameters for diffusive fluxes
% Sum of the 4 connectivity face area over the distance to the neighbour per
% element:
fel = zeros(mesh.NEV,1); % For the volume elements 
fb = zeros(mesh.NEV,1); % Boundary vector (For the boundary faces)

% Global matrix of the faces areas over the distance to the neighbour:
Fmat = sparse(mesh.NEV,mesh.NEV);
dtMat = sparse(mesh.NEV,mesh.NEV);
% Fmat = zeros(mesh.NEV,mesh.NEV);
wtBC = waitbar(0,'Attribution of the boundary conditions');
BCcount = 0;
for i = 1:mesh.NEV
    % Extraction of the nodes coordinates 
  Nc1 = mesh.POS(Velements(i).nodes(1),:);
  Nc2 = mesh.POS(Velements(i).nodes(2),:);
  Nc3 = mesh.POS(Velements(i).nodes(3),:);
  Nc4 = mesh.POS(Velements(i).nodes(4),:);
  
  % Face 1 %%%%%%%%
  if Velements(i).neighbours(1) == 0 % Condition to recognize a boundary
      % Addition to the global boundary surface
      S = S + norm(cross((Nc3-Nc1),(Nc2-Nc1)))/2;
      BCcount = BCcount + 1;
      if mesh.ENT{1,1} == 0 % Condition for mesh without any BC data
          fb(i) = norm(cross((Nc3-Nc1),(Nc2-Nc1)))*Abn/2;
      else
          for j = 1:mesh.NEB
              % Enable to find which triangle element has 3 nodes in common
              % with the considered volume element
              if sort(Velements(i).nodes(1:3)) == sort(Belements(j).nodes)
                  fb(i) = norm(cross((Nc3-Nc1),(Nc2-Nc1)))*Belements(j).Abn/2 ;
              end
          end
      end

  else
      Fmat(i,Velements(i).neighbours(1)) = norm(cross((Nc3-Nc1),(Nc2-Nc1)))/...
          (2*norm(Cell_center(i,:)-Cell_center(Velements(i).neighbours(1),:)));
      fel(i) = Fmat(i,Velements(i).neighbours(1));  
  end
  
  % Face 2 %%%%%%%%
  if Velements(i).neighbours(2) == 0
      S = S + norm(cross((Nc4-Nc2),(Nc3-Nc2)))/2;
      BCcount = BCcount + 1;
      if mesh.ENT{1,1} == 0
          fb(i) = fb(i) + norm(cross((Nc4-Nc2),(Nc3-Nc2)))*Abn/2;
      else
          for j = 1:mesh.NEB
              if sort(Velements(i).nodes(2:4)) == sort(Belements(j).nodes)
                  fb(i) = fb(i) + norm(cross((Nc4-Nc2),(Nc3-Nc2)))*Belements(j).Abn/2;
              end
          end
      end
  else
      Fmat(i,Velements(i).neighbours(2)) = norm(cross((Nc4-Nc2),(Nc3-Nc2)))/...
          (2*norm(Cell_center(i,:)-Cell_center(Velements(i).neighbours(2),:)));
      fel(i) = fel(i) + Fmat(i,Velements(i).neighbours(2));
  end
  
  % Face 3 %%%%%%%%
  if Velements(i).neighbours(3) == 0   
      S = S + norm(cross((Nc4-Nc1),(Nc2-Nc1)))/2;
      BCcount = BCcount + 1;
      if mesh.ENT{1,1} == 0
          fb(i) = fb(i) + norm(cross((Nc4-Nc1),(Nc2-Nc1)))*Abn/2;
      else
          for j = 1:mesh.NEB
              if sort(Velements(i).nodes([1:2 4])) == sort(Belements(j).nodes)
                  fb(i) = fb(i) + norm(cross((Nc4-Nc1),(Nc2-Nc1)))*Belements(j).Abn/2;
              end
          end
      end
  else
      Fmat(i,Velements(i).neighbours(3)) = norm(cross((Nc4-Nc1),(Nc2-Nc1)))/...
          (2*norm(Cell_center(i,:)-Cell_center(Velements(i).neighbours(3),:)));
      fel(i) = fel(i) + Fmat(i,Velements(i).neighbours(3));
  end
  
  % Face 4 %%%%%%%%
  if Velements(i).neighbours(4) == 0
      S = S + norm(cross((Nc4-Nc1),(Nc3-Nc1)))/2;
      BCcount = BCcount + 1;
      if mesh.ENT{1,1} == 0
          fb(i) = fb(i) + norm(cross((Nc4-Nc1),(Nc3-Nc1)))*Abn/2;
      else
          for j = 1:mesh.NEB
              if sort(Velements(i).nodes([1 3:4])) == sort(Belements(j).nodes)
                  fb(i) = fb(i) + norm(cross((Nc4-Nc1),(Nc3-Nc1)))*Belements(j).Abn/2;
              end
          end
      end     
  else
      Fmat(i,Velements(i).neighbours(4)) = norm(cross((Nc4-Nc1),(Nc3-Nc1)))/...
          (2*norm(Cell_center(i,:)-Cell_center(Velements(i).neighbours(4),:)));
      fel(i) = fel(i) + Fmat(i,Velements(i).neighbours(4));
  end
  waitbar(i / mesh.NEV);
end
close(wtBC)
if mesh.ENT{1,1} ~= 0
    figure(1)
    subplot(3, 2, 1)
    text(0,0.6,'BC absorption coefficient attributed to faces','fontsize',14);
end
% Considering that the matrix Fmat will be full of zeros, Matlab sparse
% operation enable to save storage and time. 
% Fmat = sparse(Fmat);

% The line below is possible to use only if a Nvidia GPU is installed on
% the computer
% Fmat = gpuArray(Fmat);

%% Sources and receivers definition (excitation with Interrupted source) %%

% The part below locate the center of the mesh no matter of the coordinate
% system used 
dLsx = round(abs(max(mesh.POS(:,1))) - abs(min(mesh.POS(:,1))));
dLsy = round(abs(max(mesh.POS(:,2))) - abs(min(mesh.POS(:,2))));
dLsz = round(abs(max(mesh.POS(:,3))) - abs(min(mesh.POS(:,3))));
dLx = round(abs(max(mesh.POS(:,1)) - min(mesh.POS(:,1))));
dLy = round(abs(max(mesh.POS(:,2)) - min(mesh.POS(:,2))));
dLz = round(abs(max(mesh.POS(:,3)) - min(mesh.POS(:,3))));

% Determination of the source point and receivers list
% Here source point set in the center of the mesh
if dLsx == dLx
    spoint = [floor(dLx/2),floor(dLy/2),floor(dLz/2)];
else
    spoint = [dLsx,dLsy,floor(dLsz/2)];
end
% Receivers step discretisation
del = 15;
% Building of the receiver list
Rpel = zeros(1,del);
for i = 1:del
    Rpel(i) = floor(i*mesh.NEV/del);
end

% "Middle value of the receiver list" 
rpel = Rpel(ceil(del/2));
% rpel2 = ceil(mesh.NEV/2);
snode = 1;
% Source element list initialisation for future interpolation
sel_list =[];
lsum = 0; 
li = [];
Vs = 0;
for j = 1:mesh.NN
    if norm(mesh.POS(j,:)-spoint) < 1.5*lc
        if norm(mesh.POS(j,:) - spoint) < norm(mesh.POS(snode,:) - spoint)
            snode = j;
        end
    end
end
for i = 1:mesh.NEV
    % Loop for the receiver element attribution when wanted to be close to
    % a position
    % Here position set to [dLx/4,dLy/4,dLz/4]
%     if norm(Cell_center(i,:)-[dLx/4,dLy/4,dLz/4]) < lc
%         if norm(Cell_center(i,:)-[dLx/4,dLy/4,dLz/4]) < norm(Cell_center(rpel,:)-[dLx/4,dLy/4,dLz/4])
%             rpel = i;
%         end
%     end
%     if norm(Cell_center(i,:)-[dLx,dLy,dLz]) < lc
%         if norm(Cell_center(i,:)-[dLx,dLy,dLz]) < norm(Cell_center(rpel2,:)-[dLx,dLy,dLz])
%             rpel2 = i;
%         end
%     end
    if  ismember(snode,Velements(i).nodes)
        sel_list = [sel_list;i]; %#ok<*AGROW>
        lsum = lsum + norm(Cell_center(i) - mesh.POS(snode,:));
        li = [li;norm(Cell_center(i) - mesh.POS(snode,:))];
        Vs = Vs + Vcell(i);
    end
end

Ws = 10^-2;   % Point power [W]
w1 = Ws/Vs;
% Duration of the simulation [s]
recording_steps = ceil(recording_time/dt);
% Interrupted source, time before switching of [s]
sourceon_steps = ceil(sourceon_time/dt);

s1 = w1.*ones(1,sourceon_steps);               
source = [s1 zeros(1,recording_steps-sourceon_steps)];%for source on time, interrupted source

s = zeros(mesh.NEV,1);
% Interpolation of the source point to the surronding elements. 
s(sel_list) = li*source(1)/lsum;
%% Diffusion parameters %%

V = sum(Vcell); %Total volume of the mesh 
if strcmp(filename,'sphere') || strcmp(filename,'NewTransfSphere') || strcmp(filename,'sphereTransfinite') 
    Vr = 4/3*pi*(dLx/2)^3;
    Sr = 4*pi*(dLx/2)^2;
    mfp = 4*Vr/Sr;    % Mean-free-path for a sphere
elseif strcmp(filename,'cube') || strcmp(filename,'Cube_attractors')
    mfp = 4*V/S;    % Mean-free-path for a sphere
else
    mfp = 4*V/S;
end
mft = mfp/c;    % Mean-free-time
D = mfp*c/3;    % Diffusion coefficient
m_atm = 1e-5; % Atmospheric attenuation [m^-1]

%% Time loop initialisation %%

w_new = zeros(mesh.NEV,1);
w = w_new;
w_old = w;
gamma_zero = dt*((D.*fel)+fb)./Vcell;
 
tic
%% Time loop %%
cpt10 = 0;
cpt40 = 0;
t_10 = 0;
t_40 = 0;
w_rec = zeros(del,recording_steps); 
wtTL = waitbar(0,'Numerical scheme processing');
w_data = cell(recording_steps,1);
Res = cell(recording_steps,1); %Residual value of the numerical scheme
for t_steps = 1:recording_steps
    time = t_steps*dt;
    % Scheme implementation:
    w_new = (w_old.*(1-gamma_zero) - c*m_atm*2*dt*w +2*dt*s +...
        2*dt*D*(Fmat*w)./Vcell)./(1+gamma_zero);

    Res{t_steps} = w_new.*(1+gamma_zero) - w_old.*(1-gamma_zero) +...
        c*m_atm*2*dt*w -2*dt*s - 2*dt*D*(Fmat*w)./Vcell;

    w_old = w;
    w = w_new;
    s(sel_list) = li*source(t_steps)/lsum;
    %The for loop below is slowing down the computation depending on
    % del, the discretisation of the receivers
    for rec_step = 1:del
        w_rec(rec_step,t_steps) = w_new(Rpel(rec_step));
    end
    w_data{t_steps} = w_new; %Comment this line for experiment on efficiency (Storing costs!)
    if time >= sourceon_time 
        if w_new(rpel) < w_data{sourceon_steps}(rpel)*10^(-1) && cpt10 == 0
            t_10 = time-sourceon_time;
            cpt10 = 1;
        end
        if w_new(rpel) < w_data{sourceon_steps}(rpel)*10^(-4) && cpt40 == 0
            t_40 = time-sourceon_time;
            cpt40 = 1;
        end
    end
    waitbar(t_steps / recording_steps);
end
toc
close(wtTL)
%% Extraction of the results %
% Numeric residual computation along time steps
Res_t = zeros(1,recording_steps);
for t_step = 1:recording_steps
    Res_t(t_step) = norm(Res{t_step});
end
if strcmp(filename,'sphere') || strcmp(filename,'NewTransfSphere')|| strcmp(filename,'sphereTransfinite') 
    figure(1)
    subplot(3, 2, 1)
    text(0,0.5,['Surface approximation error: ',num2str(100*abs(Sr-S)/Sr),' %'],'fontsize',14);
    text(0,0.4,['Volume approximation error: ',num2str(100*abs(Vr-V)/Vr),' %'],'fontsize',14);
    eS = 100*abs(Sr-S)/Sr;
    eV = 100*abs(Vr-V)/Vr;
end
% Safety if in case of too short recording time 
if t_10 == 0 || t_40 == 0
    figure(1)
    subplot(3, 2, 1)
    text(0,0.3,'Reverberation time RT30 > recording time','fontsize',14');
else
    RT30 = 2*(t_40-t_10);
    figure(1)
    subplot(3, 2, 1)
    text(0,0.3,['Reverberation time RT30 = ',num2str(2*(t_40-t_10)),' s'],'fontsize',14);
    if strcmp(filename,'sphere') || strcmp(filename,'NewTransfSphere') || strcmp(filename,'sphereTransfinite') 
        RTSab = 6*log(10)*4*Vr/(c*Sr*alpha);
        RTEyrin = 6*log(10)*4*Vr/(-c*Sr*log(1-alpha));
        RTXiang = 6*log(10)*4*Vr*(2-alpha)/(c*Sr*2*alpha);
        eSab = abs(((RTSab-RT30)/RTSab)*100);
        eEyr = abs(((RTEyrin-RT30)/RTEyrin)*100);
        eXia = abs(((RTXiang-RT30)/RTXiang)*100);
        eFem = abs(((2.58-RT30)/2.58)*100);
        text(0,0.15,['error_{sab} : ',num2str(abs(((RTSab-RT30)/RTSab)*100)),' % ;',...
            'error_{Eyr} : ',num2str(abs(((RTEyrin-RT30)/RTEyrin)*100)),' % '],'fontsize',14);
        text(0,0,['error_{Xiang} : ',num2str(abs(((RTXiang-RT30)/RTXiang)*100)),' % ;',...
            'error_{FEM} : ',num2str(abs(((2.58-RT30)/2.58)*100)),' %'],'fontsize',14); 
    end
    if strcmp(filename,'cube') || strcmp(filename,'Cube_attractors')
        RTSab = 6*log(10)*4*V/(c*S*alpha);
        RTEyrin = 6*log(10)*4*V/(-c*S*log(1-alpha));
        RTXiang = 6*log(10)*4*V*(2-alpha)/(c*S*2*alpha);
        text(0,0.2,['error_{sab} : ',num2str(abs(((RTSab-RT30)/RTSab)*100)),' % ;',...
            'error_{Eyr} : ',num2str(abs(((RTEyrin-RT30)/RTEyrin)*100)),' % ;',...
            'error_{Xiang} : ',num2str(abs(((RTXiang-RT30)/RTXiang)*100)),' % ;',...
            'error_{FEM} : ',num2str(abs(((2.6-RT30)/2.6)*100)),' %'],'fontsize',18);
    end
end
%line(2*mft,1);

% figure(1)
% subplot(3,2,6)
% plot((1:t_steps)*dt,w_rec,(1:t_steps)*dt,w_rec2)
% title('w evolution at the receivers','fontsize',10','fontsize',18')
% xlabel('time in s','fontsize',18')
% ylabel('Sound energy density at (Nxc/4,Nyc/4,Nzc/4)','fontsize',12')
for j = 1:del
figure(1)
subplot(3,2,6)
plot((1:t_steps)*dt,10*log10(abs(w_rec(j,:)*rho*c^2/(2e-5)^2)))
hold on
end
axis([0 recording_time 0 max(10*log10(abs(w_rec(j,:)*rho*c^2/(2e-5)^2)))+10]);
title('SPL evolution at the receivers','fontsize',10','fontsize',18')
xlabel('time in s','fontsize',18')
ylabel('Sound pressure level [dB] at receiver elements','fontsize',12)
hold off

%% 3D time-depending plot %%
prompt = 'Do you want a 3D plot? [Y/N](default value N): ';
plotmesh = newid(prompt,'3D plot choice');
plotmesh = num2str(plotmesh{1});
if isempty(plotmesh)
    str = 'N';
end
if plotmesh == 'Y'
    T = zeros(mesh.NEV,4);
    for i = 1:mesh.NEV
        T(i,:) = Velements(i).nodes;
    end
    fig=figure('Renderer','zbuffer');
    y=[];
    for t_steps = 1:100:recording_steps
        interp = scatteredInterpolant(Cell_center,w_data{t_steps});
        wint = interp(mesh.POS(:,1),mesh.POS(:,2),mesh.POS(:,3));
        F = boundary_faces(T);
        figure(1)
        subplot(25,2,(21:2:47))
        t = tsurf(F,mesh.POS,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.0,'EdgeAlpha',0.0);
        hold on;
        s = tsurf([1 1 1],mesh.POS,'EdgeColor','none',fphong);
        hold off
        title('Sound density energy evolution in a plan','fontsize',18)
        caxis([0 max(w_data{sourceon_steps})]);
        axis equal;
        scale = colorbar; % Create a colorbar
        scale.Label.String = 'Sound energy density in J.m^{-3}';
        [U,G,J,BC] = slice_tets(mesh.POS,T,[0 0 1 (-spoint(3)+0.2)]);
        set(s,'Vertices',U,'Faces',G,'CData',BC*wint);
        colormap hot
        light;
        lighting phong;
        drawnow;
        subplot(25,2,49);
        hold on 
        y=[y,1];
        bar(y);
        set(gca,'xtick',0:100:recording_steps,'ytick',[]);
        axis([0.5 recording_steps 0 1]);
        Fprocess(t_steps)=getframe(fig);
        hold off
        axis off
    end
else
    t_vec = (1:t_steps)*dt;
    figure(1)
    subplot(25,2,(21:2:47))
    plot(t_vec,Res_t);
    title('Numeric residual evolution','fontsize',18)
    xlabel('time in s','fontsize',18)
    ylabel('Residual value','fontsize',12)
    axis([0 recording_time 0 inf]);
    legend('Norm(r(el))')
    Nel = mesh.NEV;
%     AnalysisFile = strcat({'FullTransfAnalysis_'},{num2str(mesh.NEV)},{'.mat'});
%     save(AnalysisFile{1},'t_vec','Res_t','Nel','eFem','eXia','eS','eV')
end