%% octahedral elongation and compression calculation code
% authors:
% date:
% Input: structure file in xyz format
% outputs: elongation and compression along three different axis (x,y,z).
%%
clear
clc
tic
%% ! Read input: elongation !  atomic positions
filename = '38_3_CsPbCl3_Pnma.xyz';   % Input
f = fopen(filename);

natom = str2num(fgetl(f));
system = fgetl(f);

%% filename = system;
output = fopen('output_38_3.txt','a+');

poscart = [];
for j = 1:natom
    atom = fgetl(f);
    poscart = [poscart;str2double(atom(5:16)),str2double(atom(17:28)),str2double(atom(29:end)),j];
end
%% plot the original octahedra from xyz file.
figure(1)
hold on
scatter3(poscart(:,1),poscart(:,2),poscart(:,3),'cyan','Filled')

%% axis properties
grid minor
ax = gca;
ax.Box = 'on';
ax.LineWidth = 1;
ax.XLabel.String = 'x';
ax.YLabel.String = 'y';
ax.ZLabel.String = 'z';
%% find central atoms
sn = 38;        % compound serial no.
% Inputs: read the description
th = 0*ones(1,75);   % zeros(1,nooct)   %[0 45 0 45]
fi = 0*ones(1,75);   % zeros(1,nooct)   %[0 45 0 45]
si = 0*ones(1,75);   % zeros(1,nooct)   %[0 45 0 45]

n1 = 1;        % flag
for k = [1:4]   %% set index of central atoms from xyz file
    catom = poscart(k,:);

    %% measure bond length
    cenoct = catom(1:end);

    len1 = [];
    for i = 1:natom
        p = norm(cenoct(1:3) - poscart(i,1:end-1));
        len1 = [len1;p,i];
    end
    %% Design of octahedra
    p1 =  sortrows(len1);
    p2 = p1(1:7,:);

    posoct = [];
    for i = 2:7
        posoct = [posoct;poscart(p2(i,end),1:end)];
    end
    %% rotation of octahedra
    for i = 1:6
        m1 = -cenoct(1:3)+posoct(i,1:3);
        posoct1(i,1:3) = rotation_mat_3d(m1,th(k),fi(k),si(k));
        posoct(i,1:3) = cenoct(1:3)+posoct1(i,1:3);
    end
    %% mechanism to get x y z atoms
    diff = [];
    for i1 = 1:6
        diff = [diff;cenoct(1:3)-posoct(i1,1:3),i1];
    end
    minz = sortrows(diff,3);
    d56(1,:) = posoct(minz(1,end),:);
    d56(2,:) = posoct(minz(end,end),:);

    for i2 = 1:length(diff)
        if minz(1,end) == diff(i2,end)
            diff(i2,:) = [];
            break
        end
    end
    for i2 = 1:length(diff)
        if minz(end,end) == diff(i2,end)
            diff(i2,:) = [];
            break
        end
    end
    miny = sortrows(diff,2);
    d34(1,:) = posoct(miny(1,end),:);
    d34(2,:) = posoct(miny(end,end),:);
    for i2 = 1:length(diff)
        if miny(1,end) == diff(i2,end)
            diff(i2,:) = [];
            break
        end
    end
    for i2 = 1:length(diff)
        if miny(end,end) == diff(i2,end)
            diff(i2,:) = [];
            break
        end
    end
    minx = sortrows(diff,1);

    d12(1,:) = posoct(minx(1,end),:);
    d12(2,:) = posoct(minx(end,end),:);

%% Plotting of Octahedra after rotation and tilting
    figure(1)
    hold on
    scatter3(d12(:,1),d12(:,2),d12(:,3),200,'r','Filled')
    scatter3(d34(:,1),d34(:,2),d34(:,3),200,'g','Filled')
    scatter3(d56(:,1),d56(:,2),d56(:,3),200,'b','Filled')

    plot3([cenoct(1),d12(1,1)],[cenoct(2),d12(1,2)],[cenoct(3),d12(1,3)],'k')
    plot3([cenoct(1),d12(2,1)],[cenoct(2),d12(2,2)],[cenoct(3),d12(2,3)],'k')
    
    plot3([cenoct(1),d34(1,1)],[cenoct(2),d34(1,2)],[cenoct(3),d34(1,3)],'k')
    plot3([cenoct(1),d34(2,1)],[cenoct(2),d34(2,2)],[cenoct(3),d34(2,3)],'k')
    
    plot3([cenoct(1),d56(1,1)],[cenoct(2),d56(1,2)],[cenoct(3),d56(1,3)],'k')
    plot3([cenoct(1),d56(2,1)],[cenoct(2),d56(2,2)],[cenoct(3),d56(2,3)],'k')

    l0 = mean(p2(2:7));

    %% Octahedral distortion calculation
    lambda_oct12 = 0;

    for i = 1:2
        l = norm(cenoct(1:3) - d12(i,1:end-1));
        lambda_oct12 = lambda_oct12 + (l/l0)^2/2;
    end

    lambda_oct34 = 0;

    for i = 1:2
        l = norm(cenoct(1:3) - d34(i,1:end-1));
        lambda_oct34 = lambda_oct34 + (l/l0)^2/2;
    end

    lambda_oct56 = 0;

    for i = 1:2
        l = norm(cenoct(1:3) - d56(i,1:end-1));
        lambda_oct56 = lambda_oct56 + (l/l0)^2/2;
    end
%% distotion parameters plot
    figure(2)
    hold on
    scatter(sn+n1*0.1,lambda_oct12,200,'s','Filled','MarkerEdgeColor','r','MarkerFaceColor','r')
    scatter(sn+n1*0.1,lambda_oct34,200,'p','Filled','MarkerEdgeColor','g','MarkerFaceColor','g')
    scatter(sn+n1*0.1,lambda_oct56,200,'^','Filled','MarkerEdgeColor','b','MarkerFaceColor','b')
%% save data in a file
    fprintf(output,'%i %f %f %f %s %s\n',sn+n1*0.1,lambda_oct12,lambda_oct34,lambda_oct56,system,filename);
    n1 = n1+1;
end

fclose(output);

%% axis properties
ax = gca;
ax.Box = 'on';
ax.LineWidth = 2;
ax.XLabel.String = 'No. of inequivalent octahedra';
ax.YLabel.String = 'Octahedral elongation and compression';
ax.XTickLabel = [1,2,3,4];
ax.YLim = [0.997, 1.002];

toc
%% end of the code