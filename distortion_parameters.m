%% distortion parameters calculation code
% authors:
% date:
% Input: structure file in xyz format
% outputs: sigma_variance, octahedral_distortion, theta_ab
%%
clear
clc
tic
%% ! Input : example !  atomic positions form .xyz file
filename = '38_3_CsPbCl3_Pnma.xyz';   % Input
f = fopen(filename);

natom = str2num(fgetl(f));
system = fgetl(f);

%% output filename: example
output = fopen('output_38_3_dist_parameters.txt','a+');

poscart = [];
for j = 1:natom
    atom = fgetl(f);
    poscart = [poscart;str2double(atom(5:16)),str2double(atom(17:28)),str2double(atom(29:end)),j];
end
%% original structure plot
figure(1)
hold on
scatter3(poscart(:,1),poscart(:,2),poscart(:,3),'blue','Filled')

grid minor
ax = gca;

ax.Box = 'on';
ax.LineWidth = 1;
ax.XLabel.String = 'x';
ax.YLabel.String = 'y';
ax.ZLabel.String = 'z';
%% calculations of distortion parameters: (octahedral bond length distortion, octahedral angle variance, tilting and rotation angle)
%% find central atoms
noct = 4;
n1 = 1;        % flag
for k = [1:4]   %% set index of central atoms from xyz file
    % central atom
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
%% bond length plotting from the B-atom to six nearest X-atom to complete representation of the octahedron
    plot3([cenoct(1),posoct(1,1)],[cenoct(2),posoct(1,2)],[cenoct(3),posoct(1,3)],'k')
    plot3([cenoct(1),posoct(2,1)],[cenoct(2),posoct(2,2)],[cenoct(3),posoct(2,3)],'k')

    plot3([cenoct(1),posoct(3,1)],[cenoct(2),posoct(3,2)],[cenoct(3),posoct(3,3)],'k')
    plot3([cenoct(1),posoct(4,1)],[cenoct(2),posoct(4,2)],[cenoct(3),posoct(4,3)],'k')

    plot3([cenoct(1),posoct(5,1)],[cenoct(2),posoct(5,2)],[cenoct(3),posoct(5,3)],'k')
    plot3([cenoct(1),posoct(6,1)],[cenoct(2),posoct(6,2)],[cenoct(3),posoct(6,3)],'k')

%% bond length distortion
l0 = mean(p2(2:7));
delta_d = 0;
for i = 1:6
    l = norm(cenoct(1:3) - posoct(i,1:3));
delta_d = delta_d + ((l-l0)/l0)^2/6;
end
octahedral_distortion(n1) = delta_d;

    %%
    catom1(n1,:) = cenoct(1:3);
    octaatom1(:,:,n1)= posoct(:,1:3);

    %% angle variance
    ang1 = [];

    for i = 1:6
        for j = i:6
            if i == j
                continue
            end
            norm(cenoct(1:3)-posoct(j,1:3))
            ang1 = [ang1;acosd(dot(cenoct(1:3)-posoct(i,1:3),cenoct(1:3)-posoct(j,1:3))/(norm(cenoct(1:3)-posoct(i,1:3))*norm(cenoct(1:3)-posoct(j,1:3))))];
        end
    end

    ang11 = sort(ang1,'ascend');
    %
    stheta_oct1 = 0;
    for i = 1:12
        stheta_oct1 = stheta_oct1+(ang11(i)-90)^2/11;
    end

    fprintf(output,'sigma variance = %f \n',stheta_oct1);
    n1 = n1+1;
end

%% tilting and rotation angle
k3 = 1;
for i = 1:noct
    for  j = i:noct
        if i == j
            continue
        end
        c1 = catom1(i,:);
        c2 = catom1(j,:);

        o1 = octaatom1(:,:,i);
        o2 = octaatom1(:,:,j);

        for k1 = 1:6
            for k2 = 1:6
                l = o1(k1,:) == o2(k2,:);
                if l(1) == 1 && l(2) == 1 && l(3) == 1
                    theta_ab(k3) = acosd(dot(c1-o1(k1,1:3),c2-o2(k2,1:3))/(norm(c1-o1(k1,1:3))*norm(c2-o2(k2,1:3))));

                end
            end
        end
        k3 = k3+1;
    end
end
%% save data in a text file

fprintf(output,'octahedral_distortion = %f \n',octahedral_distortion);

fprintf(output,'theta_ab = %f \n',theta_ab(:));

fclose(output);

tic
%% end of code
