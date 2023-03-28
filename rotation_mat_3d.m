function m2 = rotation_mat(m1,th,fi,si)

th = th*(pi/180); % along x axis
fi = fi*(pi/180); % alpng z axis
si = si*(pi/180); %  along y axis

 
rmat = [cos(fi)*cos(si)-sin(fi)*sin(si)*sin(th), sin(fi)*cos(th), cos(fi)*sin(si)+sin(fi)*cos(si)*sin(th);
    -sin(fi)*cos(si)-cos(fi)*sin(si)*sin(th), cos(fi)*cos(th), -sin(fi)*sin(si)+cos(fi)*cos(si)*sin(th);
    -sin(si)*cos(th),-sin(th),cos(si)*cos(th)];
   
   m2 = rmat*m1';
