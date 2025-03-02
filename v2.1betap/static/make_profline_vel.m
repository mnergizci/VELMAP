function[profile]=make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof,outdir)
%========================================================
% function[profile]=make_profline_vel(trim,fitmodel,vcmmodel,invenu,prof,outdir)
%
%  forward calculation for velocities on a profile
%
% INPUT:
%  trim:     triangular mesh
%  fitmodel: fitted velocity field
%  vcmmodel: vcm of fitted velocity field
%  invenu:   inversion parameters
%  prof:     structure of profile (x0,y0,x1,y1,swath,step,id)
%  outdir:   output directory (optional)
% 
% OUTPUT:
%   profile: all datapoint along/across the profile                            
%
% NOTE:
%  vcm format: [vcm_x1x1 ... vcm_x1xn vcm_x1y1 ... vcm_x1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_xnx1 ... vcm_xnxn vcm_xny1 ... vcm_xnyn]
%              [vcm_y1x1 ... vcm_y1xn vcm_y1y1 ... vcm_y1yn]
%              [ ...     ...   ...      ...    ...   ...   ]
%              [vcm_ynx1 ... vcm_ynxn vcm_yny1 ... vcm_ynyn]
%
%  profile(:,1): distance along the profile
%  profile(:,2): distance across the profile
%  profile(:,3): velx (along the profile)
%  profile(:,4): vely (across the profile)
%  profile(:,5): velu (vertical for 3d vel only)
%  profile(:,6): std_velx
%  profile(:,7): std_vely
%  profile(:,8): std_velu (for 3d vel only)
%
% Hua Wang, 10/01/2010
%
% 13/03/2016 HW: using km as unit of swath/step
% 27/03/2010 HW: support 3D velocity
%========================================================

%get the transform matrix
%[ cos(alpha)  sin(alpha) ]
%[-sin(alpha)  cos(alpha) ]
sx = (prof.x1-prof.x0);  %x coordinate
sy = (prof.y1-prof.y0);  %y coordinate
l = sqrt(sx*sx+sy*sy);
coef = [sx/l, sy/l; -sy/l, sx/l];

%coordinate of the profile points
if prof.step==0
  error('Step must be positive!');
end
%using km as the unit of step, 13/03/2016, HW
npt=ceil(l/prof.step*100);
for i=1:npt
  pt(i).lon=prof.x0+(i-1)*sx/npt;
  pt(i).lat=prof.y0+(i-1)*sy/npt;
end

%remove points outside of the mesh using tidygps.m
%[pt]=tidygps(trim,pt); %JRW
npt=length(pt);

%coordinate transform, unit degree
%xy=[pt.lon'-prof.x0,pt.lat'-prof.y0];
%coordinate transform, unit km
xy=ll2utm([[pt.lon]',[pt.lat]'],[prof.x0,prof.y0]);
pxy = coef*xy';
pxy = pxy';

%forward calculation for the profile point
nvtx=length(trim.x);
invs=sum(invenu);
gpsmat=designgps(trim,pt,invenu);
%gpsmat=gpsmat(1:invs*npt,1:invs*nvtx);
velfit=gpsmat*fitmodel(1:invs*nvtx);
vcmfit=gpsmat*vcmmodel(1:invs*nvtx,1:invs*nvtx)*gpsmat';

%project the velocities onto the profile
vel=reshape(velfit,npt,invs);
pvel = coef*vel(:,1:2)';
pvel = pvel';

%calculate error bar of the projected velocities
for i=1:npt
  tvcm=[vcmfit(i,i),vcmfit(i,i+npt);vcmfit(i,i+npt),vcmfit(i+npt,i+npt)];
  vcmx=coef*tvcm*coef';
  stdx(i,1:2)=sqrt(diag(vcmx));
end

%U-components
if invs==3
  pvel=[pvel vel(:,3)];
  stdx=[stdx sqrt(diag(vcmfit(2*npt+1:3*npt,2*npt+1:3*npt)))];
end

%save the profile
profile=full([pxy,pvel,stdx]);

outfile=strcat('prof',num2str(prof.id));
if nargin>5
  outfile=strcat(outdir,outfile);
end
save(char(outfile),'profile','-ASCII');
