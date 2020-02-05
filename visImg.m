function [imgBrightness,g] = visImg(imgComplex,roomSize,voxelSize)

xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
% Combination of all of these to get coordinates for all the voxels
xyzVoxelCoord = combvec(xVoxel,yVoxel,zVoxel)';

nx = length(xVoxel);
ny = length(yVoxel);
nz = length(zVoxel);

imgBrightness = (abs(imgComplex).^2);
maxBrightness = max(imgBrightness(:)); 
minBrightness = min(imgBrightness(:));
imgBrightness = (imgBrightness-minBrightness)/(maxBrightness-minBrightness);

% Visualizing reconstructed image
x = reshape(xyzVoxelCoord(:,1),nx,ny,nz);
y = reshape(xyzVoxelCoord(:,2),nx,ny,nz);
z = reshape(xyzVoxelCoord(:,3),nx,ny,nz);

xslice = roomSize(1,1):voxelSize(1):roomSize(1,2);
yslice = roomSize(2,1):voxelSize(2):roomSize(2,2);
zslice = roomSize(3,1):voxelSize(3):roomSize(3,2);

p = [2 1 3];
x = permute(x, p);
y = permute(y, p);
z = permute(z, p);
imgBrightness = reshape(imgBrightness,[nx ny nz]);
imgBrightness = permute(imgBrightness, p);

%gradient detection instead of threshold
imgthresh=0.6;
imgBrightness(imgBrightness<imgthresh) = 0;
[gx1,gy1,gz1] = gradient(imgBrightness);
g1=sqrt(gx1.^2 +gy1.^2 +gz1.^2);
[gx2,gy2,gz2] = gradient(g1);
g2=sqrt(gx2.^2 +gy2.^2 +gz2.^2);
%g=g1;

  imgthresh=0.8;
 imgBrightness(imgBrightness<imgthresh) = 0;
numx=length(xVoxel);
numy=length(yVoxel);
numz=length(zVoxel);
g3=zeros(numx,numy,numz);
g3=g1./imgBrightness;
for nx=1:numx
    for ny=1:numy
        for nz=1:numz
           
            if imgBrightness(nx,ny,nz)==0 
                g3(nx,ny,nz)=0;
           
            end
        end
    end
end
g=g1;
% imgBrightnessf = csaps({xVoxel,yVoxel,zVoxel},imgBrightness,0.9999);
% imgSpline = fnval(imgBrightnessf,{xVoxel,yVoxel,zVoxel}) ;
% figure
% fnplt(imgBrightnessf)
% [gx,gy,gz] = gradient(imgSpline);
% g=sqrt(gx.^2 +gy.^2 +gz.^2);
% gthresh=0.15;
% g(g<gthresh) = 0;
% g=g./imgBrightness;
% for nx=1:numx
%     for ny=1:numy
%         for nz=1:numz
%             thres1=0.3;   thres2=0.05; 
%             if imgSpline(nx,ny,nz)< thres1 &&  g(nx,ny,nz)> thres2
%                 g(nx,ny,nz)=0;
%                 imgSpline(nx,ny,nz)=0;
%             end 
%             end
%         end
%    end
% gthresh=0.3;
% g(g<gthresh) = 0;
% % figure
% % fnplt(imgBrightnessf)
% maxg=max(g(:));

% [gx,gy,gz] = gradient(imgBrightness);
% g=sqrt(gx.^2 +gy.^2 +gz.^2);
% gthresh=0.1;
% g(g<gthresh) = 0;
% g=g./imgBrightness;
% numx=length(xVoxel);
% numy=length(yVoxel);
% numz=length(zVoxel);
% for nx=1:numx
%     for ny=1:numy
%         for nz=1:numz
%             thres1=0.3;   thres2=0.05; 
%             if imgBrightness(nx,ny,nz)< thres1 &&  g(nx,ny,nz)> thres2
%                 g(nx,ny,nz)=0;
%                 imgBrightness(nx,ny,nz)=0;
%             end 
%             end
%         end
%    end

% gthresh=0.1;
% g(g<gthresh) = 0;
% maxg=max(g(:));







figure
%h = slice(x,y,z,imgBrightness,[],[],[0.1]);
h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
xlim([roomSize(1,1),  roomSize(1,2)])
ylim([roomSize(2,1),  roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')
a = alphamap('rampup',256);
imgThresh = 100;
a(1:imgThresh)=0;
alphamap(a); 
%title('11a6 R20    Zmin-50 Zmax 150   Xcenter600 Ycenter600 ')
hold on
figure
h = slice(x,y,z,g1,xslice,yslice,zslice);
%h = slice(x,y,z,g1,[],[],[0.1]);
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
xlim([roomSize(1,1),  roomSize(1,2)])
ylim([roomSize(2,1),  roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')
a = alphamap('rampup',256);
imgThresh = 180;
a(1:imgThresh)=0;
alphamap(a); 

% figure
% h = slice(x,y,z,g2,xslice,yslice,zslice);
% %h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
% xlabel('x (m)','FontSize',14)
% ylabel('y (m)','FontSize',14)
% zlabel('z (m)','FontSize',14)
% xlim([roomSize(1,1),  roomSize(1,2)])
% ylim([roomSize(2,1),  roomSize(2,2)])
% zlim([roomSize(3,1), roomSize(3,2)])
% 
% set(h, 'EdgeColor','none',...
%     'FaceColor','interp',...
%     'FaceAlpha','interp');
% alpha('color')
% a = alphamap('rampup',256);
% imgThresh = 150;
% a(1:imgThresh)=0;
% alphamap(a); 


% figure
% h = slice(x,y,z,g3,xslice,yslice,zslice);
% %h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
% xlabel('x (m)','FontSize',14)
% ylabel('y (m)','FontSize',14)
% zlabel('z (m)','FontSize',14)
% xlim([roomSize(1,1),  roomSize(1,2)])
% ylim([roomSize(2,1), roomSize(2,2)])
% zlim([roomSize(3,1), roomSize(3,2)])
% 
% set(h, 'EdgeColor','none',...
%     'FaceColor','interp',...
%     'FaceAlpha','interp');
% alpha('color')
% a = alphamap('rampup',256);
% imgThresh = 180;
% a(1:imgThresh)=0;
% alphamap(a); 
end