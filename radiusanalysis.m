close all
clear all
clc

%% Raw images

img_bf = mat2gray(double(imread('Particle images/bf_particle.jpg')));
img_fl = double(imread('Particle images/fluo_particle.jpg'));

figure;
set(gcf,'color','white');
subplot(2,2,1);
imshow(img_bf,[]);
subplot(2,2,2);
imshow(img_fl,[]);
[XX,YY] = meshgrid(1:size(img_bf,2),1:size(img_bf,1));
% set(gcf,'Position',get(0,'ScreenSize'));

% [centers,radii] = imfindcircles(img_bf,[300,700],'ObjectPolarity','bright');

% subplot(2,2,1); hold on;
% viscircles(centers,radii,'EdgeColor','r');

bw_bf = imbinarize(mat2gray(img_bf));
bw_fl = imbinarize(mat2gray(img_fl));
xc = 767; yc = 700; 
R = 600;
% [xc,yc] = ginput(1);

% 
figure;
subplot(2,2,1);
imshow(img_bf,[]); hold on; set(gcf,'color','white');
plot(xc,yc,'rx','linewidth',2);
subplot(2,2,3);
imshow(img_fl,[]); hold on;
plot(xc,yc,'mx','linewidth',2);

set(gcf,'Position',get(0,'ScreenSize'));

%% Visualisation of interpolation
% theta = linspace(0,2*pi,200);
% 
% % for i = 1:length(theta)
% for i = 48
%     xq = linspace(xc,xc+R.*cos(theta(i)),100);
%     yq = linspace(yc,yc+R.*sin(theta(i)),100);
%     r = sqrt((xq-xc).^2 + (yq-yc).^2);
%     vq_bf = interp2(XX,YY,img_bf,xq,yq);
%     vq_fl = interp2(XX,YY,img_fl,xq,yq);
%     subplot(2,2,1);
%     p1 = plot(xq,yq,'-r','linewidth',2);
%     subplot(2,2,3);
%     p2 = plot(xq,yq,'-m','linewidth',2);
%     subplot(2,2,[2,4]);
%     plot(r,vq_bf,'r','linewidth',1.5);
%     hold on;
%     plot(r,vq_fl,'m','linewidth',1.5);
%     pause(0.05);
%     cla;
%     delete(p1); delete(p2);
% 
% 
%     % Find max of bf radius
%     [~,ind] = max(vq_bf);
%     rdata(i) = r(ind);
% 
% end
% 
% figure;

%%
close all
theta = linspace(0,2*pi,200);
dx = 5;
clear xb yb rdata
for i = 1:length(theta)
% for i = 10
    100*i/length(theta)
% for i = 48
    xq = linspace(xc,xc+R.*cos(theta(i)),100);
    yq = linspace(yc,yc+R.*sin(theta(i)),100);
    r = sqrt((xq-xc).^2 + (yq-yc).^2);
    vq_bf = interp2(XX,YY,img_bf,xq,yq);
    vq_fl = interp2(XX,YY,img_fl,xq,yq);

    % Find max of bf radius
    [~,ind] = max(vq_bf);
    previous = mean(vq_bf(1:ind-1));
    next = mean(vq_bf(ind:end));
    ind2 = find(vq_bf(ind:end)<0.4,1,'first');
    ind2 = ind2 + ind-1;
    % Find first point where the radius drops
%     try
        rdata(i) = r(ind2);
        xb(i) = xq(ind2); 
        yb(i) = yq(ind2);
%     catch
%         xb(i) = NaN;
%         yb(i) = NaN;
%         rdata(i) = NaN;
%     end
    % Peak fit at this value

end
%%
close all
inds = 1:length(theta);
% rdata = rdatabackup;
d = abs(rdata - mean(rdata));
inds(inds(d > 3*std(rdata))) = [];
xb = xb(inds);
yb = yb(inds);
rdata = rdata(inds);

figure;
plot(r,vq_bf,'k','linewidth',2); hold on;
plot(r(ind2),vq_bf(ind2),'rx');
%
figure;
histogram(rdata,'normalization','pdf');

disp(['Mean radius: ' num2str(mean(rdata)) ', standard deviation: ' num2str(std(rdata))])
%%
figure;
imshow(img_bf,[]);
hold on
plot(xb,yb,'r-','linewidth',1.5);

x_filt = sgolayfilt(xb,2,15); y_filt = sgolayfilt(yb,2,15);
plot(x_filt,y_filt,'-y','linewidth',1.5);

%%
IN = inpolygon(XX,YY,xb,yb);

fluorescence_all = img_fl(IN);
fluorescence_mean = mean(fluorescence_all);

figure;
histogram(fluorescence_all,'normalization','pdf');

% 




% subplot(2,2,2);
% imshow(bw_bf)
% 
% subplot(2,2,3);
% BWout = ~bw_bf;
% BWout = imfill(bw_bf,'holes');
% % imshow(BWout);
% 
% statOuter = regionprops(BWout,{'EquivDiameter','Centroid'});
% outerRadius = statOuter.EquivDiameter/2;
% 
% % Measure the inner radius
% BWin = imclearborder(img_bf);
% BWin = imopen(BWin, strel('disk',5)); % Remove noise
% BWin = ~BWin;
% 
% statInner = regionprops(BWin,{'EquivDiameter','Centroid'});
% innerRadius = statInner.EquivDiameter/2;
% % Show the result
% subplot(2,2,4);
% imshow(bw_bf);
% hold on
% viscircles(statOuter.Centroid, outerRadius,'Color','r')