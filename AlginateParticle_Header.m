close all
clear all
clc

matlabrc;
addpath('peripherals');

%%

DataDir = 'D:/Elisa/';
filename = 'IBIDI_061124'; % Leave extension off here
OutputDir = ['Outputs/' filename '/']; mkdir(OutputDir);
znum = [6,7,7,6,7,7]; % Which z slice to work with
Positions = [1,3,5,6]; % Which positions you want
Rline = 600; % Starting interploation line length, larger is better
THRES = 0.4; % Outer radius cutoff threshold
ORDER = 2; WINDOW = 15; % sgolay parameters
dR = 40; % Length in pixels of thickness of ring for edge analysis
Ntheta = 400;

bfr = BioformatsImage([DataDir filename '.nd2']);
NumT = bfr.sizeT; % Number of time points
NumXY = bfr.seriesCount; % Number of XY points
disp(['Channel 1: ' bfr.channelNames{1,1} ', Channel 2: ' bfr.channelNames{1,2}]);

fluorChannel = 1; brightfieldChannel = 2; % Default, check per file!
pixel_size = bfr.pxSize;

%% Setup for saving

DataOutDir = [OutputDir '/ProcessedData/']; mkdir(DataOutDir);
mkdir([DataOutDir 'Data/']);
ImageFigDir = [DataOutDir 'Figs/']; mkdir(ImageFigDir);
ImagePNGDir = [DataOutDir 'PNGS/']; mkdir(ImagePNGDir);

params.znum = znum; params.filename = filename;
params.positions = Positions;
params.THRES = THRES;
params.sgolayparams = [ORDER,WINDOW];
params.dR = dR;
params.NumT = NumT;
params.NumXY = NumXY;
params.fluorChannel = fluorChannel; 
params.brightfieldChannel = brightfieldChannel;

%% Get initial positions

% If running for the first time, get rough positions of the middle
if ~exist([OutputDir 'ParticleOrigins.mat'])
    fh = figure; set(fh,'color','white'); 
    disp('Click once somewhere in the middle of the particle')
    for ixy = 1:NumXY
        img_bf = getPlane(bfr,znum,brightfieldChannel,1,ixy);
        imshow(img_bf,[]);
        hold on;
        [XC(ixy),YC(ixy)] = ginput(1);
        clf;
    end
    save([OutputDir 'ParticleOrigins.mat'],'XC','YC','filename');
else
    load([OutputDir 'ParticleOrigins.mat'],'XC','YC','filename');
end

%% Setup saving directories
for i = 1:NumXY
    mkdir([DataOutDir 'Position_' sprintf('%02d',i) '/'])
end

%% Main

WB = waitbar(0,'Processing....');
for i = 116:NumT
    
    
    tic
% for i = 100
    
    for ixy = Positions
        waitbar(i/NumT,WB,['Time point ' num2str(i) ' of ' num2str(NumT) ', Position ' num2str(ixy) ' of ' num2str(NumXY)]);
%     for ixy = 1
        % Get starting position of this particle
        xc = XC(ixy); yc = YC(ixy);
        output.filename = filename;

        % Get images and setup meshgrid
        img_bf = mat2gray(double(getPlane(bfr,znum(ixy),brightfieldChannel,i,ixy)));
        img_fl = double(getPlane(bfr,znum(ixy),fluorChannel,i,ixy));
        [XX,YY] = meshgrid(1:size(img_bf,2),1:size(img_bf,1));

        % Find outer radius
        theta = linspace(0,2*pi,Ntheta);
        for itheta = 1:length(theta)
            xq = linspace(xc,xc+Rline.*cos(theta(itheta)),200);
            yq = linspace(yc,yc+Rline.*sin(theta(itheta)),200);
            r = sqrt((xq-xc).^2 + (yq-yc).^2);
            vq_bf = interp2(XX,YY,img_bf,xq,yq);
            vq_fl = interp2(XX,YY,img_fl,xq,yq);
    
            % Find max of bf radius
            [~,ind] = max(vq_bf);
            previous = mean(vq_bf(1:ind-1));
            next = mean(vq_bf(ind:end));
            ind2 = find(vq_bf(ind:end)<THRES,1,'first');
            if size(ind2,2) == 0
                ind2 = ind;
            else
                ind2 = ind2 + ind-1;
            end
            
            % Get rough radius and xy boundary
            rdata(itheta) = r(ind2);
            xb(itheta) = xq(ind2); 
            yb(itheta) = yq(ind2);

        end

        % Remove any strong outliers
        inds = 1:length(theta);
        d = abs(rdata - median(rdata));
        inds(inds(d > 60)) = [];
        xb = xb(inds);
        yb = yb(inds);
        rdata = rdata(inds);
        theta_b = theta(inds);

        output.rdata = rdata;
        output.theta_b = theta_b;

        % Filtered boundary
        if length(xb) < WINDOW
            if ixy == Positions(end)
                sgtitle([filename ' frame:' num2str(i)],'Interpreter','None');
                save([DataOutDir 'Data/frame_' sprintf('%04d',i) '.mat'],...
                    'output','filename','params')
                saveas(gcf,[ImageFigDir 'frame_' sprintf('%04d',i) '.fig'])
                set(gcf,'Position',get(0,'ScreenSize'));
                saveas(gcf,[ImagePNGDir 'frame_' sprintf('%04d',i) '.png'])
                close(gcf);
                toc
            end
            continue
        end

        x_filt = sgolayfilt(xb,ORDER,WINDOW); 
        y_filt = sgolayfilt(yb,ORDER,WINDOW);
        output.boundary_filt = [x_filt;y_filt];

        % Best fit circle to filtered boundary
        [xc_f,yc_f,R_f,a] = circfit(x_filt,y_filt);
        output.origin_fit = [xc_f,yc_f];
        output.radius_fit = R_f;

        % Find radius to filtered boundary with fitted centre
        rdata_filt = sqrt((xc_f - x_filt).^2 + (yc_f - y_filt).^2);
        output.rdata_filt = rdata_filt;

        % Find pixels in polygon
        IN = inpolygon(XX,YY,x_filt,y_filt);
        fluor_all = img_fl(IN);
        output.fluor_mean = mean(fluor_all(:));

        % Make new boundary inside to define annulus
        xb_inner = xc + (rdata-dR).*cos(theta_b);
        yb_inner = yc + (rdata-dR).*sin(theta_b);
        xb_inner_filt = sgolayfilt(xb_inner,ORDER,WINDOW); 
        yb_inner_filt = sgolayfilt(yb_inner,ORDER,WINDOW);
        output.innerboundary_filt = [xb_inner_filt;yb_inner_filt];

        % Find pixels in inner area
        IN_INNER = inpolygon(XX,YY,xb_inner,yb_inner);
        annulus = img_fl(IN == 1 & IN_INNER ~= 1);
        output.annulus_mean = mean(annulus(:));     

        % Areas
        output.particle_area = area(polyshape(x_filt,y_filt));
        output.annulus_area = output.particle_area - area(polyshape(xb_inner_filt,yb_inner_filt));
        
        % Saving
        save([DataOutDir 'Position_' sprintf('%02d',i) '/frame_' sprintf('%04d',i) '.mat'],...
                'output','filename','params');

        %% Figure making
        figure(1);
        set(gcf,'color','white'); box on; hold on;
        subplot(2,3,ixy);
        imshow(img_bf,[]); hold on;
        pfilt = plot(x_filt,y_filt,'m-','linewidth',1.5);
        pb = plot(x_filt,y_filt,'r-','linewidth',1.5);
        pA = plot(xb_inner_filt,yb_inner_filt,'y-','linewidth',1.5);
        title(['Position ' num2str(ixy)]);

        if ixy == Positions(end)
            sgtitle([filename ' frame:' num2str(i)],'Interpreter','None');
%             save([DataOutDir 'Data/frame_' sprintf('%04d',i) '.mat'],...
%                 'output','filename','params')
            saveas(gcf,[ImageFigDir 'frame_' sprintf('%04d',i) '.fig'])
            set(gcf,'Position',get(0,'ScreenSize'));
            saveas(gcf,[ImagePNGDir 'frame_' sprintf('%04d',i) '.png'])
            close(gcf);
            toc
        end
        
        
    end % End of looping over time
end % End of looping over XY


