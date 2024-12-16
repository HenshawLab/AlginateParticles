close all
clear all
clc

%%

% indir = 'D:/Elisa/Outputs/IBIDI_061124/ProcessedData/Data/';
outmain = 'D:/Elisa/Outputs/IBIDI_061124/ProcessedData/';

% for i = 1:6
%     mkdir([outmain 'Position_' sprintf('%02d',i) '/']);
% end
% 
% for i = 1:151
%     i
% % for i = 1
%     try
%         load([indir 'frame_' sprintf('%04d',i) '.mat']);
%     catch
%         continue
%     end
%     OUT = output;
%     for ix = 1:6
%         output = OUT(ix);
%         output.position = ix;
%         save([outmain 'Position_' sprintf('%02d',ix) '/frame_' sprintf('%04d',i) '.mat']);
%         clear output
%     end
% end

datadir = 'D:/Elisa/Outputs/IBIDI_061124/ProcessedData/';

radius_average = NaN(151,6);
radius_stderr = radius_average;
radius_fit = radius_average;
area_fullparticle = radius_average;
area_annulus = radius_average;

for i = 1:151
    for ix = 1:6
        try
        load([outmain 'Position_' sprintf('%02d',ix) '/frame_' sprintf('%04d',i) '.mat']);
        radius_average(i,ix) = mean(output.rdata_filt);
        radius_stderr(i,ix) = std(mean(output.rdata_filt))./sqrt(numel(output.rdata_filt));
        radius_fit(i,ix) = output.radius_fit;
        fluor_mean(i,ix) = output.fluor_mean;
        fluor_ring(i,ix) = output.annulus_mean;
        area_fullparticle(i,ix) = output.particle_area;
        catch
            continue
        end
%         area
    end
end