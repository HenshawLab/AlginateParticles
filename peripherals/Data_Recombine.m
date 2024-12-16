

radius_average = NaN(NumT,NumXY);
radius_stderr = radius_average;
radius_fit = radius_average;
area_fullparticle = radius_average;
area_annulus = radius_average;

for i = 1:NumT
    for ix = 1:NumXY
        try
        load([DataOutDir 'Position_' sprintf('%02d',i) '/frame_' sprintf('%04d',i) '.mat']);
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

save([OutputDir '/ProcessedData/CollatedResults.mat'],...
    'radius_average','radius_stderr','radius_fit','fluor_mean',...
    'fluor_ring','area_fullparticle');