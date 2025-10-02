function msExtractSFPs(ms)
%msExtractSFPs Extracts spatial footprints to perform chronic re-alignment
% Converts spatial footprints from m,k,n (UCLA) to n,m,k (Ziv's lab) where
% n is the number of neurons, k pixels in x axis, and m pixels in y axis
%
% Author: Guillaume Etter
% Contact: etterguillaume@gmail.com

for cell_i = 1:size(ms.SFPs,3);
    SFP_temp = ms.SFPs(:,:,cell_i);
    SFP_temp(SFP_temp<0.5*max(max(SFP_temp))) = 0; % This is to sharpen footprints, based on Ziv lab method
    SFP(cell_i,:,:) = SFP_temp;
end

figure;
img = max(permute(SFP,[2 3 1]),[],3);
quantile(img(img>0), 0:0.1:1, 'all')
imagesc(img, [0, quantile(img(img>0), 0.9)]);

fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPositionMode = 'auto'; fig.PaperPosition = [0 0 8 8]; 
print([ms.dirName '\msCam\dFF_f1-18000_source_extraction\ROI.png'], '-dpng', '-r300');

save([ms.dirName '/SFP_cut2.mat'],'SFP','-v7.3');

end

