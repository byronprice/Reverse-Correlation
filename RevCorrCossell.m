% RevCorr_Cossell.m

fileName = {'13_05_22'};
nExps = length(fileName);
[csCaData, sParameters, tfImageStack] = ...
   Load_Cossell_etal_2015_NaturalImageResponses(fileName);

for ii=1:nExps
    CaData = csCaData{ii};
    nImgPlanes = size(CaData,2);
    for jj=1:nImgPlanes
        spikeData = CaData(jj).spikes;
        nROIs = size(spikeData,1);
        nFrames = size(spikeData,3);
        nImageSlot = size(spikeData,2); % 1:60 and 1861:1920 are spontaneous activity
                        % 61:1860 are evoked activity
        spontaneousSlots = sParameters.vnSpontaneousFrames;
        evokedSlots = setdiff(1:max(spontaneousSlots),spontaneousSlots);
        imStackIndex = sParameters.vnPresentedImage;
    end
end