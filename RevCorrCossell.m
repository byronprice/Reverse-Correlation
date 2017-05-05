% RevCorr_Cossell.m

fileName = {'13_05_22','13_05_20','13_05_19','13_05_14'};
nExps = length(fileName);
[csCaData, sParameters, tfImageStack] = ...
   Load_Cossell_etal_2015_NaturalImageResponses(fileName);

convertedStack = zeros(size(tfImageStack,3),size(tfImageStack,1)*size(tfImageStack,2));
for ii=1:size(tfImageStack,3)
   tempIm = tfImageStack(:,:,ii);
   convertedStack(ii,:) = tempIm(:)';
end

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
        
        numNeurons = size(spikeData,1);
        numIms = length(evokedSlots);
        reducedSpikeData = sum(spikeData(:,:,3:5),3)./3;
        
        DIM = [size(tfImageStack,1),size(tfImageStack,2)];
        L = zeros(DIM(1)*DIM(2),DIM(1)*DIM(2));
        bigCount = 1;
        for kk=1:DIM(2)
            for ll=1:DIM(1)
                tempMat = zeros(DIM(1),DIM(2));
                tempMat(ll,kk) = 4;
                if kk > 1 && kk < DIM(2)
                    tempMat(kk-1,ll) = -1;
                    tempMat(kk+1,ll) = -1;
                end
                if ll > 1 && ll < DIM(1)
                    tempMat(kk-1,ll) = -1;
                    tempMat(kk+1,ll) = -1;
                end
                L(bigCount,:) = tempMat(:)';bigCount = bigCount+1;
            end
        end
        
        F = zeros(numNeurons,DIM(1)*DIM(2));
        bigLambda = logspace(0,4,10);
        train = round(numIms*0.7);
        for kk=1:numNeurons
            r = reducedSpikeData(kk,61:60+train)';
            constraints = [r;zeros(DIM(1)*DIM(2),1)];
            tempA = convertedStack(61:60+train,:);
            
            tempF = zeros(length(bigLambda),DIM(1)*DIM(2));
            RMS = zeros(length(bigLambda),1);
            for ll=1:length(bigLambda)
                A = [tempA;bigLambda(ll).*L];
                tempF(ll,:) = A\constraints;clear A;
                rTest = reducedSpikeData(kk,61+train:numIms)';
                A = convertedStack(61+train:numIms,:);
                RMS(ll) = norm(rTest-A*tempF(ll,:)')./sqrt(DIM(1)*DIM(2));
            end
             clear tempA A;
             [~,bestMap] = min(RMS);
             F(kk,:) = tempF(bestMap,:);
             figure();imagesc(reshape(F(kk,:),[DIM(1),DIM(2)]));
             title(sprintf('Im Plane %d - Neuron %d',jj,kk));
        end
        newFileName = strcat(fileName{ii},'-',num2str(jj));
        save(newFileName,'F','bigLambda');
    end
end