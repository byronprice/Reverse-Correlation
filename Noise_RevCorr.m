function [] = Noise_RevCorr(AnimalName,NoiseType)
%Noise_RevCorr.m
%   Display a series of colored noise stimuli to infer the receptive fields
%    See eg Smyth et al. 2003 Receptive Field Organization ...
%
% Briefly, if the receptive field for a given neuron is described by a
%  p-by-1 vector (p being the number of pixels), f, and the stimuli described
%  by a s-by-p (s being the number of stimuli presented) matrix, S, and the 
%  response of that neuron described by an s-by-1 vector, then we can write
%  the responses in terms of the stimuli and the receptive field as 
%  r = S * f .  We want to infer the spatial receptive field structure
%  f, so we need S^-1 * r = S^-1 * S * f => f = S^-1 * r ... So, the
%  following code will present a series of white noise stimuli with each
%  stimulus lasting flipInterval msec, followed by flipInterval msec of flat 
%  grey.  The stimuli will be output as the matrix S, along with the 
%  timeStamps for when they were generated.
% Images are displayed on the screen in Image Coordinates
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Optional Inputs
%       NoiseType - 'white' or 'pink' or 'brown' ... defaults to pink
%
%       see file NoiseVars.mat for more changeable presets
%
%OUTPUT: file named 'NoiseStimDate_AnimalName.mat' , e.g. 
%          NoiseStim20160718_12345.mat
%
%        S - matrix sized numStimuli-by-numPixels that represent each of
%          the stimuli presented, try 
%          image = reshape(S(1,:),[width,height]); to view one of the white
%          noise stimuli
%        effectivePixels - effective width (and height) of the display window 
%          in pixels, chosen to fill most of the screen, with the exception
%          of a small strip on each horizontal edge of the screen
%        DistToScreen - 25 cm for now
%
% Created: 2016/03/04, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/07/25
% By: Byron Price

cd('~/CloudStation/ByronExp/NoiseRetino')
load('NoiseVars.mat');

switch nargin
    case 1
        NoiseType = 'pink';
end

% for motion-contingent display / interaction with recording computer
startEXP = 254;
endEXP = 255;

tcpipClient = tcpip('128.197.59.169',30000,'NetworkRole','client');
bufferSize = 50000; % bytes, (we won't need this much)
set(tcpipClient,'InputBufferSize',bufferSize);
set(tcpipClient,'Timeout',5);
fopen(tcpipClient);

% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;
usb = usb1208FSPlusClass;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

TimeEstimate = numStimuli*(flipInterval+WaitTime+0.5)/60+4*0.5;
fprintf('\nEstimated time is %3.2f minutes.',TimeEstimate);
WaitSecs(1);

% Choose screen with maximum id - the secondary display:
screenid = max(Screen('Screens'));

% % Open a fullscreen onscreen window on that display, choose a background
% % color of 127 = gray with 50% max intensity; 0 = black;255 = white
background = 127;
[win,~] = Screen('OpenWindow', screenid,background);

gammaTable = makeGrayscaleGammaTable(gama,0,255);
Screen('LoadNormalizedGammaTable',win,gammaTable);

% Query window size in pixels
[w_pixels,h_pixels] = Screen('WindowSize', win);
minPix = min(w_pixels,h_pixels);
maxPix = max(w_pixels,h_pixels)-200;

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
maxPix = maxPix-mod(maxPix,screenPix_to_effPix);
minPix = minPix-mod(minPix,screenPix_to_effPix);

effectivePixels = [maxPix/screenPix_to_effPix,minPix/screenPix_to_effPix];

gridSpacing = screenPix_to_effPix*conv_factor/10;
spatialSampleFreq = 1/gridSpacing; % in units of 1/cm

% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
%     spaceExp = 2;
elseif strcmp(NoiseType,'pink') == 1
    beta = -2;
%     spaceExp = 2;
elseif strcmp(NoiseType,'brown') == 1
    beta = -3;
%     spaceExp = 2;
else 
    display('NoiseType must be ''white'', ''pink'' or ''brown'' as a string.')
    return;
end

Grey = 127;
DIM = [effectivePixels(1),effectivePixels(2)];
S = zeros(numStimuli,DIM(1)*DIM(2),'uint8');
unbiasedS = zeros(numStimuli,DIM(1)*DIM(2),'uint8');
% below pink noise from Jon Yearsley, 1/f noise generate spatial data

for ii=1:numStimuli
    [X,Y] = spatialPattern([DIM(2),DIM(1)],beta);
    X = X-min(min(X));
    X = (X./max(max(X))).*255;
    Z = X(:);
    meanVal = mean(Z);difference = meanVal-Grey;
    %figure();imagesc(reshape(Y-difference,DIM));
    S(ii,:) = Z-difference;
    
    Y = Y-min(min(Y));
    Y = (Y./max(max(Y))).*255;
    Z = Y(:);
    meanVal = mean(Z);difference = meanVal-Grey;
    %figure();imagesc(reshape(Y-difference,DIM));
    unbiasedS(ii,:) = Z-difference;
end
S = uint8(S);
unbiasedS = uint8(unbiasedS);
% Sdisplay = S;
% Sdisplay(S<60) = 0;Sdisplay(S>=60 & S<196) = 127;Sdisplay(S>=196) = 255;

shuffled = randperm(numStimuli);
test = shuffled(1:350);
train = shuffled(350+1:end);

% load('NaturalImageSet.mat');
% numNaturalImages = size(NaturalImSet,1);
% for ii=1:length(test)
%    index = random('Discrete Uniform',numNaturalImages);
%    temp = NaturalImSet{index};
%    S(test(ii),:) = temp(:);
% end

wLow = round((w_pixels-maxPix)/2);
wHigh = round(w_pixels-wLow);
hLow = round((h_pixels-minPix)/2);
hHigh = round(h_pixels-hLow);
destRect = [wLow hLow wHigh hHigh];

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);

WaitTimes = WaitTime+unifrnd(0,1,[numStimuli,1]);

usb.startRecording;WaitSecs(1);usb.strobeEventWord(0);
WaitSecs(30);
usb.strobeEventWord(startEXP);WaitSecs(1);
tt = 1;
vbl = Screen('Flip',win);
while tt <= numStimuli
    if tcpipClient.BytesAvailable > 0
        data = fread(tcpipClient,tcpipClient.BytesAvailable/8,'double');
        if sum(data) > 0
            WaitSecs(2);
        else
            % Convert it to a texture 'tex':
            Img = reshape(S(tt,:),[DIM(2),DIM(1)]);
            tex = Screen('MakeTexture',win,Img);
            Screen('DrawTexture',win, tex,[],destRect,[],0); % 0 is nearest neighbor
            % 1 is bilinear filter
            vbl = Screen('Flip',win);usb.strobeEventWord(1);
            vbl = Screen('Flip',win,vbl-ifi/2+flipInterval);usb.strobeEventWord(2);
            vbl = Screen('Flip',win,vbl-ifi/2+WaitTimes(tt));
            Screen('Close',tex);
            tt = tt+1;
        end
    end
    if mod(tt,round(numStimuli/3)) == 0
        Screen('Flip',win);usb.strobeEventWord(0);
        WaitSecs(30);
    end
end

Screen('Flip',win);usb.strobeEventWord(0);
usb.strobeEventWord(endEXP);
WaitSecs(1);
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);
% NoiseType = 'pinkHC';
Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseStim%s%d_%d.mat',NoiseType,Date,AnimalName);
save(filename,'S','numStimuli','flipInterval','effectivePixels',...
    'DistToScreen','screenPix_to_effPix','minPix','NoiseType',...
    'conv_factor','WaitTimes','beta','DIM','spatialSampleFreq',...
    'maxPix','test','train','unbiasedS');
end

function gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
% Generates a 256x3 gamma lookup table suitable for use with the
% psychtoolbox Screen('LoadNormalizedGammaTable',win,gammaTable) command
% 
% gammaTable = makeGrayscaleGammaTable(gamma,blackSetPoint,whiteSetPoint)
%
%   gamma defines the level of gamma correction (1.8 or 2.2 common)
%   blackSetPoint should be the highest value that results in a non-unique
%   luminance value on the monitor being used (sometimes values 0,1,2, all
%   produce the same black pixel value; set to zero if this is not a
%   concern)
%   whiteSetPoint should be the lowest value that returns a non-unique
%   luminance value (deal with any saturation at the high end)
% 
%   Both black and white set points should be defined on a 0:255 scale

gamma = max([gamma 1e-4]); % handle zero gamma case
gammaVals = linspace(blackSetPoint/255,whiteSetPoint/255,256).^(1./gamma);
gammaTable = repmat(gammaVals(:),1,3);
end

