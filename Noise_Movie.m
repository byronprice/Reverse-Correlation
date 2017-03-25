function [] = Noise_Movie(AnimalName,NoiseType)
%Noise_Movie.m
%   Display a movie of noise stimuli to infer the receptive fields
%    of neurons using reverse correlation or a GLM point process model
%    See Smyth et al. 2003 Receptive Field Organization ...
%
%
%INPUT: AnimalName - unique identifier for the animal as a number, e.g.
%            12345
%       Optional Inputs
%       NoiseType - 'white' or 'pink' or 'brown' ... defaults to pink
%
%       see file MovieNoiseVars.mat for more changeable presets
%
%OUTPUT: file named 'NoiseStimDate_AnimalName.mat' , e.g. 
%          NoiseStim20160718_12345.mat
%
%        S - matrix sized numStimuli-by-numPixels that represent each of
%          the stimuli presented, try 
%          image = reshape(S(1,:),[width,height]); to view one of the white
%          noise stimuli
%        effectivePixels - effective width (and height) of the display window 
%          in pixels, first chosen as the minimum of the width and height 
%          of the current display ... the display is forced to be square
%          and the width in effectivePixels is 1/4 the width in true screen
%          pixels
%        DistToScreen - 25 cm for now
%
% Created: 2017/03/24, 24 Cummington, Boston
%  Byron Price
% Updated: 2017/03/24
% By: Byron Price

cd('~/CloudStation/ByronExp/NoiseRetino')
load('MovieNoiseVars.mat');

switch nargin
    case 1
        NoiseType = 'pink';
end
% Acquire a handle to OpenGL, so we can use OpenGL commands in our code:
global GL;
usb = usb1208FSPlusClass;

% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

numStimuli = 1000; %movieTime_Seconds*movie_FrameRate;

fprintf('\nEstimated time is %3.2f minutes.',movieTime_Seconds/60);
WaitSecs(10);

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
maxPix = max(w_pixels,h_pixels)-250;

% screen size in millimeters and a conversion factor to get from mm to pixels
[w_mm,h_mm] = Screen('DisplaySize',screenid);
conv_factor = (w_mm/w_pixels+h_mm/h_pixels)/2;

% perform unit conversions
degPerPix = atand((1*conv_factor)/(DistToScreen*10));

horzDegrees = atand((maxPix*conv_factor)/(DistToScreen*10));
vertDegrees = atand((minPix*conv_factor)/(DistToScreen*10));

% a bit confusing, we want the stimuli produced to have a certain number of
%  effective pixels, which project to larger squares of on-screen pixels
screenPix_to_effPix = 30;
maxPix = maxPix-mod(maxPix,screenPix_to_effPix);
minPix = minPix-mod(minPix,screenPix_to_effPix);
numPixels = maxPix*minPix;

effectivePixels = [maxPix/screenPix_to_effPix,minPix/screenPix_to_effPix];

% GENERATION OF NOISE
if strcmp(NoiseType,'white') == 1
    beta = 0;
elseif strcmp(NoiseType,'pink') == 1
    beta = -3;
elseif strcmp(NoiseType,'brown') == 1
    beta = -6;
else 
    fprintf('NoiseType must be ''white'', ''pink'' or ''brown'' as a string.\n');
    return;
end

Grey = 127;
%S = zeros(effectivePixels(1),effectivePixels(2),numStimuli,'uint8');

% below white/pink/brown noise from Jon Yearsley
DIM = [effectivePixels(1),effectivePixels(2),numStimuli];

u = linspace(0,horzDegrees,DIM(1));v = linspace(0,vertDegrees,DIM(2));
t = linspace(0,1/movie_FrameRate,numStimuli);
[U,V,T] = meshgrid(u,v,t);
S_f = single((U.^2+V.^2+T.^2).^(beta/2));
clear U V T;
S_f(S_f==inf) = 1;
% S_f = S_f.^0.5;
noise = randn([DIM(2),DIM(1),DIM(3)],'single');
Y = fftn(noise);
Y = Y.*S_f;
X = ifftn(Y);
%X = X.*conj(X);
X = real(X);
X = X-min(min(min(X)));
X = (X./max(max(max(X)))).*255;
meanVal = mean(mean(mean(X)));difference = meanVal-Grey;
S = X-difference;
S = uint8(S);

Priority(9);
% Retrieve monitor refresh duration
ifi = Screen('GetFlipInterval', win);
% flipIntervals = flipInterval+exprnd(0.1,[numStimuli,1]);
% WaitTimes = WaitTime+exprnd(0.25,[numStimuli,1]);
flipInterval = 1/movie_FrameRate;

usb.startRecording;
WaitSecs(5);
tt = 1;
vbl = Screen('Flip', win);
while tt <= numStimuli
    % Convert it to a texture 'tex':
    Img = uint8(kron(single(S(:,:,tt)),ones(screenPix_to_effPix)));
    tex = Screen('MakeTexture',win,Img);
    Screen('DrawTexture',win, tex);
    vbl = Screen('Flip',win,vbl-ifi/2+flipInterval);usb.strobe;
    tt = tt+1;
    Screen('Close',tex);
end
WaitSecs(2);
usb.stopRecording;
% Close window
Screen('CloseAll');
Priority(0);

Date = datetime('today','Format','yyyy-MM-dd');
Date = char(Date); Date = strrep(Date,'-','');Date = str2double(Date);
filename = sprintf('NoiseMovieStim%s%d_%d.mat',NoiseType,Date,AnimalName);
save(filename,'S','numStimuli','movie_FrameRate','effectivePixels','movieTime_Seconds',...
    'DistToScreen','screenPix_to_effPix','minPix','NoiseType','degPerPix');
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