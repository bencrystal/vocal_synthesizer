
%% Record Sample
close all; clc; clear;
tSeconds = 0;


while (tSeconds <= 0) || (tSeconds > 8)
    tSeconds = input('How many seconds (1-8) would you like your recorded sample to be? ')
end





% disp(tSeconds + ' has been selected.');

Fs = 16000; %Sets Sampling Frequency
nBits = 16; %8, 16, or 24 bits per sample 
nChannels = 1; %Mono audio
seconds = tSeconds; %Number of seconds recorded

disp('Get Ready!');
disp('3');
pause(1);
disp('2');
pause(1);
disp('1');
pause(1);

tone = audiorecorder(Fs, nBits, nChannels);
disp('Speak Now');
recordblocking(tone,seconds)
disp('Recording Done');
voice = getaudiodata(tone); %Contains 24000 data points for the audio file

soundsc(voice(1:length(voice)),Fs);

    
%% Pitch Modulation 

%(pvoc)
n = 1024; %FFT Size
hop = n/4; %25 percent window overlap for smooth reconstruction

f = n;
w = n;
x = voice;
%(stft)



%X = stft(voice, n, n, hop);

% expect x as a row
if size(x,1) > 1
  x = x';
end

s = length(x);

if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,f);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = f/2;   % midpoint of win
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen)); %hanning equation
    win = zeros(1, f);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
  end
  else
  win = w;
end

w = length(win);
% now can set default hop
if hop == 0
  hop = floor(w/2);
end

c = 1;

% pre-allocate output array
d = zeros((1+f/2),1+fix((s-f)/hop));

for b = 0:hop:(s-f)
  u = win.*x((b+1):(b+f));
  t = fft(u); %window by window FFTs
  d(:,c) = t(1:(1+f/2))'; %fills in output array
  c = c+1;
end;

  % Otherwise, no plot, but return STFT
  D = 1.0*d;



%% (pvoc)

% Calculate the new timebase samples
[rows, cols] = size(D);

%% SCALING FACTOR
r2 = 8/9; %major 2
rd5 = 7071/10000;%1/sqrt(2); %diminished 5
ra5 = 5/8;%augmented 5
rd7 = 3/5; %dim 7
r = 4/5; % MAJOR 3 -- timescale signal r times faster, CHANGEABLE (.8 for major thirds, .6667 for fifths, .5 for octaves
rm = 5/6;%261.6256/311.1270; % minor 3
r4 = 3/4;
r5 = 2/3;
rM7 = 261.6256/493.8833;
rm7 = 261.6256/466.1638;
r8 = 1/2;
%%
td7 = 0:rd7:(cols-2);
t2 = 0:r2:(cols-2);
td5 = 0:rd5:(cols-2);
ta5 = 0:ra5:(cols-2);
t = 0:r:(cols-2);
tm = 0:rm:(cols-2);
t4 = 0:r4:(cols-2);
t5 = 0:r5:(cols-2);
tM7 = 0:rM7:(cols-2);
tm7 = 0:rm7:(cols-2);
t8 = 0:r8:(cols-2);
% Have to stay two cols off end because (a) counting from zero, and 
% (b) need col n AND col n+1 to interpolate


%% NEW PVSAMPLE SECTION





% Generate the new spectrogram
%DEFINE HOP?
[rows,cols] = size(d); %size of stft array (b) MAY NEED TO RENAME D? rows2?
N = 2*(rows-1); %fft size

if hop == 0
  % default value
  hop = N/2;
end

% Expected phase advance in each bin
dphi = zeros(1,N/2+1);
%%
dphi(2:(1 + N/2)) = (2*pi*hop)./(N./(1:(N/2)));

% Phase accumulator
% Preset to phase of first frame for perfect reconstruction
% in case of 1:1 time scaling
ph = angle(d(:,1));

% Append a 'safety' column on to the end of b to avoid problems 
% taking *exactly* the last frame (i.e. 1*b(:,cols)+0*b(:,cols+1))
d = [d,zeros(rows,1)];

ocol = 1;

%% PVSAMPLE Time Shifting
cd7 = zeros(rows, length(td7));
c2 = zeros(rows, length(t2));
cd5 = zeros(rows, length(td5));
ca5 = zeros(rows, length(ta5));
c = zeros(rows, length(t));
cm = zeros(rows, length(tm));
c4 = zeros(rows, length(t4));
c5 = zeros(rows, length(t5));
cM7 = zeros(rows, length(tM7));
cm7 = zeros(rows, length(tm7));
c8 = zeros(rows, length(t8));


for tt = td7
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  cd7(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = t2
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  c2(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = td5
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  cd5(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = t
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  c(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = ta5
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  ca5(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = t4
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  c4(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = t5
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  c5(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = tm
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  cm(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = tM7
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  cM7(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = tm7
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  cm7(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end

for tt = t8
  % Grab the two columns of b
  
    dcols = d(:,floor(tt)+[1 2]);
    %dcols = d(:,[1 2]); %Fixes indexing error?
  
  tf = tt - floor(tt);
  dmag = (1-tf)*abs(dcols(:,1)) + tf*(abs(dcols(:,2)));
  % calculate phase advance
  dp = angle(dcols(:,2)) - angle(dcols(:,1)) - dphi';
  % Reduce to -pi:pi range
  dp = dp - 2 * pi * round(dp/(2*pi));
  % Save the column
  c8(:,ocol) = dmag .* exp(j*ph);
  % Cumulate phase, ready for next frame
  ph = ph + dphi' + dp;
  ocol = ocol+1;
end



% c is the interpolated stft array based on rate r 
% t is the new time values





% X2 = pvsample(X, t, hop); % X = D?

%% (istft) D > C
% Invert to a waveform


ftsize = n; %second input of istft function
h = hop; %set default hop size for istft

s2 = size(c2);
if s2(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
cols2 = s2(2);

sd7 = size(cd7);
if sd7(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsd7 = sd7(2);

sd5 = size(cd5);
if sd5(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsd5 = sd5(2);

sa5 = size(ca5);
if sa5(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsa5 = sa5(2);
% BREAKER

s = size(c);
if s(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
cols = s(2);

sm = size(cm);
if sm(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsm = sm(2);

s4 = size(c4);
if s4(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
cols4 = s4(2);

s5 = size(c5);
if s5(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
cols5 = s5(2);

sm7 = size(cm7);
if sm7(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsm7 = sm7(2);

sM7 = size(cM7);
if sM7(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
colsM7 = sM7(2);

s8 = size(c8);
if s8(1) ~= (ftsize/2)+1
  error('number of rows should be fftsize/2+1')
end
cols8 = s8(2);



if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,ftsize);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = ftsize/2;
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, ftsize);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    % 2009-01-06: Make stft-istft loop be identity for 25% hop
    win = 2/3*win;
  end
else
  win = w;
end

w = length(win);

% now can set default hop
if h == 0 
  h = floor(w/2);
end

xlend7 = ftsize + (colsd7-1)*h;
xd7 = zeros(1,xlend7);


xlen2 = ftsize + (cols2-1)*h;
x2 = zeros(1,xlen2);

xlend5 = ftsize + (colsd5-1)*h;
xd5 = zeros(1,xlend5);

xlena5 = ftsize + (colsa5-1)*h;
xa5 = zeros(1,xlena5);
%

xlen = ftsize + (cols-1)*h;
x = zeros(1,xlen);

xlenm = ftsize + (colsm-1)*h;
xm = zeros(1,xlenm);

xlen4 = ftsize + (cols4-1)*h;
x4 = zeros(1,xlen4);

xlen5 = ftsize + (cols5-1)*h;
x5 = zeros(1,xlen5);

xlenm7 = ftsize + (colsm7-1)*h;
xm7 = zeros(1,xlenm7);

xlenM7 = ftsize + (colsM7-1)*h;
xM7 = zeros(1,xlenM7);

xlen8 = ftsize + (cols8-1)*h;
x8 = zeros(1,xlen8);


for b = 0:h:(h*(colsd7-1))
  ft = cd7(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xd7((b+1):(b+ftsize)) = xd7((b+1):(b+ftsize))+px.*win;
end;



for b = 0:h:(h*(cols2-1))
  ft = c2(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x2((b+1):(b+ftsize)) = x2((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(colsd5-1))
  ft = cd5(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xd5((b+1):(b+ftsize)) = xd5((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(colsa5-1))
  ft = ca5(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xa5((b+1):(b+ftsize)) = xa5((b+1):(b+ftsize))+px.*win;
end;


for b = 0:h:(h*(cols-1))
  ft = c(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x((b+1):(b+ftsize)) = x((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(colsm-1))
  ft = cm(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xm((b+1):(b+ftsize)) = xm((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(cols4-1))
  ft = c4(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x4((b+1):(b+ftsize)) = x4((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(cols5-1))
  ft = c5(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x5((b+1):(b+ftsize)) = x5((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(colsm7-1))
  ft = cm7(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xm7((b+1):(b+ftsize)) = xm7((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(colsM7-1))
  ft = cM7(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  xM7((b+1):(b+ftsize)) = xM7((b+1):(b+ftsize))+px.*win;
end;

for b = 0:h:(h*(cols8-1))
  ft = c8(:,1+b/h)';
  ft = [ft, conj(ft([((ftsize/2)):-1:2]))];
  px = real(ifft(ft));
  x8((b+1):(b+ftsize)) = x8((b+1):(b+ftsize))+px.*win;
end;


% y = istft(X2, n, n, hop)';
% function x = istft(d, ftsize, w, h)

%% Resampling and Output (Chord Ratios) Rounded values!!!
voiced7 = resample(xd7,3,5);
voice2 = resample(x2,8,9);
voiced5 = resample(xd5,7071,10000);
voicea5 = resample(xa5,5,8);
voiceM3 = resample(x,4,5); 
%voicem3 = resample(xm,261.6256,311.1270);
voicem3 = resample(xm,5,6); %similar ratio to above
voice4 = resample(x4,3,4);
voice5 = resample(x5,2,3);
%voiceM7 = resample(xM7,261.6256,493.8833);
voiceM7 = resample(xM7,8,15);
%voicem7 = resample(xm7,261.6256,466.1638);
voicem7 = resample(xm7,9,16);
voice8 = resample(x8,1,2);

%% Output needs to be length of shortest if they all are same duration
%M3 records as the shortest version so this keeps all durations the same,
%and the +1 is to stop indices from reaching below 0

chordType = 'empty'; % What's the desired chord type?
while (ismember(chordType, 'empty')) %while the user has not selected MAJOR or minor
    chordType = input('What type of chord would you like to hear? ','s');
end

if ismember(chordType, {'a5','aug5','augmented triad','aug triad','Augmented Triad','AUGTRIAD','Augmented triad','augmentedtriad'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voiceM3(end-length(voice2)+1:end),Fs);
soundsc(voicea5(end-length(voice2)+1:end),Fs);

% elseif ismember(chordType, {'d7','dim7','diminished seventh','dim seventh','Dimished Seventh','DIMSEVENTH','Diminished seventh','diminishedseventh'})
% soundsc(voice(end-length(voice)+1:end),Fs);
% soundsc(voicem3(end-length(voice)+1:end),Fs);
% soundsc(voiced5(end-length(voice)+1:end),Fs);


elseif ismember(chordType, {'d5','dim5','diminished triad','dim triad','Dimished Triad','DIMTRIAD','Diminished triad','diminishedtriad'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voicem3(end-length(voice2)+1:end),Fs);
soundsc(voiced5(end-length(voice2)+1:end),Fs);


elseif ismember(chordType, {'sus2','Sus2','Sus 2','sus 2','suspended 2','suspended second','SUS2','SUS 2'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voice2(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'M3','Major3','MAJOR3','major3','M 3','Major 3','MAJOR 3','major 3'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voiceM3(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'m3','Minor3','MINOR3','minor3','m 3','Minor 3','MINOR 3','minor 3'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voicem3(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'p4','P4','p 4','P 4','perfect4','Perfect4','PERFECT4','Perfect 4','4'})
soundsc(voice(end-length(voiceM3)+1:end),Fs);
soundsc(voice4(end-length(voiceM3)+1:end),Fs);

elseif ismember(chordType, {'p5','P5','p 5','P 5','perfect5','Perfect5','PERFECT5','Perfect 5','5'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'major','Major','MAJOR','major triad','M','Major Triad','MAJOR TRIAD'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voiceM3(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'minor','Minor','MINOR','minor triad','m','Minor Triad','MINOR TRIAD'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voicem3(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'m7', 'Minor7','Minor 7','minor7', 'minor 7'}) 
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voicem3(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);
soundsc(voicem7(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'Major7', 'M7', 'Major 7' , 'M 7' , 'MAJOR7', 'MAJOR 7'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voiceM3(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);
soundsc(voiceM7(end-length(voice2)+1:end),Fs);


elseif ismember(chordType, {'sus4', 'Suspended 4', 'sus 4', 'Sus 4', 'SUS4', 'SUS 4'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voice4(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);


elseif ismember(chordType, {'half dim 7', 'half diminished 7th', 'Half Dim 7', 'halfdim7', 'Halfdim7', 'HALFDIM7'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voicem3(end-length(voice2)+1:end),Fs);
soundsc(voiced5(end-length(voice2)+1:end),Fs);
soundsc(voicem7(end-length(voice2)+1:end),Fs);

elseif ismember(chordType, {'dim 7', 'diminished 7th', 'Dim 7', 'dim7', 'Dim7', 'DIM7'})
soundsc(voice(end-length(voiced7)+1:end),Fs);
soundsc(voicem3(end-length(voiced7)+1:end),Fs);
soundsc(voiced5(end-length(voiced7)+1:end),Fs);

soundsc(voiced7(end-length(voiced7)+1:end),Fs);


elseif ismember(chordType, {'dom 7', 'dominant 7th', 'Dom 7', 'dom7', 'Dom7', 'DOM7'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voiceM3(end-length(voice2)+1:end),Fs);
soundsc(voice5(end-length(voice2)+1:end),Fs);
soundsc(voicem7(end-length(voice2)+1:end),Fs);


else ismember(chordType, {'Octave', 'OCTAVE', '8', 'p8', 'P8', 'octave'})
soundsc(voice(end-length(voice2)+1:end),Fs);
soundsc(voice8(end-length(voice2)+1:end),Fs);
end 



%% Example for Presentation

soundsc(voice(end-length(voiceM3)+1:end),Fs);
pause(3)
%and an octave
soundsc(voice(end-length(voiceM3)+1:end),Fs);
soundsc(voice8(end-length(voiceM3)+1:end),Fs);
pause(3)
%and a major 7th chord
soundsc(voice(end-length(voiceM3)+1:end),Fs);
soundsc(voicem3(end-length(voiceM3)+1:end),Fs);
soundsc(voice5(end-length(voiceM3)+1:end),Fs);
soundsc(voicem7(end-length(voiceM3)+1:end),Fs);






%% Available Notes
soundsc(voice(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voicem3(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voiceM3(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voice4(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voice5(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voicem7(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voiceM7(end-length(voiceM3)+1:end),Fs);
pause(.5)
soundsc(voice8(end-length(voiceM3)+1:end),Fs);

%%
soundsc(voice(end-length(voiceM3)+1:end),Fs);
soundsc(voicem3(end-length(voiceM3)+1:end),Fs);
soundsc(voice5(end-length(voiceM3)+1:end),Fs);
%%
% pause(3)
% 
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% pause(.5)
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voice4(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voicem3(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% pause(.5)
% soundsc(voice8(end-length(voiceM3)+1:end),Fs);
% pause(.25)
% soundsc(voicem7(end-length(voiceM3)+1:end),Fs);
% pause(.5)
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% soundsc(voiceM3(end-length(voiceM3)+1:end),Fs);
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% soundsc(voiceM7(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% soundsc(voice8(end-length(voiceM3)+1:end),Fs);


%%
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% soundsc(voice4(end-length(voiceM3)+1:end),Fs);
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% pause(1)
% soundsc(voice(end-length(voiceM3)+1:end),Fs);
% soundsc(voiceM3(end-length(voiceM3)+1:end),Fs);
% soundsc(voice5(end-length(voiceM3)+1:end),Fs);
% pause(.5)
% soundsc(voice8(end-length(voiceM3)+1:end),Fs);

%% OLD
% if strcmp(chordType, 'M3' | 'Major3' | 'MAJOR3')
%     print('hi')
% end
% if chordType == 'M3' | chordType == 'Major3' | chordType == 'MAJOR3' | chordType == 'major3' | chordType == 'M 3' | chordType == 'Major 3' | chordType == 'MAJOR 3' | chordType == 'major 3'
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voiceM3(1:length(voice8)),Fs);
% elseif chordType == 'm3' | chordType == 'Minor3' | chordType == 'MINOR3' | chordType == 'minor3' | chordType == 'm 3' | chordType == 'Minor 3' | chordType == 'MINOR 3' | chordType == 'minor 3'
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voicem3(1:length(voice8)),Fs);
% elseif chordType == 'p4' | chordType == 'p 4' | chordType == 'P4' | chordType == 'P 4' | chordType == 'perfect4' | chordType == 'Perfect 4' | chordType == 'perfect 4' | chordType == 'Perfect4' 
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voice4(1:length(voice8)),Fs);
% elseif chordType == 'p5' | chordType == 'p 5' | chordType == 'P5' | chordType == 'P 5' | chordType == 'perfect5' | chordType == 'Perfect 5' | chordType == 'perfect 5' | chordType == 'Perfect5' 
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voice5(1:length(voice8)),Fs);
% elseif chordType == 'M' | chordType == 'Major' | chordType == 'MAJOR' | chordType == 'major' | chordType == 'M triad' | chordType == 'Major Triad' | chordType == 'MAJOR TRIAD' | chordType == 'major triad'
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voiceM3(1:length(voice8)),Fs);
% soundsc(voice5(1:length(voice8)),Fs);
% elseif chordType == 'm' | chordType == 'Minor' | chordType == 'MINOR' | chordType == 'minor' | chordType == 'm triad' | chordType == 'Minor' | chordType == 'MINOR TRIAD' | chordType == 'minor triad'
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voicem3(1:length(voice8)),Fs);
% soundsc(voice5(1:length(voice8)),Fs);
% elseif chordType == 'm7' | chordType == 'Minor7' | chordType == 'Minor 7' |chordType == 'minor7' | chordType == 'minor 7' 
% soundsc(voice(1:length(voiceM3)),Fs);
% soundsc(voicem3(1:length(voiceM3)),Fs);
% soundsc(voice5(1:length(voiceM3)),Fs);
% soundsc(voicem7(1:length(voiceM3)),Fs);
% elseif chordType == 'Major7' | chordType == 'M7' | chordType == 'Major 7' | chordType == 'M 7' | chordType == 'MAJOR7' | chordType == 'MAJOR 7'
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voiceM3(1:length(voice8)),Fs);
% soundsc(voice5(1:length(voice8)),Fs);
% soundsc(voiceM7(1:length(voice8)),Fs);
% else chordType == 'Octave' | chordType == 'OCTAVE' | chordType == '8' | chordType == 'p8' | chordType == 'P8' | chordType == 'octave' 
% soundsc(voice(1:length(voice8)),Fs);
% soundsc(voice8(1:length(voice8)),Fs);
% end 
% 
% 
% 
% 
% 
% 
% %pause(.02)
% %soundsc(voice(1:length(voice3rd)),Fs);
% %pause(.3)
% %soundsc(voiceM3rd(1:length(voice3rd)),Fs);
% %pause(3)
% %soundsc(voice3rd(1:length(voice3rd))+voice(1:length(voice3rd)),Fs);
% %soundsc(voice5th(1:length(voice8ve)),Fs);
% %pause(1)
%soundsc(voice8ve(1:length(voice8ve)),Fs);

soundsc(voice(end-length(voiced7)+1:end),Fs);
pause(.5)
soundsc(voicem3(end-length(voiced7)+1:end),Fs);
pause(.5)
soundsc(voiced5(end-length(voiced7)+1:end),Fs);
pause(.5)
soundsc(voiced7(end-length(voiced7)+1:end),Fs);



%% For Presentation
soundsc(voice(end-length(voiced7)+1:end),Fs);
    pause(1.5)
    %%
    
soundsc(voice(end-length(voiced7)+1:end),Fs);
soundsc(voicem3(end-length(voiced7)+1:end),Fs);
    pause(1.5)
    %%
    
soundsc(voice(end-length(voiced7)+1:end),Fs);
soundsc(voiceM3(end-length(voiced7)+1:end),Fs);
soundsc(voice5(end-length(voiced7)+1:end),Fs);
    pause(3)
soundsc(voice(end-length(voiced7)+1:end),Fs);
soundsc(voicem3(end-length(voiced7)+1:end),Fs);
soundsc(voiced5(end-length(voiced7)+1:end),Fs);
soundsc(voiced7(end-length(voiced7)+1:end),Fs);
    pause(1.5)
    %%
soundsc(voice(end-length(voiced7)+1:end),Fs);
soundsc(voiceM3(end-length(voiced7)+1:end),Fs);
soundsc(voice5(end-length(voiced7)+1:end),Fs);
soundsc(voice8(end-length(voiced7)+1:end),Fs);







