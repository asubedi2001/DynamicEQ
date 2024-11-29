clc; clear; close all

% MATLAB toolkits required:
% - Communications Toolbox
% - Audio Toolbox

% human hearing limits: 20 hz to 20khz
minFrequency = 0;
maxFrequency = 20000;

% how many samples to process at a time
frameLength = 1024;

% audio reader:
songReader = dsp.AudioFileReader('../audio/rag_crop.wav', 'SamplesPerFrame', frameLength);
sampleRate = songReader.SampleRate;
% writer -- outputs to speakers:
deviceWriter = audioDeviceWriter('SampleRate', songReader.SampleRate);

% used to display sound waves while playing:
scope = timescope( ...
    'SampleRate', songReader.SampleRate, ...
    'TimeSpan', 2, ...
    'BufferLength', songReader.SampleRate*2*2, ...
    'YLimits', [-1,1], ...
    'TimeSpanOverrunAction', "Scroll");

% plot
freqFigure = figure();
outputFreqX = zeros(1);
outputFreqY = zeros(1);
outputPlot = plot(outputFreqX, outputFreqY, 'DisplayName', 'Output');
hold on
inputFreqX = zeros(1);
inputFreqY = zeros(1);
inputPlot = plot(inputFreqX, inputFreqY, 'DisplayName', 'Input');
xscale log
title("Frequency Domain Graph")
xlabel("Frequency (Hz)")
ylabel("Amplitude")
outputPlot.XDataSource = 'outputFreqX';
outputPlot.YDataSource = 'outputFreqY';
inputPlot.XDataSource = 'inputFreqX';
inputPlot.YDataSource = 'inputFreqY';
legend
hold off

% code from 
% https://stackoverflow.com/questions/61369933/how-to-build-a-bandpass-filter-in-matlab-with-the-butter-function
% https://www.mathworks.com/matlabcentral/answers/1890802-quickest-way-to-do-a-bandpass-filter
f_low = 1; % Hz
f_high = 800; % Hz
f_sampling = songReader.SampleRate;
f_nrm_low   = f_low /(f_sampling/2);
f_nrm_high  = f_high /(f_sampling/2);
% determine filter coefficients:
[z,p,k] = butter(4,[f_nrm_low f_nrm_high],'bandpass');
% convert to zero-pole-gain filter parameter (recommended)
sos = zp2sos(z,p,k);

inputEq = graphicEQ('Structure','Cascade', 'SampleRate', songReader.SampleRate);
inputEq.Gains = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

gainsFigure = figure();
gains = zeros(10, 1);
gainsPlot = plot(gains, 'DisplayName', 'Gains');
hold on
title("Equalizer Gains")
gainsPlot.YDataSource = 'inputEq.Gains';
%ylim([-1 inf])
hold off

bassOctaves = logical([1, 1, 1, 1, 0, 0, 0, 0, 0, 0]);
midOctaves =  logical([0, 0, 0, 0, 1, 1, 1, 1, 0, 0]);
trebOctaves = logical([0, 0, 0, 0, 0, 0, 0, 0, 1, 1]);
bassFreqRange = [22, 355];
midFreqRange = [355, 5623];
trebFreqRange = [5623, 22387];

bassInOutFreqDiffs = [];
bassOutFreqSums = [];
midInOutFreqDiffs = [];
midOutFreqSums = [];
trebInOutFreqDiffs = [];
trebOutFreqSums = [];

frameCount = 0;
framePeriod = 20;
redrawPeriod = 35;
while ~isDone(songReader) && frameCount < 900 % use frame count for early stopping TODO: remove
    audioFrame = songReader();
    outputFrame = inputEq(audioFrame);
    noiseFrame = sosfilt(sos, wgn(frameLength,1, -0.5));
    %noiseFrame = wgn(frameLength,1,-25);
    inputFrame = outputFrame + noiseFrame;

    deviceWriter(inputFrame);
    scope(inputFrame)

    [outputFreqX, outputFreqY] = perform_and_plot_fft(outputFrame, sampleRate, maxFrequency);
    [inputFreqX, inputFreqY] = perform_and_plot_fft(inputFrame, sampleRate, maxFrequency);

    % bass:
    bassFreqs = outputFreqX(:) >= bassFreqRange(1) & outputFreqX(:) < bassFreqRange(2); % outputFreqX should be same as inputFreqX
    if isempty(bassInOutFreqDiffs)
        bassInOutFreqDiffs = zeros(sum(bassFreqs), 1);
        bassOutFreqSums = zeros(sum(bassFreqs), 1);
    end
    bassInOutFreqDiffs = bassInOutFreqDiffs + (inputFreqY(bassFreqs) - outputFreqY(bassFreqs));
    bassOutFreqSums = bassOutFreqSums + outputFreqY(bassFreqs);

    % mid:
    midFreqs = outputFreqX(:) >= midFreqRange(1) & outputFreqX(:) < midFreqRange(2); % outputFreqX should be same as inputFreqX
    if isempty(midInOutFreqDiffs)
        midInOutFreqDiffs = zeros(sum(midFreqs), 1);
        midOutFreqSums = zeros(sum(midFreqs), 1);
    end
    midInOutFreqDiffs = midInOutFreqDiffs + (inputFreqY(midFreqs) - outputFreqY(midFreqs));
    midOutFreqSums = midOutFreqSums + outputFreqY(midFreqs);

    % treb:
    trebFreqs = outputFreqX(:) >= trebFreqRange(1) & outputFreqX(:) < trebFreqRange(2); % outputFreqX should be same as inputFreqX
    if isempty(trebInOutFreqDiffs)
        trebInOutFreqDiffs = zeros(sum(trebFreqs), 1);
        trebOutFreqSums = zeros(sum(trebFreqs), 1);
    end
    trebInOutFreqDiffs = trebInOutFreqDiffs + (inputFreqY(trebFreqs) - outputFreqY(trebFreqs));
    trebOutFreqSums = trebOutFreqSums + outputFreqY(trebFreqs);


    if mod(frameCount, framePeriod) == 0 && frameCount ~= 0
        disp('Calculating new gains')
        % TODO: Reduce gains if significantly higher than noise, but ensure
        % always >= 0
        % TODO: Set max limits for gains to prevent distortion

        % bass:
        bassDiffs = bassInOutFreqDiffs / framePeriod;
        bassOutSums = bassOutFreqSums / framePeriod;
        if sum(bassDiffs) > sum(bassOutSums)
            disp('Updating bass gains')
            inputEq.Gains(bassOctaves) = inputEq.Gains(bassOctaves) + 1;
            disp(inputEq.Gains)
        end
        bassInOutFreqDiffs(:) = 0;
        bassOutFreqSums(:) = 0;

        % mid
        midDiffs = midInOutFreqDiffs / framePeriod;
        midOutSums = midOutFreqSums / framePeriod;
        if sum(midDiffs) > sum(midOutSums)
            disp('Updating mid gains')
            inputEq.Gains(midOctaves) = inputEq.Gains(midOctaves) + 1;
            disp(inputEq.Gains)
        end
        midInOutFreqDiffs(:) = 0;
        midOutFreqSums(:) = 0;

        % treb:
        trebDiffs = trebInOutFreqDiffs / framePeriod;
        trebOutSums = trebOutFreqSums / framePeriod;
        if sum(trebDiffs) > sum(trebOutSums)
            disp('Updating treb gains')
            inputEq.Gains(trebOctaves) = inputEq.Gains(trebOctaves) + 1;
            disp(inputEq.Gains)
        end
        trebInOutFreqDiffs(:) = 0;
        trebOutFreqSums(:) = 0;

    end

    if mod(frameCount, redrawPeriod) == 0
        refreshdata
        drawnow
    end

    frameCount = frameCount + 1;
end

function [freqX, freqY] = perform_and_plot_fft(signal, sampleRate, maxFrequency)
    numSamples = numel(signal);
    [freqX, freqY, scale] = perform_fft(signal, numSamples, sampleRate);
    cutoff = maxFrequency/scale;
    freqX = freqX(1:floor(cutoff));
    freqY = freqY(1:floor(cutoff));
end

function [freqX, freqY, scale] = perform_fft(songSignal, numSamples, sampleRate)
    % based on https://www.mathworks.com/help/matlabmobile/ug/acquire-and-analyze-audio-data.html
    fftData = fft(songSignal);
    signalLength = numSamples;
    freqY = abs(fftData/signalLength); % abs = sqrt of sum of squares
    freqY = freqY(1:floor(signalLength/2)+1);
    freqY(2:end-1) = 2*freqY(2:end-1);
    freqX = sampleRate*(0:(signalLength/2))/signalLength;
    scale = sampleRate/signalLength;
end

release(songReader)
release(deviceWriter)
release(scope)