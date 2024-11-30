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

% create frame buffer
bufferSize = 8;
frameBuffer = zeros(frameLength, bufferSize);
bufferIdx = 1;

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

bassBand = struct( ...
    'octaves', logical([1, 1, 1, 1, 0, 0, 0, 0, 0, 0]), ...
    'freqRange', [22, 355], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'movingAvg', [] ...
    );

midBand = struct( ...
    'octaves', logical([0, 0, 0, 0, 1, 1, 1, 1, 0, 0]), ...
    'freqRange', [355, 5623], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'movingAvg', [] ...
    );

trebBand = struct( ...
    'octaves', logical([0, 0, 0, 0, 0, 0, 0, 0, 1, 1]), ...
    'freqRange', [5623, 22387], ...
    'inOutFreqDiffs', [], ...
    'outFreqSums', [], ...
    'movingAvg', [] ...4
    );

bands = [bassBand midBand trebBand];

frameCount = 0;
framePeriod = 20;
redrawPeriod = 35;
while ~isDone(songReader) && frameCount < 900 % use frame count for early stopping TODO: remove
    audioFrame = songReader();
    outputFrame = inputEq(audioFrame);
    noiseFrame = sosfilt(sos, wgn(frameLength,1, -0.5));
    %noiseFrame = wgn(frameLength,1,-25);
    inputFrame = outputFrame + noiseFrame;

    % update buffer (and index) according to current frame
    frameBuffer(:, bufferIdx) = inputFrame;
    bufferIdx = mod(bufferIdx, bufferSize) + 1;

    deviceWriter(inputFrame);
    scope(inputFrame)

    [outputFreqX, outputFreqY] = perform_and_plot_fft(outputFrame, sampleRate, maxFrequency);
    [inputFreqX, inputFreqY] = perform_and_plot_fft(inputFrame, sampleRate, maxFrequency);

    for i = 1:numel(bands)
        bandFreqs = outputFreqX(:) >= bands(i).freqRange(1) & outputFreqX(:) < bands(i).freqRange(2); % outputFreqX should be same as inputFreqX
        if isempty(bands(i).inOutFreqDiffs)
            bands(i).inOutFreqDiffs = zeros(sum(bandFreqs), 1);
            bands(i).outFreqSums = zeros(sum(bandFreqs), 1);
            bands(i).movingAvg = zeros(sum(bandFreqs), 1);
        end
        bands(i).inOutFreqDiffs = bands(i).inOutFreqDiffs + (inputFreqY(bandFreqs) - outputFreqY(bandFreqs));
        bands(i).outFreqSums = bands(i).outFreqSums + outputFreqY(bandFreqs);

        % update moving average
        if frameCount == 0
            % initialize first frame value
            bands(i).movingAvg = bands(i).inOutFreqDiffs;
        elseif frameCount <= bufferSize
            % if buffer is not full update b += (difference - b) / numFrames
            bands(i).movingAvg = bands(i).movingAvg + (bands(i).inOutFreqDiffs - bands(i).movingAvg) / frameCount;
        else
            % if buffer full, update according to the last value 
            oldestFrameIdx = mod(bufferIdx, bufferSize) + 1;
            oldestFreqDiff = frameBuffer(:, oldestFrameIdx);
            % b += (difference - (difference bufferSize ago)) / bufferSize
            bands(i).movingAvg = bands(i).movingAvg + (bands(i).inOutFreqDiffs - oldestFreqDiff(bandFreqs)) / bufferSize;
        end
    end

    if mod(frameCount, framePeriod) == 0 && frameCount ~= 0
        disp('Calculating new gains')
        % TODO: Reduce gains if significantly higher than noise, but ensure
        % always >= 0
        % TODO: Set max limits for gains to prevent distortion

        % define minimum and maximum value thresholds for gain adjustment
        minValue = 0; min(cellfun(@min, {bands.movingAvg})); % minimum value moving average has
        maxValue = max(cellfun(@max, {bands.movingAvg})); % maximum value moving average has

        for i = 1:numel(bands)
            % there are probably better / more 'correct' ways to normalize this
            bandMovingAvgNorm = (bands(i).movingAvg - minValue) / (maxValue - minValue);
            bandMovingAvgNorm = max(0, min(1, bandMovingAvgNorm));
            bandGain = bandMovingAvgNorm * 10;
            disp('Updating gains for band')
            bands(i).gain = mean(bandGain);
            disp(bands(i).gain);
            inputEq.Gains(bands(i).octaves) = bands(i).gain;

            bands(i).inOutFreqDiffs(:) = 0;
            bands(i).outFreqSums(:) = 0;
        end
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