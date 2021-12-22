function InhomogenousPoissonTrain = IPT_forModel(NSpines, SpikeS, duration, F, strength)
% Written by Nikolai Dembrow, latest version 2021
% This function generates a series of random (poisson distributed) 
% stimulus/spike times for the desired number of inputs. Note that output
% times are in milliseconds.
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The inputs-->
% NSpines: number of spines/presynaptic neurons you are generating
% SpikeS: Mean Rate (Hz)
% duration: in **seconds***
% F: carrier frequency (Hz), for sinuisoidal oscillation around the mean  
% strength: amplitude +/- of the sinuisoidal oscillation of the carrier
% around the mean rate (Hz)
% The outputs--> write two text files: 
%          1. ('FreqOut.txt') of the Mean Input Frequency for the duration 
%              requested, sampled at 20 KHz
%          2. ('StimRun.txt'), a two-column text file with each stimulus 
%              time (**in milliseconds**) and the input # the stim is from 
%   Also make an array of the output with the same numbers in two called 
%    "InhomogenousPoissonTrain" and plot ISIs 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

meanrate = SpikeS;

rng ('shuffle');

dt = 5e-5; % this model assumes a sampling interval of 20KHz
carrierfreq = F;
carrierstength = meanrate*(strength/SpikeS);
% the constraint of the experimental device (for model this is set to the
% sampling interval), and inactivated below

minmove = 5e-5; 

times = dt:dt:duration;
freq = meanrate+(carrierstength*sin(carrierfreq*2*pi*times));

% make a figure, plot the mean rate
figure(1) 
plot (times,freq)
hold on

% make empty arrays to deposit values in
spikebank=[]; 
ISIbank = [];
tstart = [];

for train = 1: NSpines
    vt = rand(size(times));     % make random numbers for each spot in time

% One of the easiest ways to generate a Poisson spike train is to rely 
% on the approximation: 
% Prob {1 spike during (t, t + timeStepS) } = r(t) * timeStepS 
% where r(t) is the instantaneous firing rate (in spikes per second) and timeStepS  
% is the time step. This assumes the probability of more than one spike in
% the time frame is non-existent (safe for dt < 1 ms)
    
    spikes = (freq*dt) > vt;      
    slots = find (spikes);
    refrspikes = spikes;
    % zero out the refractory spikes
        for delete=1:numel(slots)
            counter = slots(delete)+1;
            refrspikes(counter:counter+39) = 0;  % refractory period is set to 40*.05 = 2 ms
        end
    
    TotalSpikes = sum(spikes);
    TotalRefrSpikes = sum(refrspikes);
    spiketimes = dt*(find(spikes));
    refrspiketimes = (dt*(find(refrspikes)))';
    spineref = zeros(length(refrspiketimes),1);
    spineref = spineref+train;
    thistrain = [spineref, refrspiketimes];  
    spikebank = [spikebank; thistrain];
    
    rasters = ones(size(refrspiketimes))*max(freq);   
    ISIs = spiketimes(2:end)-spiketimes(1:end-1);
    rISIs =refrspiketimes(2:end)-refrspiketimes(1:end-1);
    spineref = (zeros(length(rISIs),1))+train;
    thistrain = [spineref, rISIs];
    ISIbank = [ISIbank; thistrain];
    
end


    ticMargin = 0.01;
    ticHeight = ( ((max(freq, [], 'all')/NSpines)) * (1-ticMargin));

    for i=1:length(spikebank)
   
        yOffset = ticMargin + ((spikebank(i,1)-1)*ticHeight); 
        line([(spikebank(i,2)), (spikebank(i,2))], [yOffset, yOffset + ticHeight],'Color', 'red','LineWidth',1);   
    end
    
    tstart = spikebank;
    
    if (length (tstart)==0)
        SpikeOutput=[];
        SpineOutput=[];
        adjmean_rate=[];
        adjCV_int=[];
        mean_count=[];
        var_count=[];
        return
    end
    
    % Now adjust the event times such that they work within the constraint
    % of the experimental device (minimum movement can be set minmove = 5e-3 above)
    
    [~, arraypos] = sort(tstart(:,2));  % sort the spikes in time
    SortedSpike = tstart(arraypos,:);
    StimBank = SortedSpike;
    StimTimes = SortedSpike(:,2);    
    StimTimes =  dt * round (StimTimes/dt);
    StimIntervals = (StimTimes (2:end)-StimTimes(1:end-1));
    fprintf ('number with too small of an interval % d\n', length (find(StimIntervals<minmove)));
    fprintf ('out of % d\n', length(StimIntervals)); % length (find(StimIntervals<minmove))/
    AdjSpikeTimes = StimTimes;
    StimIntervals = [0 ; StimIntervals(:)];
    
    crib = find(StimIntervals<minmove); % all the ISIs that are too small 
    IntervalstoChange = StimIntervals(crib);
    Difference = (ones(length(crib),1) *minmove)- StimIntervals(crib);
    Difference = dt * round(Difference/dt);
    Changeit = StimIntervals<minmove;
    
    spikeIntervals = [];
    adjspikeIntervals = [];
    
        for scan = 2:1:length(crib)
            shifter = Difference(scan)+dt;
          %  if (minmove == dt)
          %  shifter = 0; 
          %  end
            AdjSpikeTimes((crib(scan)):end) = AdjSpikeTimes((crib(scan)):end)+shifter;
        end

    AdjSpikeTimes = dt * round(AdjSpikeTimes /dt);
    adjspikeIntervals = (AdjSpikeTimes (2:end)-AdjSpikeTimes(1:end-1));
    fprintf ('number with too small of an interval after adjustment % d\n', length (find(adjspikeIntervals<minmove)));
    
    AdjSet = [StimBank(:,1), AdjSpikeTimes];
    [~, arraypos] = sort(AdjSet(:,1)); 
    SortbySpine = AdjSet(arraypos,:);
    adjspikeIntervals = [0; adjspikeIntervals];
    cribtoo = adjspikeIntervals<minmove; 
        
    NewStimTimes = SortbySpine(:,2);
    NewStimIntervals = SortbySpine(2:length(SortbySpine),:)- SortbySpine(1:(length(SortbySpine)-1),:);

 % plot the "shifted intervals" in blue 
    for i=1:length(spikebank)
        yOffset = ticMargin + ((AdjSet(i,1)-1)*ticHeight); 
        line([(AdjSet(i,2)), (AdjSet(i,2))], [yOffset, yOffset + ticHeight],'LineWidth',1, 'Color', 'b'); 
    end

    
    Negate = (find(NewStimIntervals(:,1)~=0));
    NewStimIntervals(Negate,2) = 0;
    
    adjISIs = (NewStimIntervals(:,2));
    
    newISIs = adjISIs(adjISIs~=0);
    
    ISIs = ISIbank(:,2);
    ISIs = dt * round(ISIs / dt);
    newISIs = dt * round(ISIs / dt);
    hold off
    mean_int = mean(ISIs(:));
    mean_rate =  1./mean_int;
    adjmean_int = mean(newISIs(:));
    adjmean_rate = 1./adjmean_int;
    CV_int = sqrt((var(ISIs(:)))/mean_int);
    adjCV_int = sqrt((var(newISIs(:)))/adjmean_int);
    
    mean_count = mean(TotalSpikes);
    var_count = var(TotalSpikes);

    fprintf ('Mean rate = %d Hz\n', mean_rate);
    fprintf ('Mean rate after adjustments = %d Hz\n', adjmean_rate);
    fprintf ('Interval CV = %d \n', CV_int);
    fprintf ('Interval CV after adjustments = %d \n', adjCV_int);
    fprintf ('Mean count = %d with %d variance\n', mean_count,var_count)
    
% NOW LETs Plot those beautiful Intervals!
    figure(2)
    binSize = .001;                                            % 1 ms bins 
    x = [binSize:binSize:1]; 
    intervalDist = hist(ISIs(ISIs < 1), x);
    intervalDist = intervalDist / sum(intervalDist) / binSize; % normalize by dividing by spike number 
    
   adjintervalDist = hist(newISIs(newISIs < 1), x);
   adjintervalDist = adjintervalDist / sum(adjintervalDist) / binSize;

    
    bar(x, intervalDist); 
    hold on;
    axis([0 (.500) 0 max(intervalDist)*1.1]);
    x=x+.5;
    bar(x, adjintervalDist, 'r');
    axis([0 (.500) 0 max(adjintervalDist)*1.1]);     
    hold off;
    
    InhomogenousPoissonTrain = AdjSet; % 2 column array with spine number
    % and the spike time in ms
    
% NOW LETs write the spike times to a text file!
    FreqOut = [times; freq];
    FreqOut = FreqOut';
    dlmwrite(strcat('StimRun.txt'),AdjSet,'delimiter', ',')
    dlmwrite(strcat('FreqOut.txt'),FreqOut,'delimiter', ',')
    
end