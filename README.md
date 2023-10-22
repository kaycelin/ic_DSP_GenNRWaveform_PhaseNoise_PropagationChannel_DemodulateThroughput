# ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput
## Script
- example_nrWaveformGen_EVM_Throughput.m
- nrWaveformGen.m

## Main Section:
- History
- Test Setup
- Phase Noise Model
- Propagation Channel Model
- PUSCH Waveform Generator
- PDSCH Waveform Generator
- Test Model Waveform Generator
- Apply Phase Noise Model
- Apply Propagation Channel Model
- Assign Power to Waveform
- Apply AWGN
- Plot Waveform
- Apply LNA
- Measurements: EVM and Throughput

## Test Setup
- select waveform format
```js
if 1
    waveformFormat      = 'testModel'; % nr test model waveform
else
    waveformFormat      = 'pdsch'; % nr pdsch waveform
    waveformFormat      = 'pusch'; % nr pusch waveform
end
```
- set channel bandwidth and scs
```js
bwMHz                   = 100;
scskHz                  = 30;
```
- configure sampling rate
```js
sampleRate              = 245.76e6;
carrierFreq_nSubframes_sampleRate_cell  = {0, [], sampleRate};
```
- configure duplex mode, only be used for test model generator
```js
if 1
    dm = 'TDD';
else
    dm = 'FDD';
end
```
- set RF impariment
```js
isPhaseNoise            = 1;
isPropagationChannel    = 1;
isAWGN                  = 0;
isPwr                   = 1;
```
- set RF chain device
```js
isLNA                   = 1;
```
- set number of antennas
```js
NTxAnts = 1;
NRxAnts = 2;
```
- others settings
```js
fnum                    = 0928; % plot
isDebug                 = 0;    % debug mode: demodulation
isSave                  = 1;    % save waveform
```

## Phase Noise Model
- set phase noise model
```js
sim.pn.FcHz     = 30e9;
sim.pn.PNModel  = 'A'; % 'A' (TDoc R1-163984 Set A), 'B' (TDoc R1-163984 Set B), 'C' (TR 38.803)
foffsetLog      = (4:0.1:log10(sampleRate/2)); % Model offset from 1e4 Hz to sr/2 Hz
foffset         = 10.^foffsetLog;         % Linear frequency offset
pn_PSD          = hPhaseNoisePoleZeroModel(foffset,sim.pn.FcHz,sim.pn.PNModel); % dBc/Hz
```
- gen. phase noise model
```js
pnoise              = comm.PhaseNoise('FrequencyOffset',foffset,'Level',pn_PSD,'SampleRate',sampleRate);
pnoise.RandomStream = "mt19937ar with seed";
```
- Visualize spectrum mask of phase noise
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/a399b433-0c08-40c2-800f-cab22c8f6e4f" width="60%">

## Propagation Channel Model
- set propagation channel model
```js
if 1
    ch.DelayProfile = 'CDL';
else
    ch.DelayProfile = 'TDL';
end
```
- Constructed the CDL or TDL channel model object
```js
if contains(ch.DelayProfile,'CDL','IgnoreCase',true)

    channel = nrCDLChannel; % CDL channel object

    % Turn the overall number of antennas into a specific antenna panel
    % array geometry. The number of antennas configured is updated when
    % nTxAnts is not one of (1,2,4,8,16,32,64,128,256,512,1024) or nRxAnts
    % is not 1 or even.
    [channel.TransmitAntennaArray.Size,channel.ReceiveAntennaArray.Size] = ...
        hArrayGeometry(NTxAnts,NRxAnts);
    nTxAnts = prod(channel.TransmitAntennaArray.Size);
    nRxAnts = prod(channel.ReceiveAntennaArray.Size);
    NTxAnts = nTxAnts;
    NRxAnts = nRxAnts;
else
    channel = nrTDLChannel; % TDL channel object

    % Set the channel geometry
    channel.NumTransmitAntennas = NTxAnts;
    channel.NumReceiveAntennas  = NRxAnts;
end
```
- Assign simulation channel parameters and waveform sample rate to the object
```js
channel.DelaySpread         = 300e-9;
channel.MaximumDopplerShift = 5;
channel.SampleRate          = sampleRate;
```
## PUSCH Waveform Generator
## PDSCH Waveform Generator
## Test Model Waveform Generator
- set fixed reference channel
```js
ulnrref = "G-FR1-A5-7"; % FRC
```
- plot the waveform
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/7683c10d-6ec8-4be2-b6a4-e673e2ce0d36" width="60%">


## Apply Phase Noise Model
```js
txWaveform_PN = txWvCls.setPhaseNoise(txWaveform, pnoise, sampleRate);
```
## Apply Propagation Channel Model
```js
txWaveform_Ch = txWvCls.setPropagationChannel(txWaveform, channel, sampleRate);
```
## Assign Power to Waveform
- set power dBm
```js
pwrdBm = -80;
```
## Apply AWGN
- set SNRdB
```js
SNRdB = 20;
```
## Plot Waveform
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/133d4c1e-45d4-43bc-93d5-dd1c9285ef10" width="60%">
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/51b10e9f-fd3b-40f1-a844-a52eed37312a" width="60%">

## Apply LNA 
```js
[rxWaveform_LNA(:,k),pwrdBm_rxWv2,evmLNA] = applyRxLink(rxWaveform(:,k), {}, {sampleRate, bwInband}, [1018,size(rxWaveform,2),1,k]);
```
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/66aefe0f-7be5-43f1-9617-faefaf3070a6" width="60%">

## Measurements: EVM and Throughput
- Throughput results
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/31334c11-3bf1-4e76-9498-e974bb716acb" width="60%">


- EVM results
<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/44ae970d-bc4c-46dc-8c55-472c4c6c687c" width="60%">

<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/f81dece9-6013-4335-ac9e-c31592e3b761" width="60%">

<img src="https://github.com/kaycelin/ic_DSP_GenNRWaveform_PhaseNoise_PropagationChannel_DemodulateThroughput/assets/87049112/a614860a-8618-4994-bdf7-c9732921a892" width="60%">








  
  
