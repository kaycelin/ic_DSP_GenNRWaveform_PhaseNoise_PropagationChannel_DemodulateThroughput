%% Histroy
% 2023-09-21, change class name from waveformGen to nrWaveformeGen
% 2023-09-21, refer to https://www.mathworks.com/help/5g/ug/5g-nr-tm-waveform-generation.html
% 2023-09-23, refer to https://www.mathworks.com/help/5g/ug/nr-pusch-throughput.html

%% 5G NR-TM and FRC Waveform Generation
% To generate 5G NR test models (NR-TMs) and uplink and downlink fixed reference channels (FRCs).
% To specify the NR-TM or FRC name, the channel bandwidth, the subcarrier spacing, and the duplexing mode.

%% Introduction
% The NR test models (NR-TM), for the purpose of base station (BS) RF testing
% The downlink fixed reference channels (FRC), for user equipment (UE) input testing.

% The NR-TMs for FR1 are defined in TS 38.141-1 Section 4.9.2, and the NR-TMs for FR2 are defined in TS 38.141-2 Section 4.9.2.
% The cases of base station (BS) RF testing:
% - BS output power
% - Timing alignment error (TAE)
% - Occupied bandwidth emissions
% - Adjacent channel leakage ratio (ACLR)
% - Operating band unwanted emissions
% - Transmitter spurious emissions
% - Transmitter intermodulation

% The physical downlink shared channel (PDSCH) FRC for FR1 are defined in TS 38.101-1 Annex A.3, and for FR2 are defined in TS 38.101-2 Annex A.3.
% The cases of UE test, including:
% - UE receiver requirements
% - Maximum UE input level testing

% The physical uplink shared channel (PUSCH) FRC for FR1 and FR2 are defined in TS 38.104 Annex A.
% The cases of BS RF reception tests, including:
% - Reference sensitivity
% - Adjacent channel selectivity (ACS)
% - In-band and out-of-band blocking
% - Receiver intermodulation
% - In-channel selectivity
% - Dynamic range
% - Performance requirements

% This application uses the MATLAB class hNRReferenceWaveformGenerator.
% hNRReferenceWaveformGenerator class provides properties:
% - fr1bandwidthtable, fr2bandwidthtable:                           Bandwidth configuration table
% - fr1testmodels, fr2testmodels, fr1downlinkfrc, fr2downlinkfrc:   The Release 15 and Release 16 test model and FRC lists
% - generateWaveform:                                               Baseband waveform generation
% - displayResourceGrid:                                            Resource grid visualization

%% Script structure
% Functions:
% - nrWaveformGen: main (class)
% - installParams: intall parameters (priviate)
% - waveformGen:   nr waveform generator (static)


classdef nrWaveformGen < handle

    properties % Test case
        wvFormat = 'pusch';
    end

    properties % configuration
        PUSCH_config
        PDSCH_config
        TM_config
    end

    properties % Carrier settings
        CarrierFreq     = 0;
        NSubframes      = [];
        SampleRate      = [];
        PhaseNoise      = [];
        ChannelModel    = []
    end

    properties(Access=private) % Simulation setup
        fnum            = 0;
        isDebug         = 0;
        isSave          = 0;
    end

    properties % waveform
        Waveform
        WaveformInfo
        CarrierBw
        ChannelBw
        SubcarrierSpacing
    end

    properties % results
        Summary
    end

    methods
        function [obj, waveform, wv_config] = nrWaveformGen(wvFormat, carrier_pxsch_pxschExtension_sim_pn_ch_cell ,ref_bw_scs_dm_ncellid_cell, ...
                carrierFreq_nSubframes_sampleRate_cell, phaseNoise_channel_cell, isTDD, fnum, isDebug, isSave) % 2023-09-21
            if ~exist('isTDD','var')||isempty(isTDD)
                isTDD = 1;
            end
            if ~exist('isDebug','var')||isempty(isDebug)
                isDebug = 0;
            end
            if ~exist('fnum','var')||isempty(fnum)
                isFnum = 0;
                fnum = [];
            else
                isFnum = 1;
            end
            if ~exist('isSave','var')||isempty(isSave)
                isSave = 1;
            end
            tddSeqPlt = [];

            if ~exist('carrier_pxsch_pxschExtension_sim_pn_ch_cell','var')||isempty(carrier_pxsch_pxschExtension_sim_pn_ch_cell)
                carrier_pxsch_pxschExtension_sim_pn_ch_cell = [];
            end
            if ~exist('ref_bw_scs_dm_ncellid_cell','var')||isempty(ref_bw_scs_dm_ncellid_cell)
                ref_bw_scs_dm_ncellid_cell = [];
            end
            if ~exist('carrierFreq_nSubframes_sampleRate_cell','var')||isempty(carrierFreq_nSubframes_sampleRate_cell)
                carrierFreq_nSubframes_sampleRate_cell = [];
            end
            if ~exist('phaseNoise_channel_cell','var')||isempty(phaseNoise_channel_cell)
                phaseNoise_channel_cell = [];
            end

            % install parameters
            obj.installParams(wvFormat, carrier_pxsch_pxschExtension_sim_pn_ch_cell, ref_bw_scs_dm_ncellid_cell, carrierFreq_nSubframes_sampleRate_cell, phaseNoise_channel_cell);

            % generate waveform by No. of carrier
            switch wvFormat
                case {'pdsch'}
                    % generate waveform
                    [waveform, wv_config, obj.WaveformInfo] = obj.pdschTx([], [], [], [], [], fnum);

                    if isDebug*0
                        % demodulate waveform
                        obj.pdschRx(waveform,[],fnum);
                    end

                    % export
                    obj.NSubframes = obj.PDSCH_config.sim.NFrames * 10;

                case {'pusch'}
                    [waveform, wv_config, obj.WaveformInfo] = obj.puschTx;

                    % export
                    obj.NSubframes = obj.PUSCH_config.sim.NFrames * 10;

                    % plot
                    wvInfo = 'pusch-' + string(obj.PUSCH_config.carrier.NSizeGrid) + 'RBs-' +...
                        'scs' + string(obj.PUSCH_config.carrier.SubcarrierSpacing) + 'kHz-' +...
                        string(obj.PUSCH_config.pusch.Modulation) + '-' +...
                        'targetCodeRate ' + string(obj.PUSCH_config.puschExtension.TargetCodeRate);

                case {'TMs','FRCs','UplinkFrc','testModel','TestModel'}
                    [waveform, wv_config, obj.WaveformInfo] = obj.TestModelTx([],[],[],[],[],[],fnum);

                    % debug
                    if isDebug*0
                        obj.TestModelRx(waveform,[],[],fnum);
                    end

                    % export - plot and wvInfo
                    if strcmpi(obj.TM_config.dm,'TDD')
                        nr = nrSpec;
                        nrFrc = nr.nrUlFRC(obj.TM_config.nrRef);
                        fr = convertStringsToChars(obj.TM_config.nrRef);
                        idx = strfind(obj.TM_config.nrRef,"FR");
                        tddSeqPlt = {nrFrc.ScskHz, fr(idx:idx+2)};
                    end
            end

            if ~isempty(obj.PhaseNoise)
                for k=1:size(waveform,2)
                    [waveform(:,k), pnoise, pnInfo] = obj.setPhaseNoise(waveform(:,k), obj.PhaseNoise, obj.SampleRate);
                end
                obj.WaveformInfo = obj.WaveformInfo+"-"+pnInfo;
                obj.PhaseNoise = pnoise;
            end

            if ~isempty(obj.ChannelModel)
                % apply channel model
                [waveform,pathGains,sampleTimes,maxChDelay,chMlInfo] = obj.setPropagationChannel(waveform, obj.ChannelModel, obj.SampleRate);
                obj.WaveformInfo = obj.WaveformInfo+"-"+chMlInfo;
                if 1
                    obj.Summary.channelModelPathGains = pathGains;
                end
            end

            if isTDD
                % For TDD we only care about the DL slots/symbols, but for FDD we exploit all the slots
                usefulWvfm      = (abs(waveform)>eps);
            else
                usefulWvfm = 1:numel(waveform);
            end

            % normalize
            waveform = waveform./(rms(waveform(usefulWvfm)));
            if 0
                % Normalize the waveform to fit the dynamic range of the nonlinearity.
                tmWaveform = tmWaveform/max(abs(tmWaveform),[],'all');
            end

            if isDebug
                % demodulate waveform
                switch wvFormat
                    case 'pdsch'
                        obj.pdschRx(waveform,[],fnum);
                    case 'pusch'
                        obj.puschRx(waveform,[],fnum);
                    case 'TMs'
                        obj.TestModelRx(waveform,[],[],fnum);
                end
            end

            % plot spectrum and time
            if isFnum
                if 1 % plot spectrum
                    fs = obj.SampleRate;

                    fnum_legend_title_color_cell                = {[fnum,2,1,1] , obj.WaveformInfo, 'Spectrum', 'y'};
                    fc_freqSpan_rbw_cell                        = {0, fs, 10e3};
                    unit_winType_cell                           = {'PSD', []};
                    if 1
                        isSemiLogx_positiveF_nyquistZone_cell  = {0,0,0};
                        [spectrumResults] = plotSpectrum(waveform, fs, fnum_legend_title_color_cell,...
                            fc_freqSpan_rbw_cell, unit_winType_cell, isSemiLogx_positiveF_nyquistZone_cell);
                    else
                        waveCfg_sampleRate_cell(k,:)        = {obj.refWavegenClass.Config, fs};
                        carrierBw_chOffset_cell(k,:)        = {[obj.CarrierBw], [obj.ChannelBw]};
                        nCA_chGapBw_isTDD_cell(k,:)         = {1, 0, obj.dm};
                        plotACLR(waveform, waveCfg_sampleRate_cell(k,:), carrierBw_chOffset_cell(k,:), nCA_chGapBw_isTDD_cell(k,:),...
                            fnum_legend_title_color_cell, fc_freqSpan_rbw_cell, unit_winType_cell);
                    end
                end

                if 1 % plot time
                    sweepTime_dt = [obj.NSubframes*1e-3, 1e-6];
                    fc_ZeroSpan_rbw_cell = {0,0};

                    [~, tddUlDlConfigTable] = plotSpectrum(waveform, fs,...
                        {[[fnum,2,1,2], 2, 1, 2], obj.WaveformInfo, 'Time Domain', 'y'}, fc_ZeroSpan_rbw_cell(:,:), unit_winType_cell, ...
                        [], [], [], sweepTime_dt, 'stem', tddSeqPlt,1);
                end
            end

            % export
            obj.Waveform = waveform;

            if isSave
                filename = obj.WaveformInfo;
                fprintf("\nsave to file %s.mat \n", filename)
                try
                    save(filename,'obj')
                catch
                    object = struct(obj);
                    save(filename,'object')
                end
            end
        end

    end % end-methods

    methods(Access=protected)
        function installParams(obj, wvFormat, carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell, ref_bw_scs_dm_ncellid_cell, carrierFreq_nSubframes_sampleRate_cell, phaseNoise_channel_cell, nCA, k)
            try obj.wvFormat = wvFormat; end
            try obj.nCA     = nCA; end
            try obj.CarrierFreq  = carrierFreq_nSubframes_sampleRate_cell{1}; end
            try obj.NSubframes   = carrierFreq_nSubframes_sampleRate_cell{2}; end
            try obj.SampleRate   = carrierFreq_nSubframes_sampleRate_cell{3}; end
            try obj.PhaseNoise   = phaseNoise_channel_cell{1}; end
            try obj.ChannelModel = phaseNoise_channel_cell{2}; end

            switch wvFormat
                case {'pusch','PUSCH'}
                    try obj.PUSCH_config.carrier = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{1}; end
                    try obj.PUSCH_config.pusch = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{2}; end
                    try obj.PUSCH_config.puschExtension = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{3}; end
                    try obj.PUSCH_config.sim = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{4}; end
                    try obj.PUSCH_config.phaseNoise = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{5}; end
                    try obj.PUSCH_config.channel = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{6}; end

                case {'pdsch','PDSCH'}
                    try obj.PDSCH_config.carrier = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{1}; end
                    try obj.PDSCH_config.pdsch = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{2}; end
                    try obj.PDSCH_config.pdschExtension = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{3}; end
                    try obj.PDSCH_config.sim = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{4}; end
                    try obj.PDSCH_config.phaseNoise = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{5}; end
                    try obj.PDSCH_config.channel = carrier_pxsch_pxschExtension_sim_pn_ch_cfgCell{6}; end

                case {'TMs','FRCs','testModel','TestModel'}
                    try obj.TM_config.nrRef        = string(ref_bw_scs_dm_ncellid_cell{1}); end
                    try obj.TM_config.bw           = string(ref_bw_scs_dm_ncellid_cell{2}); end
                    try obj.TM_config.scs         = string(ref_bw_scs_dm_ncellid_cell{3}); end
                    try obj.TM_config.dm           = string(ref_bw_scs_dm_ncellid_cell{4}); end
                    try obj.TM_config.ncellid         = double(ref_bw_scs_dm_ncellid_cell{5}); end
            end
        end
    end % end-methods

    methods(Access=public)
        function [puschWaveform,wv_config_cell,wvInfo, obj] = puschTx(obj, carrier, pusch, puschExtension, sim, sampleRate, channel, fnum)
            if 1 % validate configuration
                if ~exist('carrier','var')||isempty(carrier)
                    carrier = obj.PUSCH_config.carrier;
                elseif ~strcmpi(class(carrier), 'nrCarrierConfig')
                    error('nrWaveformGen.pdschGen, check input of carrier')
                end
                if ~exist('pusch','var')||isempty(pusch)
                    pusch = obj.PUSCH_config.pusch;
                elseif ~strcmpi(class(pusch), 'nrPUSCHConfig')
                    error('nrWaveformGen.puschTx, check input of pusch')
                end
                if ~exist('puschExtensionCfg','var')||isempty(puschExtension)
                    puschExtension = obj.PUSCH_config.puschExtension;
                elseif ~isfield(puschExtension, 'TargetCodeRate')
                    error('nrWaveformGen.pdschGen, check input of pdschExtension.TargetCodeRate')
                end
                if ~exist('sim','var')||isempty(sim)
                    sim = obj.PUSCH_config.sim;
                elseif ~isfield(sim, 'NFrames')
                    error('nrWaveformGen.puschGen, check input of sim.NFrames')
                end
                if ~exist('sampleRate','var')||isempty(sampleRate)
                    sampleRate = obj.SampleRate;
                end
                if ~exist('fnum','var')||isempty(fnum)
                    isFnum = 0;
                    fnum = [];
                else
                    isFnum = 1;
                end
            end

            % Set up redundancy version (RV) sequence for all HARQ processes
            if puschExtension.EnableHARQ
                % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
                % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
                % should be used. If self-decodability is the priority, it should be
                % taken into account that the upper limit of the code rate at which
                % each RV is self-decodable is in the following order: 0>3>2>1
                rvSeq = [0 2 3 1];
            else
                % HARQ disabled - single transmission with RV=0, no retransmissions
                rvSeq = 0;
            end

            % Create UL-SCH encoder system object to perform transport channel encoding
            encodeULSCH = nrULSCH;
            if 1 % install ulsch
                encodeULSCH.MultipleHARQProcesses = true;
                encodeULSCH.TargetCodeRate = puschExtension.TargetCodeRate;
            end

            % Create UL-SCH decoder System object to perform transport channel decoding
            % Use layered belief propagation for LDPC decoding, with half the number of
            % iterations as compared to the default for belief propagation decoding
            if strcmpi(obj.wvFormat,'pusch')
                decodeULSCH = nrULSCHDecoder;
                decodeULSCH.MultipleHARQProcesses = true;
                decodeULSCH.TargetCodeRate = puschExtension.TargetCodeRate;
                decodeULSCH.LDPCDecodingAlgorithm = puschExtension.LDPCDecodingAlgorithm;
                decodeULSCH.MaximumLDPCIterationCount = puschExtension.MaximumLDPCIterationCount;
            end

            if 0
                % Get information about the baseband waveform after OFDM modulation step
                wvInfo = nrOFDMInfo(carrier,'SampleRate',sampleRate);
                waveformInfo_original = nrOFDMInfo(carrier);
            end

            % Specify the fixed order in which we cycle through the HARQ process IDs
            harqSequence = 0:puschExtension.NHARQProcesses-1;

            % Initialize the state of all HARQ processes
            if strcmpi(obj.wvFormat,'pusch')
                harqEntity = HARQEntity(harqSequence,rvSeq);
            else
                harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);
            end

            % Total number of slots in the simulation period
            NSlots = sim.NFrames * carrier.SlotsPerFrame;
            NSlotSymbs = carrier.SymbolsPerSlot;
            NSymbSubcars = 12;
            NRBs = carrier.NSizeGrid;

            % Initialize the grid for specified number of frames
            txGrid = zeros(NRBs*NSymbSubcars, NSlots*NSlotSymbs, pusch.NumLayers);

            if 0 % Initialization cell
                dmrsSymCell = cell(1,NSlots);
                dmrsIndCell = cell(1,NSlots);
                ptrsSymbolsCell = cell(1,NSlots);
                ptrsIndicesCell = cell(1,NSlots);
            end
            rng('default')
            for nSlot = 0:NSlots-1
                % Update the carrier slot numbers for new slot
                carrier.NSlot = nSlot;

                % Calculate the transport block sizes for the transmission in the slot
                [puschIndices,puschIndicesInfo] = nrPUSCHIndices(carrier,pusch);
                trBlkSizes = nrTBS(pusch.Modulation,pusch.NumLayers,numel(pusch.PRBSet),puschIndicesInfo.NREPerPRB,puschExtension.TargetCodeRate,puschExtension.XOverhead);

                % HARQ processing
                if strcmpi(obj.wvFormat,'pusch')
                    % If new data for current process and codeword then create a new UL-SCH transport block
                    if harqEntity.NewData
                        trBlk = randi([0 1],trBlkSizes,1);
                        setTransportBlock(encodeULSCH,trBlk,harqEntity.HARQProcessID);
                        if 1 % 2023-10-04, removed
                            % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
                            if harqEntity.SequenceTimeout
                                resetSoftBuffer(decodeULSCHLocal,harqEntity.HARQProcessID);
                            end
                        end
                    end
                else
                    % If new data for current process and codeword then create a new DL-SCH transport block
                    if harqEntity.NewData(cwIdx)
                        trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                        setTransportBlock(encodeULSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                        if 0 % 2023-10-04, removed
                            % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
                            if harqEntity.SequenceTimeout(cwIdx)
                                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                            end
                        end
                    end
                end

                % Encode the UL-SCH transport blocks
                codeTrBlocks = encodeULSCH(pusch.Modulation,pusch.NumLayers, ...
                    puschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

                % Create resource grid for a slot
                puschGrid = nrResourceGrid(carrier, sim.NTxAnts);

                % PUSCH modulation and precoding
                puschSymbols = nrPUSCH(carrier, pusch, codeTrBlocks);
                if strcmpi(obj.wvFormat,'pusch')
                    % Implementation-specific PUSCH MiMO precoding and
                    % mapping. This MIMO precoding step is in addition to
                    % any codebook based MIMO precoding done during PUSCH
                    % modulation above
                    if strcmpi(pusch.TransmissionScheme,'codebook')
                        % Codebook based MIMO precoding, F precodes between
                        % PUSCH transmit antenna ports and transmit
                        % antennas
                        F = eye(pusch.NumAntennaPorts,sim.NTxAnts);
                    else
                        % Non-codebook based MIMO precoding, F precodes
                        % between PUSCH layers and transmit antennas
                        F = eye(pusch.NumLayers,sim.NTxAnts);
                    end
                    [~,puschAntIndices] = nrExtractResources(puschIndices,puschGrid);
                    puschGrid(puschAntIndices) = puschSymbols * F;
                else
                    pdschGrid(puschIndices) = pdschSymbols;
                end

                % PUSCH DM-RS precoding and mapping
                dmrsSymbols = nrPUSCHDMRS(carrier,pusch);
                dmrsIndices = nrPUSCHDMRSIndices(carrier,pusch);
                if strcmpi(obj.wvFormat,'pusch')
                    for p = 1:size(dmrsSymbols,2)
                        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
                        puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
                    end
                else
                    pdschGrid(dmrsIndices) = dmrsSymbols;
                end

                if strcmpi(obj.wvFormat,'pdsch')
                    % PDSCH PT-RS precoding and mapping
                    ptrsSymbols = nrPDSCHPTRS(carrier,pusch);
                    ptrsIndices = nrPDSCHPTRSIndices(carrier,pusch);
                    pdschGrid(ptrsIndices) = ptrsSymbols;
                end

                if 0
                    % export reference gird with DM-RS symbols into Cell
                    dmrsSymCell{nSlot+1} = dmrsSymbols;
                    dmrsIndCell{nSlot+1} = dmrsIndices;
                    ptrsSymbolsCell{nSlot+1} = ptrsSymbols;
                    ptrsIndicesCell{nSlot+1} = ptrsIndices;
                end

                % Creat txGrid by pdschGrid
                txGrid(:, nSlot*carrier.SymbolsPerSlot+1: (nSlot+1)*carrier.SymbolsPerSlot) = puschGrid;
            end

            % OFDM modulation
            puschWaveform = nrOFDMModulate(carrier, txGrid, 'SampleRate', sampleRate);

            % export - configuration
            wv_config_cell = {carrier,pusch,puschExtension,sim};

            % export to object - txGrid for demodulation
            obj.Summary.txGrid = txGrid;
            obj.Waveform = puschWaveform;
            obj.SubcarrierSpacing = carrier.SubcarrierSpacing;
            obj.CarrierBw = carrier.NSizeGrid * obj.SubcarrierSpacing * 1e3 * 12;
            obj.ChannelBw;

            % export - description
            waveformInfo = nrOFDMInfo(carrier,'SampleRate',sampleRate);
            wvInfo = ['pusch-'+string(carrier.NSizeGrid)+'RBs-'+'scskHz'+string( obj.SubcarrierSpacing)+'-fs'+string(waveformInfo.SampleRate/1e6)+'MHz-'+string(sim.NFrames)+'Frames'];
            wvInfo = strrep(wvInfo,'.','p');

            if isFnum
                plotSpectrum(puschWaveform, sampleRate, {fnum,wvInfo,'Spectrum'});
            end
        end
        function [txWaveform,wv_config_cell,wvInfo, obj] = pdschTx(obj, carrier, pdsch, pdschExtension, sim, sampleRate, channel, fnum)
            if 1 % validate configuration
                if ~exist('carrier','var')||isempty(carrier)
                    carrier = obj.PDSCH_config.carrier;
                elseif ~strcmpi(class(carrier), 'nrCarrierConfig')
                    error('nrWaveformGen.pdschGen, check input of carrier')
                end
                if ~exist('pdsch','var')||isempty(pdsch)
                    pdsch = obj.PDSCH_config.pdsch;
                elseif ~strcmpi(class(pdsch), 'nrPDSCHConfig')
                    error('nrWaveformGen.pdschGen, check input of pdsch')
                end
                if ~exist('pdschExtensionCfg','var')||isempty(pdschExtension)
                    pdschExtension = obj.PDSCH_config.pdschExtension;
                elseif ~isfield(pdschExtension, 'TargetCodeRate')
                    error('nrWaveformGen.pdschGen, check input of pdschExtension.TargetCodeRate')
                end
                if ~exist('sim','var')||isempty(sim)
                    sim = obj.PDSCH_config.sim;
                elseif ~isfield(sim, 'NFrames')
                    error('nrWaveformGen.pdschGen, check input of sim.NFrames')
                end
                if ~exist('sampleRate','var')||isempty(sampleRate)
                    sampleRate = obj.SampleRate;
                end
                if ~exist('fnum','var')||isempty(fnum)
                    isFnum = 0;
                    fnum = [];
                else
                    isFnum = 1;
                end
            end

            % Set up redundancy version (RV) sequence for all HARQ processes
            if pdschExtension.EnableHARQ
                % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
                % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
                % should be used. If self-decodability is the priority, it should be
                % taken into account that the upper limit of the code rate at which
                % each RV is self-decodable is in the following order: 0>3>2>1
                rvSeq = [0 2 3 1];
            else
                % HARQ disabled - single transmission with RV=0, no retransmissions
                rvSeq = 0;
            end

            % Create DL-SCH encoder system object to perform transport channel encoding
            encodeDLSCH = nrDLSCH;
            if 1 % install dlsch
                encodeDLSCH.MultipleHARQProcesses = true;
                encodeDLSCH.TargetCodeRate = pdschExtension.TargetCodeRate;
            end

            if 0
                % Get information about the baseband waveform after OFDM modulation step
                wvInfo = nrOFDMInfo(carrier,'SampleRate',sampleRate);
                waveformInfo_original = nrOFDMInfo(carrier);
            end

            % Specify the fixed order in which we cycle through the HARQ process IDs
            harqSequence = 0:pdschExtension.NHARQProcesses-1;

            % Initialize the state of all HARQ processes
            harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);

            % Total number of slots in the simulation period
            NSlots = sim.NFrames * carrier.SlotsPerFrame;
            NSlotSymbs = carrier.SymbolsPerSlot;
            NSymbSubcars = 12;
            NRBs = carrier.NSizeGrid;

            % Initialize the grid for specified number of frames
            txGrid = zeros(NRBs*NSymbSubcars, NSlots*NSlotSymbs, pdsch.NumLayers);

            if 0 % Initialization cell
                dmrsSymCell = cell(1,NSlots);
                dmrsIndCell = cell(1,NSlots);
                ptrsSymbolsCell = cell(1,NSlots);
                ptrsIndicesCell = cell(1,NSlots);
            end
            rng('default')
            for nSlot = 0:NSlots-1
                % Update the carrier slot numbers for new slot
                carrier.NSlot = nSlot;

                % Calculate the transport block sizes for the transmission in the slot
                [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
                trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschExtension.TargetCodeRate,pdschExtension.XOverhead);

                % HARQ processing
                for cwIdx = 1:pdsch.NumCodewords
                    % If new data for current process and codeword then create a new DL-SCH transport block
                    if harqEntity.NewData(cwIdx)
                        trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                        setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                        if 0 % 2023-10-04, removed
                            % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
                            if harqEntity.SequenceTimeout(cwIdx)
                                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
                            end
                        end
                    end
                end

                % Encode the DL-SCH transport blocks
                codeTrBlocks = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers, ...
                    pdschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);

                % Create resource grid for a slot
                pdschGrid = nrResourceGrid(carrier, sim.NTxAnts);

                % PDSCH modulation and precoding
                pdschSymbols = nrPDSCH(carrier, pdsch, codeTrBlocks);
                pdschGrid(pdschIndices) = pdschSymbols;

                % PDSCH DM-RS precoding and mapping
                dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
                dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
                pdschGrid(dmrsIndices) = dmrsSymbols;

                % PDSCH PT-RS precoding and mapping
                ptrsSymbols = nrPDSCHPTRS(carrier,pdsch);
                ptrsIndices = nrPDSCHPTRSIndices(carrier,pdsch);
                pdschGrid(ptrsIndices) = ptrsSymbols;

                if 0
                    % export reference gird with DM-RS symbols into Cell
                    dmrsSymCell{nSlot+1} = dmrsSymbols;
                    dmrsIndCell{nSlot+1} = dmrsIndices;
                    ptrsSymbolsCell{nSlot+1} = ptrsSymbols;
                    ptrsIndicesCell{nSlot+1} = ptrsIndices;
                end

                % Creat txGrid by pdschGrid
                txGrid(:, nSlot*carrier.SymbolsPerSlot+1: (nSlot+1)*carrier.SymbolsPerSlot) = pdschGrid;
            end

            % OFDM modulation
            txWaveform = nrOFDMModulate(carrier, txGrid, 'SampleRate', sampleRate);

            % export - configuration
            wv_config_cell = {carrier,pdsch,pdschExtension,sim};

            % export to object - txGrid for demodulation
            obj.Summary.txGrid = txGrid;
            obj.Waveform = txWaveform;
            obj.SubcarrierSpacing = carrier.SubcarrierSpacing;
            obj.CarrierBw = carrier.NSizeGrid * obj.SubcarrierSpacing * 1e3 * 12;
            obj.ChannelBw;

            % export - description
            waveformInfo = nrOFDMInfo(carrier,'SampleRate',sampleRate);
            wvInfo = ['pdsch-'+string(carrier.NSizeGrid)+'RBs-'+'scskHz'+string( obj.SubcarrierSpacing)+'-fs'+string(waveformInfo.SampleRate/1e6)+'MHz-'+string(sim.NFrames)+'Frames'];
            wvInfo = strrep(wvInfo,'.','p');

            if isFnum
                plotSpectrum(txWaveform, sampleRate, {fnum,wvInfo,'Spectrum'});
            end
        end
        function obj = puschRx(obj, rxWaveform, sim, fnum)
            if 1 % validate configuration
                if ~exist('rxWaveform','var')||isempty(rxWaveform)
                    error('nrWaveformGen.puschRx, check input of rxWaveform')
                end
                if ~exist('sim','var')||isempty(sim)
                    sim = obj.PUSCH_config.sim;
                end
                if ~exist('fnum','var')||isempty(fnum)
                    isFnum = 0;
                    fnum = [];
                else
                    isFnum = 1;
                end
            end
            if 1 % install
                carrier = obj.PUSCH_config.carrier;
                pusch = obj.PUSCH_config.pusch;
                sim = obj.PUSCH_config.sim;

                sampleRate = obj.SampleRate;

                NSlots = sim.NFrames * carrier.SlotsPerFrame;
                NSlotSymbs = carrier.SymbolsPerSlot;
            end

            if ~isempty(obj.ChannelModel)&&isfield(obj.Summary,'channelModelPathGains')
                channel = obj.ChannelModel;
                pathGains = obj.Summary.channelModelPathGains;
                isPerfectChannelEstimator = 1;
            else
                isPerfectChannelEstimator = 0;
            end
            % Signal Timing Sychronization
            if sim.isSychronization
                dmrsLayersSymCell = cell(1,NSlots);
                dmrsLayersIndCell = cell(1,NSlots);
                dmrsSlotGrid = [];
                for nSlot = 0:NSlots-1
                    indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);

                    carrier.NSlot = nSlot;
                    tmpGrid = nrResourceGrid(carrier,pusch.NumLayers);

                    % Get reference gird with DM-RS index and symbols
                    dmrsLayersSymCell{nSlot+1} = nrPUSCHDMRS(carrier,pusch);
                    dmrsLayersIndCell{nSlot+1} = nrPUSCHDMRSIndices(carrier,pusch);

                    tmpGrid(dmrsLayersIndCell{nSlot+1}) = dmrsLayersSymCell{nSlot+1};
                    dmrsSlotGrid(:,indSymbs) = tmpGrid;
                end

                % offset calculation by dmrs grids
                if isPerfectChannelEstimator
                    % Perfect synchronization. Use information provided by the
                    % channel to find the strongest multipath component
                    pathFilters = getPathFilters(channel);
                    [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
                    offset = offset+1; % 2023-10-06
                else
                    offset = nrTimingEstimate(carrier,rxWaveform(:,:),dmrsSlotGrid, 'SampleRate', sampleRate);
                end
                rxWaveform = rxWaveform(1+offset:end,:);
            else
                error('nrWaveformGen.puschRx, sim.isSychronization should be TRUE')
            end

            if sim.isOFDMDemodulate
                % Perform OFDM demodulation on the received data to recreate the
                % resource grid, including padding in the event that practical5
                % synchronization results in an incomplete slot being demodulated
                rxGrid = nrOFDMDemodulate(carrier, rxWaveform, 'SampleRate', sampleRate);
                [K, L, R] = size(rxGrid);
                if (L < carrier.SymbolsPerSlot)
                    rxGrid = cat(2, rxGrid, zeros(K, carrier.SymbolsPerSlot-L, R));
                end
            else
                error('nrWaveformGen.puschRx, sim.isOFDMDemodulate should be TRUE')
            end

            if sim.isChannelEstimation
                estChannelGridCell = cell(1,NSlots);
                pxschEqCell = cell(1,NSlots);
                csiCell = cell(1,NSlots);
                noiseEstVec = zeros(NSlots,1);

                eqSymbols = [];  % equalized symbols for constellation plot
                txSymbols = [];

                flagEvmEq = 1;
                if flagEvmEq == 1
                    txGrid = obj.Summary.txGrid;
                end

                % Get pusch indices per slot
                [puschIndices, puschIndicesInfo] = nrPUSCHIndices(carrier, pusch);

                % Channel Estimation and Equalization
                for nSlot = 0:NSlots-1
                    indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);
                    rxSlotGrid = rxGrid(:,indSymbs);
                    txSlotGrid = txGrid(:,indSymbs);

                    % Get dmrs symbols and indices
                    dmrsSymbols = dmrsLayersSymCell{nSlot+1};
                    dmrsIndices = dmrsLayersIndCell{nSlot+1};

                    % Channel estimation
                    [estChannelGrid,noiseEst] = nrChannelEstimate(rxSlotGrid,dmrsIndices,dmrsSymbols,"CDMLengths",pusch.DMRS.CDMLengths);

                    % Get PUSCH resource elements from the received grid
                    [puschRx,puschHest] = nrExtractResources(puschIndices,rxSlotGrid,estChannelGrid);
                    [puschTx,~] = nrExtractResources(puschIndices,txSlotGrid,ones(size(estChannelGrid)));

                    % Equalization
                    [puschEq, csi] = nrEqualizeMMSE(puschRx, puschHest, noiseEst);

                    % Store the equalized symbols and output them for all the slots
                    eqSymbols = [eqSymbols; puschEq];
                    txSymbols = [txSymbols, puschTx];

                    % export to cell - estChannelGrid, noiseEst, pdschEq,
                    % csi
                    estChannelGridCell{nSlot+1} = estChannelGrid;
                    pxschEqCell{nSlot+1} = puschEq;
                    csiCell{nSlot+1} = csi;
                    noiseEstVec(nSlot+1) = noiseEst;
                end

                if flagEvmEq
                    evmEq = evm(txSymbols, eqSymbols);
                    obj.Summary.ChannlEstimation_EQ = ['Channel estimation and equalization evm: ',num2str(round(evmEq,2)), '%'];
                    PLOT_Constellation(eqSymbols, obj.Summary.ChannlEstimation_EQ, fnum, pusch.Modulation, []);
                    fprintf("\n%s\n",obj.Summary.ChannlEstimation_EQ )
                end
            end

            if sim.isCPE*sim.isChannelEstimation
                puschEqCell = cell(1,NSlots);
                csiCell = cell(1,NSlots);
                eqSymbolsCpe = [];
                flagEvmCpe = 1;

                for nSlot = 0:NSlots-1
                    indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);
                    rxSlotGrid = rxGrid(:,indSymbs);

                    carrier.NSlot = nSlot;
                    % Initialize tmporary grid to store equalized symbols
                    tempGrid = nrResourceGrid(carrier, pusch.NumLayers);

                    % Get PT-RS indices and symbols per slot
                    ptrsIndices = nrPUSCHPTRSIndices(carrier, pusch);
                    ptrsSymbols = nrPUSCHPTRS(carrier, pusch);

                    % Extract PT-RS symbols from received grid and estimated channel grid
                    [ptrsRx,ptrsHest,~,~,ptrsHestIndices,ptrsLayerIndices] = nrExtractResources(ptrsIndices,...
                        rxSlotGrid,estChannelGridCell{nSlot+1},tempGrid);

                    % Equalize PT-RS symbols and map them to tempGrid
                    ptrsEq = nrEqualizeMMSE(ptrsRx, ptrsHest, noiseEstVec(nSlot+1));
                    tempGrid(ptrsLayerIndices) = ptrsEq;

                    % Estimate the residual channel at the PT-RS locations in tempGrid
                    cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols);

                    % Sum estimates across subcarriers, receive antennas, and layers. Then,
                    % get the CPR by taking the angle of the resultant sum
                    cpe = angle(sum(cpe, [1 3 4]));

                    % Get PDSCH resource elements from the received grid
                    [puschRx,puschHest] = nrExtractResources(puschIndices,rxSlotGrid,estChannelGridCell{nSlot+1});

                    % CPE Equalization
                    [pdschCPEEq, csi] = nrEqualizeMMSE(puschRx,puschHest,noiseEstVec(nSlot+1));

                    % Map the equalized PDSCH symbols to tempGrid
                    tempGrid(puschIndices) = pdschCPEEq;

                    % Correct CPE in each OFDM symbol within the range of reference PT-RS OFDM symbols
                    if numel(puschIndicesInfo.PTRSSymbolSet) > 0
                        symLoc = puschIndicesInfo.PTRSSymbolSet(1) + 1:puschIndicesInfo.PTRSSymbolSet(end)+1;
                        tempGrid(:, symLoc, :) = tempGrid(:, symLoc, :).*exp(-1i*cpe(symLoc));
                    end

                    % Extract PDSCH symbols from CPE_compensation
                    puschEq_compensation = tempGrid(puschIndices);

                    % Store the equalized symbols and output them for all the slots
                    eqSymbolsCpe = [eqSymbolsCpe; puschEq_compensation]; %#ok<AGROW>

                    % Save pxschEq to cell
                    puschEqCell{nSlot+1} = puschEq_compensation;
                    csiCell{nSlot+1} = csi;
                end

                if flagEvmCpe
                    evmCpe = evm(txSymbols, eqSymbolsCpe);
                    obj.Summary.CPE_EQ = ['Common phase error equalization evm: ',num2str(round(evmCpe,2)), '%'];
                    PLOT_Constellation(eqSymbolsCpe, obj.Summary.CPE_EQ, fnum*10, pusch.Modulation, []);
                end
            end

            if sim.isDecode*sim.isChannelEstimation*sim.isCPE % Decode to throughput
                csi = [];
                simThroughput = 0;
                maxThroughput = 0;

                % Get TargetCodeRat, XOverhead from puschExtension
                puschExtension = obj.PUSCH_config.puschExtension;

                % Set up redundancy version (RV) sequence for all HARQ processes
                if puschExtension.EnableHARQ
                    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
                    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
                    % should be used. If self-decodability is the priority, it should be
                    % taken into account that the upper limit of the code rate at which
                    % each RV is self-decodable is in the following order: 0>3>2>1
                    rvSeq = [0 2 3 1];
                else
                    % HARQ disabled - single transmission with RV=0, no retransmissions
                    rvSeq = 0;
                end

                % Specify the fixed order in which we cycle through the HARQ process IDs
                harqSequence = 0:puschExtension.NHARQProcesses-1;

                % Initialize the state of all HARQ processes
                harqEntity = HARQEntity(harqSequence,rvSeq);

                % Create DL-SCH decoder system object to perform transport channel decoding
                % Use layered belief propagation for LDPC decoding, with half the number of
                % iterations as compared to the default for belief propagation decoding
                decodeULSCH = nrDLSCHDecoder;
                decodeULSCH.MultipleHARQProcesses = true;
                decodeULSCH.TargetCodeRate = puschExtension.TargetCodeRate;
                decodeULSCH.LDPCDecodingAlgorithm = puschExtension.LDPCDecodingAlgorithm;
                decodeULSCH.MaximumLDPCIterationCount = puschExtension.MaximumLDPCIterationCount;

                for nSlot = 0:NSlots-1
                    carrier.NSlot = nSlot;

                    % Decode the pdsch physical channel
                    [rxBits, rxSymbols] = nrPUSCHDecode(carrier, pusch, puschEqCell{nSlot+1}, noiseEstVec(nSlot+1));

                    if strcmpi(obj.wvFormat,'pusch')
                        % Apply channel state information (CSI) produced by the equalizer,
                        % including the effect of transform precoding if enabled
                        if (pusch.TransformPrecoding)
                            MRB  = numel(pusch.PRBSet);
                            MSC = MRB * 12;
                            csi = nrTransformDeprecode(csiCell{nSlot+1},MRB) / sqrt(MSC);
                            csi = repmat(csi((1:MSC:end).'),1,MSC).';
                            csi = reshape(csi,size(rxSymbols));
                        end
                        csi = nrLayerDemap(csiCell{nSlot+1});
                        Qm = length(rxBits) / length(rxSymbols);
                        csi = reshape(repmat(csi{1}.',Qm,1),[],1);
                        rxBits = rxBits .* csi;
                    else
                        % Scale LLRs by CSI
                        csi = nrLayerDemap(csiCell{nSlot+1}); % CSI layer demapping
                        for cwIdx = 1:pdsch.NumCodewords
                            Qm = length(rxBits{cwIdx}) / length(rxSymbols{cwIdx}); % bits per symbol
                            csi{cwIdx} = repmat(csi{cwIdx}.', Qm, 1); % expand by each bit per symbol
                            rxBits{cwIdx} = rxBits{cwIdx} .* csi{cwIdx}(:); % scale by CSI
                        end
                    end

                    % Get transport block size
                    trBlkSizes = nrTBS(pusch.Modulation,pusch.NumLayers,numel(pusch.PRBSet),puschIndicesInfo.NREPerPRB,puschExtension.TargetCodeRate,puschExtension.XOverhead);
                    % Decode the UL-SCH transport channel
                    decodeULSCH.TransportBlockLength = trBlkSizes;
                    [decbits, blkerr] = decodeULSCH(rxBits, pusch.Modulation, pusch.NumLayers, ...
                        harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

                    % Store value to calculate throuhput
                    if ~blkerr
                        isPass = 'transmission passed';
                    else
                        isPass = 'transmission failed';
                    end

                    trBlkSize = numel(decbits);
                    simThroughput = simThroughput + (~blkerr * trBlkSize);
                    maxThroughput = maxThroughput + trBlkSize;
                    fprintf('\n(%3.2f%%) NSlot=%d, %s',100*(nSlot+1)/NSlots,nSlot,isPass);
                end

                fprintf('\n\n')
                duration_ms = sim.NFrames*10;
                obj.Summary.Throughput = ['Throughput = ' + string(simThroughput*1e-6/duration_ms) + 'Mbps (' + string(simThroughput/maxThroughput*100) + '%)'];
            end
        end
        function obj = pdschRx(obj, rxWaveform, sim, fnum)
            if 1 % validate configuration
                if ~exist('rxWaveform','var')||isempty(rxWaveform)
                    error('nrWaveformGen.pdschRx, check input of rxWaveform')
                end
                if ~exist('sim','var')||isempty(sim)
                    sim = obj.PDSCH_config.sim;
                end
                if ~exist('fnum','var')||isempty(fnum)
                    isFnum = 0;
                    fnum = [];
                else
                    isFnum = 1;
                end
            end
            if 1 % install
                carrier = obj.PDSCH_config.carrier;
                pdsch = obj.PDSCH_config.pdsch;
                sim = obj.PDSCH_config.sim;

                sampleRate = obj.SampleRate;

                NSlots = sim.NFrames * carrier.SlotsPerFrame;
                NSlotSymbs = carrier.SymbolsPerSlot;
            end

            if ~isempty(obj.ChannelModel)&&isfield(obj.Summary,'channelModelPathGains')
                channel = obj.ChannelModel;
                pathGains = obj.Summary.channelModelPathGains;
                isPerfectChannelEstimator = 1;
            else
                isPerfectChannelEstimator = 0;
            end
            % Signal Timing Sychronization
            if sim.isSychronization
                dmrsSymCell = cell(1,NSlots);
                dmrsIndCell = cell(1,NSlots);
                dmrsSlotGrid = [];
                for nSlot = 0:NSlots-1
                    indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);

                    carrier.NSlot = nSlot;
                    tmpGrid = nrResourceGrid(carrier,pdsch.NumLayers);

                    % Get reference gird with DM-RS index and symbols
                    dmrsSymCell{nSlot+1} = nrPDSCHDMRS(carrier,pdsch);
                    dmrsIndCell{nSlot+1} = nrPDSCHDMRSIndices(carrier,pdsch);

                    tmpGrid(dmrsIndCell{nSlot+1}) = dmrsSymCell{nSlot+1};
                    dmrsSlotGrid(:,indSymbs) = tmpGrid;
                end

                % offset calculation by dmrs grids
                if isPerfectChannelEstimator
                    % Perfect synchronization. Use information provided by the
                    % channel to find the strongest multipath component
                    pathFilters = getPathFilters(channel);
                    [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
                    offset = offset+1; % 2023-10-06
                else
                    offset = nrTimingEstimate(carrier,rxWaveform(:,:),dmrsSlotGrid, 'SampleRate', sampleRate);
                end
                rxWaveform = rxWaveform(1+offset:end,:);
            else
                error('nrWaveformGen.pdschRx, sim.isSychronization should be TRUE')
            end

            if sim.isOFDMDemodulate
                % Perform OFDM demodulation on the received data to recreate the
                % resource grid, including padding in the event that practical5
                % synchronization results in an incomplete slot being demodulated
                rxGrid = nrOFDMDemodulate(carrier, rxWaveform, 'SampleRate', sampleRate);
                [K, L, R] = size(rxGrid);
                if (L < carrier.SymbolsPerSlot)
                    rxGrid = cat(2, rxGrid, zeros(K, carrier.SymbolsPerSlot-L, R));
                end
            else
                error('nrWaveformGen.pdschRx, sim.isOFDMDemodulate should be TRUE')
            end

            if sim.isChannelEstimation
                estChannelGridCell = cell(1,NSlots);
                pxschEqCell = cell(1,NSlots);
                csiCell = cell(1,NSlots);
                noiseEstVec = zeros(NSlots,1);

                eqSymbols = [];  % equalized symbols for constellation plot
                txSymbols = [];

                flagEvmEq = 1;
                if flagEvmEq == 1
                    txGrid = obj.Summary.txGrid;
                end

                % Get pdsch indices per slot
                [pdschIndices, pdschIndicesInfo] = nrPDSCHIndices(carrier, pdsch);

                % Channel Estimation and Equalization
                for nSlot = 0:NSlots-1
                    indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);
                    rxSlotGrid = rxGrid(:,indSymbs);
                    txSlotGrid = txGrid(:,indSymbs);

                    % Get dmrs symbols and indices
                    dmrsSymbols = dmrsSymCell{nSlot+1};
                    dmrsIndices = dmrsIndCell{nSlot+1};

                    % Channel estimation
                    [estChannelGrid,noiseEst] = nrChannelEstimate(rxSlotGrid,dmrsIndices,dmrsSymbols,"CDMLengths",pdsch.DMRS.CDMLengths);

                    % Get PDSCH resource elements from the received grid
                    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxSlotGrid,estChannelGrid);
                    [pdschTx,~] = nrExtractResources(pdschIndices,txSlotGrid,ones(size(estChannelGrid)));

                    % Equalization
                    [pdschEq, csi] = nrEqualizeMMSE(pdschRx, pdschHest, noiseEst);

                    % Store the equalized symbols and output them for all the slots
                    eqSymbols = [eqSymbols; pdschEq];
                    txSymbols = [txSymbols, pdschTx];

                    % export to cell - estChannelGrid, noiseEst, pdschEq,
                    % csi
                    estChannelGridCell{nSlot+1} = estChannelGrid;
                    pxschEqCell{nSlot+1} = pdschEq;
                    csiCell{nSlot+1} = csi;
                    noiseEstVec(nSlot+1) = noiseEst;
                end

                if flagEvmEq
                    evmEq = evm(txSymbols, eqSymbols);
                    obj.Summary.ChannlEstimation_EQ = ['Channel estimation and equalization evm: ',num2str(round(evmEq,2)), '%'];
                    PLOT_Constellation(eqSymbols, obj.Summary.ChannlEstimation_EQ, fnum, pdsch.Modulation, []);
                    fprintf("\n%s\n",obj.Summary.ChannlEstimation_EQ )
                end
            end

            if sim.isCPE*sim.isChannelEstimation
                pdschEqCell = cell(1,NSlots);
                csiCell = cell(1,NSlots);
                eqSymbolsCpe = [];
                flagEvmCpe = 1;

                if pdsch.EnablePTRS
                    for nSlot = 0:NSlots-1
                        indSymbs = nSlot*NSlotSymbs+(1:NSlotSymbs);
                        rxSlotGrid = rxGrid(:,indSymbs);

                        carrier.NSlot = nSlot;
                        % Initialize tmporary grid to store equalized symbols
                        tempGrid = nrResourceGrid(carrier, pdsch.NumLayers);

                        % Get PT-RS indices and symbols per slot
                        ptrsIndices = nrPDSCHPTRSIndices(carrier, pdsch);
                        ptrsSymbols = nrPDSCHPTRS(carrier, pdsch);

                        % Extract PT-RS symbols from received grid and estimated channel grid
                        [ptrsRx,ptrsHest,~,~,ptrsHestIndices,ptrsLayerIndices] = nrExtractResources(ptrsIndices,...
                            rxSlotGrid,estChannelGridCell{nSlot+1},tempGrid);

                        % Equalize PT-RS symbols and map them to tempGrid
                        ptrsEq = nrEqualizeMMSE(ptrsRx, ptrsHest, noiseEstVec(nSlot+1));
                        tempGrid(ptrsLayerIndices) = ptrsEq;

                        % Estimate the residual channel at the PT-RS locations in tempGrid
                        cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols);

                        % Sum estimates across subcarriers, receive antennas, and layers. Then,
                        % get the CPR by taking the angle of the resultant sum
                        cpe = angle(sum(cpe, [1 3 4]));

                        % Get PDSCH resource elements from the received grid
                        [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxSlotGrid,estChannelGridCell{nSlot+1});

                        % CPE Equalization
                        [pdschCPEEq, csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEstVec(nSlot+1));

                        % Map the equalized PDSCH symbols to tempGrid
                        tempGrid(pdschIndices) = pdschCPEEq;

                        % Correct CPE in each OFDM symbol within the range of reference PT-RS OFDM symbols
                        if numel(pdschIndicesInfo.PTRSSymbolSet) > 0
                            symLoc = pdschIndicesInfo.PTRSSymbolSet(1) + 1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
                            tempGrid(:, symLoc, :) = tempGrid(:, symLoc, :).*exp(-1i*cpe(symLoc));
                        end

                        % Extract PDSCH symbols from CPE_compensation
                        pdschEq_compensation = tempGrid(pdschIndices);

                        % Store the equalized symbols and output them for all the slots
                        eqSymbolsCpe = [eqSymbolsCpe; pdschEq_compensation]; %#ok<AGROW>

                        % Save pxschEq to cell
                        pdschEqCell{nSlot+1} = pdschEq_compensation;
                        csiCell{nSlot+1} = csi;
                    end
                end

                if flagEvmCpe
                    evmCpe = evm(txSymbols, eqSymbolsCpe);
                    obj.Summary.CPE_EQ = ['Common phase error equalization evm: ',num2str(round(evmCpe,2)), '%'];
                    PLOT_Constellation(eqSymbolsCpe, obj.Summary.CPE_EQ, fnum*10, pdsch.Modulation, []);
                end
            end

            if sim.isDecode*sim.isChannelEstimation*sim.isCPE % Decode to throughput
                ecsi = [];
                simThroughput = 0;
                maxThroughput = 0;

                % Get TargetCodeRat, XOverhead from pdschExtension
                pdschExtension = obj.PDSCH_config.pdschExtension;

                % Set up redundancy version (RV) sequence for all HARQ processes
                if pdschExtension.EnableHARQ
                    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
                    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
                    % should be used. If self-decodability is the priority, it should be
                    % taken into account that the upper limit of the code rate at which
                    % each RV is self-decodable is in the following order: 0>3>2>1
                    rvSeq = [0 2 3 1];
                else
                    % HARQ disabled - single transmission with RV=0, no retransmissions
                    rvSeq = 0;
                end

                % Specify the fixed order in which we cycle through the HARQ process IDs
                harqSequence = 0:pdschExtension.NHARQProcesses-1;

                % Initialize the state of all HARQ processes
                harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);

                % Create DL-SCH decoder system object to perform transport channel decoding
                % Use layered belief propagation for LDPC decoding, with half the number of
                % iterations as compared to the default for belief propagation decoding
                decodeDLSCH = nrDLSCHDecoder;
                decodeDLSCH.MultipleHARQProcesses = true;
                decodeDLSCH.TargetCodeRate = pdschExtension.TargetCodeRate;
                decodeDLSCH.LDPCDecodingAlgorithm = pdschExtension.LDPCDecodingAlgorithm;
                decodeDLSCH.MaximumLDPCIterationCount = pdschExtension.MaximumLDPCIterationCount;

                for nSlot = 0:NSlots-1
                    carrier.NSlot = nSlot;

                    % Decode the pdsch physical channel
                    [rxBits, rxSymbols] = nrPDSCHDecode(carrier, pdsch, pdschEqCell{nSlot+1}, noiseEstVec(nSlot+1));

                    % Scale LLRs by CSI
                    csi = nrLayerDemap(csiCell{nSlot+1}); % CSI layer demapping
                    for cwIdx = 1:pdsch.NumCodewords
                        Qm = length(rxBits{cwIdx}) / length(rxSymbols{cwIdx}); % bits per symbol
                        csi{cwIdx} = repmat(csi{cwIdx}.', Qm, 1); % expand by each bit per symbol
                        rxBits{cwIdx} = rxBits{cwIdx} .* csi{cwIdx}(:); % scale by CSI
                    end

                    % Get transport block size
                    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschExtension.TargetCodeRate,pdschExtension.XOverhead);
                    % Decode the DL-SCH transport channel
                    decodeDLSCH.TransportBlockLength = trBlkSizes;
                    [decbits, blkerr] = decodeDLSCH(rxBits, pdsch.Modulation, pdsch.NumLayers, ...
                        harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

                    % Store value to calculate throuhput
                    if ~blkerr
                        isPass = 'transmission passed';
                    else
                        isPass = 'transmission failed';
                    end

                    trBlkSize = numel(decbits);
                    simThroughput = simThroughput + (~blkerr * trBlkSize);
                    maxThroughput = maxThroughput + trBlkSize;
                    fprintf('\n(%3.2f%%) NSlot=%d, %s',100*(nSlot+1)/NSlots,nSlot,isPass);
                end

                fprintf('\n\n')
                duration_ms = sim.NFrames*10;
                obj.Summary.Throughput = ['Throughput = ' + string(simThroughput*1e-6/duration_ms) + 'Mbps (' + string(simThroughput/maxThroughput*100) + '%)'];
            end
        end
    end

    methods (Access=public)
        function [tmWaveform, tmRef_class, wvInfo] = TestModelTx(obj,nrRef,bw,scs,dm,ncellid,sampleRate,fnum)
            %% 2023-09-21, refer to https://www.mathworks.com/help/5g/ug/5g-nr-tm-waveform-generation.html

            if 1 % validate configuration
                if ~exist('nrRef','var')||isempty(nrRef)
                    nrRef = obj.TM_config.nrRef;
                end
                if ~exist('bw','var')||isempty(bw)
                    bw = obj.TM_config.bw;
                elseif isnumeric(bw)
                    bw = string(bw)+"MHz";
                end
                if ~exist('scs','var')||isempty(scs)
                    scs = obj.TM_config.scs;
                elseif isnumeric(scs)
                    scs = string(scs)+"kHz";
                end
                if ~exist('dm','var')||isempty(dm)
                    dm = obj.TM_config.dm;
                end
                if ~exist('ncellid','var')||isempty(ncellid)
                    ncellid = obj.TM_config.ncellid;
                end
                sv            = "16.7.0"; % TS 38.141-x version (NR-TM only)
                if ~exist('sampleRate','var')||isempty(sampleRate)
                    sampleRate = obj.SampleRate;
                end
                if ~exist('fnum','var')||isempty(fnum)
                    isFnum = 0;
                    fnum = [];
                else
                    isFnum = 1;
                end

                nSubframes    = obj.NSubframes;
                carrierFreq   = obj.CarrierFreq;
            end

            if 1 % Introduction
                %                 PUSCH FRC Waveform Generation
                %                 Refer to Annex A of TS 38.104, each FRC waveform is defined by a combination of thest parameters:
                %                 * X: Frequency range number
                %                 * Y: Order of the modulation and coding scheme (MCS) used
                %                 * Z: Index of the FRC for the given MCS
                %                 The selected reference channel (rc) must follow a G-FR X - A Y - Z format.
                %                 Refer to Annex of TS 38.104, each  PUSCH FRC waveform defines a number of key parameters including:
                %                 * Frequency range
                %                 * Channel bandwidth
                %                 * Subcarrier spacing
                %                 * Code rate
                %                 * Modulation
                %                 * DM-RS configuration
                %                 Additionally, the asscoiated receiver tests introduce some additional parameters that are defined in:
                %                 * Table 8.2.1.1-1 (Conducted performance requirements for PUSCH without transform precoding)
                %                 * Table 8.2.2.1-1 (Conducted performance requirements for PUSCH with transform precoding)
                %                 * Table 11.2.2.1.1-1 (Radiated performance requirements for BS type 2-O for PUSCH without transform precoding)
                %                 * Table 11.2.2.2.1-1 (Radiated performance requirements for BS type 2-O for PUSCH with transform precoding)
                %                 FR2 waveforms are TDD and 20ms in length, and FR1 waveforms are FDD and 10ms.
                %                 The PUSCH FRC are defined with type A mapping, type B mapping, or, in some cases, either mapping type. In the latter case, type A mapping is configuraed.
                %                 FR2 waveforms without transform precoding are configured with PT-RS, otherwise PT-RS are off.
                %                 Scarmbling identities are set to 0.
                %                 Power levels for all resource elements are uniform.
                %                 The transport block data source is ITU PN9 with RV = 0 ie.e no retransmissions.
            end

            % Create generator object for the above NR-TM/PDSCH FRC reference model
            try
                tmRef_class = hNRReferenceWaveformGenerator(nrRef,bw,scs,dm,ncellid,sv);
                if tmRef_class.IsReadOnly
                    tmRef_class = makeConfigWritable(tmRef_class);
                end
            catch
                tmRef_class = hNRReferenceWaveformGenerator_kc(nrRef,bw,scs,dm,ncellid,sv);
            end

            % Assign carrier settings
            if ~isempty(nSubframes)
                tmRef_class.Config.NumSubframes        = nSubframes;
            else
                nSubframes  = tmRef_class.Config.NumSubframes;
            end
            if ~isempty(sampleRate)
                tmRef_class.Config.SampleRate          = sampleRate;
            end
            if ~isempty(carrierFreq)
                tmRef_class.Config.CarrierFrequency    = carrierFreq;
            else
                carrierFreq = tmRef_class.Config.CarrierFrequency;
            end

            % Generate waveform
            [tmWaveform,tmWaveInfo,resourceInfo] = generateWaveform(tmRef_class);

            if isempty(sampleRate)
                sampleRate = length(tmWaveform) / (nSubframes*1e-3);
            end

            % export - NSubframes and SampleRate
            obj.NSubframes = nSubframes;
            obj.SampleRate = sampleRate;
            obj.SubcarrierSpacing = tmRef_class.Config.SCSCarriers{:}.SubcarrierSpacing;
            obj.CarrierBw = tmRef_class.Config.SCSCarriers{:}.NSizeGrid * obj.SubcarrierSpacing * 1e3 * 12;
            obj.ChannelBw = tmRef_class.Config.ChannelBandwidth * 1e6;
            obj.TM_config.bw = string(tmRef_class.Config.ChannelBandwidth) + "MHz";
            obj.TM_config.scs = string(tmRef_class.Config.SCSCarriers{:}.SubcarrierSpacing) + "kHz";
            if ~isempty(tmRef_class.ConfiguredModel{4})
                obj.TM_config.dm = tmRef_class.ConfiguredModel{4};
            else
                obj.TM_config.dm = "FDD";
            end

            % export - waveform information
            duration_ms = nSubframes;
            wvInfo = string(nrRef) + "-" + string(obj.TM_config.bw) + "-" + string(obj.TM_config.scs) + "-" + string(obj.TM_config.dm) + "-" + string(duration_ms) + "ms";
            wvInfo = strrep(wvInfo,'.','p');

            % export - class
            obj.TM_config.class = tmRef_class;

        end

        function [evmInfo,eqSym,refSym,throughputResults] = TestModelRx(obj, waveform, tmCarrierConfig, cfg, fnum)
            if ~exist('fnum','var')||isempty(fnum)
                fnum = [];
                isFnum = 0;
            else
                isFnum = 1;
            end

            if ~exist("tmCarrierConfig",'var')||isempty(tmCarrierConfig)||~strcmpi(class(tmCarrierConfig),'nrULCarrierConfig')||~strcmpi(class(tmCarrierConfig),'nrDLCarrierConfig')
                tmCarrierConfig = obj.TM_config.class.Config;
            end

            if ~exist('cfg','var')||isempty(cfg)
                % Compute and display EVM measurements
                cfg = struct();
                cfg.Evm3Gpp = true;
                cfg.PlotEVM = isFnum;
                cfg.DisplayEVM = true;
            end

            [evmInfo,eqSym,refSym,throughputResults] = nrWaveformDemod_hNRPUSCHEVM(tmCarrierConfig,waveform,cfg);
            fprintf("\nemv rms %0.2f, evm peak %0.2f, %s \n", round(evmInfo.OverallEVM.RMS*100,2), round(evmInfo.OverallEVM.Peak*100,2), throughputResults.results)
            obj.Summary.Throughput = throughputResults.results;
            obj.Summary.evmRMS = round(evmInfo.OverallEVM.RMS*100,2);
            obj.Summary.evmPeak = round(evmInfo.OverallEVM.Peak*100,2);
        end
    end
    % end-methods

    methods(Access=public)
        function [rxWaveform, pnoise, pnInfo] = setPhaseNoise(obj, txWaveform, pnModel, sampleRate)
            disp('Phase Noise Model Application ...')
            if ~exist('sampleRate','var')||isempty(sampleRate)
                sampleRate = obj.SampleRate;
            else

            end
            if strcmpi(class(pnModel),'comm.PhaseNoise')
                pnoise = pnModel;
                pnInfo = 'phaseNoise';
            elseif ~isempty(sampleRate)&&isfield(pnModel,'FcHz')&&isfield(pnModel,'PNModel')
                % Phase noise level
                foffsetLog = (4:0.1:log10(sampleRate/2)); % Model offset from 1e4 Hz to sr/2 Hz
                foffset = 10.^foffsetLog;         % Linear frequency offset
                pn_PSD = hPhaseNoisePoleZeroModel(foffset,pnModel.Fc,pnModel.PNModel); % dBc/Hz
                if 1
                    % Visualize spectrum mask of phase noise
                    figure
                    semilogx(foffset,pn_PSD)
                    xlabel('Frequency offset (Hz)')
                    ylabel('dBc/Hz')
                    title('Phase noise magnitude response')
                    grid on
                end
                % Set phase noise level
                pnoise = comm.PhaseNoise('FrequencyOffset',foffset,'Level',pn_PSD,'SampleRate',sampleRate);
                pnoise.RandomStream = "mt19937ar with seed";

                switch pnModel.PNModel
                    case 'A'
                        pnInfo = 'TDocR1-163984SetA';
                    case 'B'
                        pnInfo = 'TDocR1-163984SetB';
                    case 'C'
                        pnInfo = 'TR38p803';
                end
            else
                error('nrWaveformGen.setPhaseNoise, check the input of pnModel')
            end
            rxWaveform = pnoise(txWaveform);
        end

        function [rxWaveform,pathGains,sampleTimes,maxChDelay,channelModelInfo] = setPropagationChannel(obj, txWaveform, channel, sampleRate)
            disp('Propagation Channel Model Construction ...')
            try
                channel.SampleRate = sampleRate;
            end
            if ~(strcmpi(class(channel),'nrCDLChannel')||strcmpi(class(channel),'nrTDLChannel'))
                error('nrWaveformGen.setPropagationChannel, check the input of channel')
            end
            % Get the maximum channel delay
            chInfo = info(channel);
            maxChDelay = ceil(max(chInfo.PathDelays * channel.SampleRate)) + chInfo.ChannelFilterDelay;
            txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];
            [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
            channelModelInfo = channel.DelayProfile;
        end
    end % end-method

end % end-classdef
