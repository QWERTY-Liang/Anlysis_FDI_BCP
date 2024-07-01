% Current Source Density analysis - takes the second spatial derivative
% Download the toolbox from here: http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/
% 8/4/2020 Kieran tweaked the CSD function to make it more efficient, and called it CSD_km. 
% It takes 25 sec where CSD takes 141 sec - quite a time saving when scaled up to many subjects

clear

allsubj = {'P01', 'P02', 'P03', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16'};

bigmatfolder = 'bigmats/';
nchan = 128; % number of channels

% load matrices that CSD function needs to run for the specific biosemi 128 cap:
load CSDStuff % has two matrices G and H that are required for the CSD function, and were constructed based on this particular 128-channel cap that we use in UCD (Biosemi)

for s=1:length(allsubj)       
    disp(['Running ' allsubj{s} '...'])
    load([bigmatfolder allsubj{s} '_intp'])
    erp(1:nchan,:,:) = CSD_km(erp(1:nchan,:,:), G, H);
    % Now save a new version:
    save([bigmatfolder allsubj{s} '_intpCSD'],'erp','t','modir','coh','cond','respLR','RT','blocknum','filenames','anapar','artifact','erpr','validrlock','blink','artifact');
end

