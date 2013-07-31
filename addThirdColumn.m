%% Add third column of associated word to AllEventsTimes file
load allEventTimes
%load('E:\DelayWord\Summary\Files\brocaSoundFiles\BrocaWords.mat')
wordlist=names(setdiff(1:length(names),find(strcmp(names,'slide'))));
for i=1:size(allEventTimes,1)
    if ~isempty(find(strcmp(allEventTimes{i,2},wordlist)))
        currentword=allEventTimes{i,2};
        if strmatch(allEventTimes{i-1,2},'slide')
            allEventTimes{i+1,3}=currentword;
        end
    elseif strmatch(allEventTimes{i,2},'slide')
        currentword= allEventTimes{i+1,2};
    end
    allEventTimes{i,3}=currentword;
end
%%
%cd([dtpath filesep expt{e} filesep 'Analog'])
save('allEventTimes','allEventTimes')
