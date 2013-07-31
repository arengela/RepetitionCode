function getAnalogEvents(dtpath, wpath, expt, names) 
    evnt = ECogFindEvents(dtpath,wpath,{expt},names)
    cd([dtpath filesep expt filesep 'Analog'])
    %% convert transcripts to evnts found
    clear tmp
    idx=1;
    for i=1:length(evnt)
        if evnt(i).confidence>.80
            tmp{idx,1}=evnt(i).StartTime
            tmp{idx,2}=evnt(i).name
            tmp{idx,3}=evnt(i).name
            if strcmp(evnt(i).name,{'Electronic_Chime-KevanGC-495939803'})
                tmp{idx,2}='beep';
                tmp{idx,3}=evnt(i).name
            elseif ~strcmp(evnt(i).name,'slide')
                %keyboard
                idx=idx+1;
                tmp{idx,1}=evnt(i).StopTime
                tmp{idx,2}='we'
                tmp{idx,3}=evnt(i).name
            end
            idx=idx+1;
        end
    end
    tmp2=tmp(:,[1,2]);
    E_times=cell2mat(tmp2(:,1))
    trialslog=tmp2(:,2)
    BadTimesConverterGUI3 (E_times,trialslog,sprintf('transcript_AN%d.lab',2))
    %makeCombinedEventFiles({sprintf('transcript_AN%d.lab',2)})
    try
        makeCombinedEventFiles({'transcript_AN2.lab','transcript_AN1.lab'})
    catch
            makeCombinedEventFiles({'transcript_AN2.lab'})
    end
    addThirdColumn
end