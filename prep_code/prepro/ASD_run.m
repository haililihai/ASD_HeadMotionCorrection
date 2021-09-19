% Not fully completed: usm, um, pitt2, kki
% sjh?
wd=fullfile(pwd(),'../..');
site={'TR2_200','TR2_240','TR3_120','TR3_240'};
cores=6;
parpool('local',cores);
for i=1:numel(site)
    sub_list=fullfile(wd,site{i},'sub.csv');
    subject=readtable(sub_list);
    parfor j=1:numel(subject.ParticipantID)
        if ~exist(fullfile(wd,site{i},'rest',subject.ParticipantID{j},'prepro'),'dir')
            run_prepro(site{i},'Sess1_Scan1',subject.ParticipantID{j});
            run_prepro(site{i},'Sess1_Scan2',subject.ParticipantID{j});
        end
    end
end
