function [filenames,triggers] = getFilesCriteria(region,contrast,gain_to_look_at,data_dir)
%from all the mat files in data_dir, find the ones that have units in
%region, specified gain to look at and specified contrast value
%returns filenames with sessions fulfilling the criteria, as well as trial
%number(s) of first trials that fulfill condition
files = dir(fullfile(data_dir,'*.mat'));
filenames={};
triggers = {};
for iF = 1:numel(files)
    data = load(fullfile(files(iF).folder,files(iF).name),'anatomy','trial_gain','trial_contrast');
    gaincontrastcombo = false;
    if contains(files(iF).name,'mismatch') || contains(files(iF).name,'playback') || contains(files(iF).name,'dark')
        continue
    end
    if isfield(data,'trial_gain') && isfield(data,'trial_contrast')
        if strcmp(gain_to_look_at,'baseline') || any(gain_to_look_at == 1)
            %find condition with given contrast and no gain change
            bin = data.trial_gain == gain_to_look_at & data.trial_contrast == contrast;
            trigger_tmp = strfind(bin',ones(1,16)) +6;
            trigger_tmp(trigger_tmp <10) = [];
            trigger_tmp(trigger_tmp>(numel(data.trial_gain)-11))=[];
            gaincontrastcombo = nnz(trigger_tmp)>0;
        elseif any(gain_to_look_at == 0)
            %special case: find where gain ==1 and for a given contrast
            %change
            %gaincontrastcombo = any(data.trial_gain == 1 & data.trial_contrast == contrast);
            
            trigger_tmp = strfind((data.trial_gain == 1 & data.trial_contrast == contrast)',[0 1])+1;
%             if any(trigger_tmp)
%                 keyboard
%             end
%                 
            bl_preceeds = data.trial_contrast(trigger_tmp-1)==100;
            trigger_tmp(~bl_preceeds)=[];
            gaincontrastcombo = nnz(trigger_tmp)>0;
        elseif isempty(gain_to_look_at)
            %find contrast steps (100 to contrast)
            pot = strfind(data.trial_contrast'  == contrast,[0 1])+1;
            pot = pot(data.trial_contrast(pot-1)==100);
            gaincontrastcombo = any(pot);
            trigger_tmp = pot;
         
        else
        gaincontrastcombo = any(data.trial_gain == gain_to_look_at & data.trial_contrast == contrast);
        trigger_tmp = strfind((data.trial_gain == gain_to_look_at & data.trial_contrast == contrast)',[0 1])+1;
        end
        end
    containsregion = false;
    if isfield(data,'anatomy')
        if isempty(region)
            containsregion = true; % all files containing gaincontrastcombo
        else
            containsregion = any(startsWith(data.anatomy.cluster_parent,region));
        end
    end
    
    if gaincontrastcombo && containsregion

        filenames{end+1}=fullfile(files(iF).folder,files(iF).name);
        triggers{end+1}=trigger_tmp;
    end
end
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

end

