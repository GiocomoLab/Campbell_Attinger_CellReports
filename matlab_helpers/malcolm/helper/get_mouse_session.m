function ms = get_mouse_session(session_name)
% MOUSE_SESSION Converts session name into [mouse name]_[session name]
% MGC 4/4/2020

if numel(session_name) == 1
    strsplit_this = strsplit(session_name{1},'_');
    ms = sprintf('%s_%s',strsplit_this{1},strsplit_this{2});
elseif numel(session_name)>1
    ms = cell(numel(session_name),1);
    for i = 1:numel(session_name)
        strsplit_this = strsplit(session_name{i},'_');
        ms{i} = sprintf('%s_%s',strsplit_this{1},strsplit_this{2});
    end
elseif numel(session_name)==0
    disp('WARNING: empty input');
    ms = [];
end

end

