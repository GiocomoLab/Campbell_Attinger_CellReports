function uniq_id = get_cell_id(session,cell_id)
%BUILD_CELL_ID Build cell ID from session name (session) and cell ID
%(cell_id)
% MGC 4/4/2020

if numel(session) ~= numel(cell_id)
    disp('ERROR: numel(session) does not equal numel(cell_id)');
    return;
end

if numel(session) == 1
    mouse_session = get_mouse_session(session);
    uniq_id = sprintf('%s_c%d',mouse_session,cell_id);
elseif numel(session)>1
    uniq_id = cell(numel(session),1);
    for i = 1:numel(session)
        mouse_session = get_mouse_session(session(i));
        uniq_id{i} = sprintf('%s_c%d',mouse_session,cell_id(i));
    end
elseif numel(session)==0
    disp('WARNING: empty input');
    uniq_id = [];
end

end

