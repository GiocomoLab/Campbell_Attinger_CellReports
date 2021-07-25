function session_name = read_session_file(session_file,paths)

fid = fopen(fullfile(paths.sessions,strcat(session_file,'.txt')),'r');
session_name = {};
counter = 1;
while ~feof(fid)
    session_name{counter} = fgetl(fid);
    counter = counter+1;
end
session_name = session_name';
fclose(fid);

end