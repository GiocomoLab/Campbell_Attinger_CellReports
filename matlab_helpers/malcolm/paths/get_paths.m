function paths = get_paths()

% change these paths to point to the data, helpers, and intermediate data folders:

paths = struct;
paths.data = '/Volumes/T7/data_to_publish_2'; % path to Neuropixels data
paths.tetrode_data = 'C:\Users\malcg\Dropbox (Personal)\GiocomoLab\Malcolms_VR_data'; % path to tetrode data
paths.matlab_helpers = '/Users/attialex/Downloads\Campbell_Attinger_CellReports\matlab_helpers';
paths.intermediate_data = '/Users/attialex/Downloads/Campbell_Attinger_CellReports/intermediate_data';
paths.spikes_repo = '/Users/attialex/code/spikes'; % download spikes repo from here: https://github.com/cortex-lab/spikes
paths.figs = '/Volumes/T7/figures';
end