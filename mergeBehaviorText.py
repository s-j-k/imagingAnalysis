import pandas as pd
from glob import glob
from tifffile import TiffFile
from pathlib import Path

def load_master_bhv_file(file_name):
    header = ["session",
              "trial_number",
              "tone",
              "response",
              "no_lick_period",
              "response_time",
              "delay_after_response",
              "total_trial_time_minus_resp_time",
              "lick_frame",
              "reward_frame",
              "total_trial_time",
              "tone_frame",
              "context"]
    return pd.read_csv(file_name, sep=',', header=None, names=header)

animal_id = 240
day = 4

session_path = Path(fr"O:\sjk\DATA\imagingData\meso\sk{animal_id}\behavior\{animal_id}-D{day}") #this is the data folder w suiet2p folder inside
master_bhv = []
frame_count, num_trial = 0, 0

num_tif_files = len(glob(str(session_path.joinpath("raw", "*.tif"))))

for i_block in range(num_tif_files-1):
    df = load_master_bhv_file(session_path.joinpath(F"sk{animal_id}_{day}v{i_block+1}.txt"))
    for column in ["lick_frame", "reward_frame", "tone_frame"]:
        is_valid = df[column] != 1000000.0
        df.loc[is_valid, column] += num_trial
    df['trial_number'] += num_trial
    master_bhv.append(df)

    num_frame_in_tiff = len(TiffFile(glob(str(session_path.joinpath('raw', "*.tif")))[i_block]).pages)
    frame_count += num_frame_in_tiff // 2
    num_trial += len(df)

master_bhv = pd.concat(master_bhv).reset_index(drop=True)
master_bhv.to_csv(session_path.joinpath("master_behavior_file.txt"), header=False, index=False)

# for i_block, block_path in enumerate(glob(str(session_path.joinpath("block*")))):
#     df = load_master_bhv_file(Path(block_path).joinpath("master_behavior_file.txt"))
#     for column in ["lick_frame", "reward_frame", "tone_frame"]:
#         is_valid = df[column] != 1000000.0
#         df.loc[is_valid, column] += num_trial
#     df['trial_number'] += num_trial
#     master_bhv.append(df)

#     f1 = np.load(Path(block_path).joinpath("suite2p", "plane0", "F.npy"))
#     frame_count += 2 * f1.shape[1]
#     num_trial += len(df)

# master_bhv = pd.concat(master_bhv).reset_index(drop=True)
# master_bhv.to_csv(session_path.joinpath("combined", "master_behavior_file.txt"), header=False, index=False)