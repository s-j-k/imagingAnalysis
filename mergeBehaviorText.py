from pathlib import Path

import pandas as pd
from tifffile import TiffFile

HEADER = [
    "session",
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
    "context",
]

DATA_ROOT = Path(r"O:\rockfish\su\Data\sk318_reg")


def load_master_bhv_file(file_name):
    return pd.read_csv(file_name, sep=",", header=None, names=HEADER)


def normalize_days(days):
    """Return a list whether days is one day, a list, a tuple, or a range."""
    if isinstance(days, int):
        return [days]

    normalized = [int(day) for day in days]
    if not normalized:
        raise ValueError("At least one day must be supplied.")
    return normalized


def behavior_block_number(path, animal_id, day):
    """Return the block number only for an exact sk<animal>_<day>v<block>.txt name."""
    prefix = f"sk{animal_id}_{day}v"
    filename = path.name

    # Windows globbing can make *.txt also return sidecars such as .txtframes.
    # Require the exact behavior-file suffix and an all-numeric block field.
    if not filename.lower().endswith(".txt"):
        return None
    if not filename.lower().startswith(prefix.lower()):
        return None

    block_text = filename[len(prefix):-4]
    if not block_text.isdigit():
        return None
    return int(block_text)


def merge_behavior_text(animal_id, days, data_root=DATA_ROOT):
    """Concatenate the existing behavior blocks for one day or multiple days."""
    days = normalize_days(days)

    folder_name = (
        f"D{days[0]}"
        if len(days) == 1
        else f"D{days[0]}_D{days[-1]}"
    )
    session_path = Path(data_root) / folder_name

    behavior_files = []
    for day in days:
        numbered_files = []
        for path in session_path.joinpath("behavior").glob(
            f"sk{animal_id}_{day}v*"
        ):
            block_number = behavior_block_number(path, animal_id, day)
            if block_number is not None:
                numbered_files.append((block_number, path))

        day_behavior_files = [
            path for _, path in sorted(numbered_files)
        ]
        if not day_behavior_files:
            raise FileNotFoundError(
                f"No behavior files were found for day {day} in "
                f"{session_path / 'behavior'}"
            )
        behavior_files.extend(day_behavior_files)

    # Pair only the TIFFs needed by the behavior files. Any final baseline TIFF
    # is therefore ignored, and day 1 will stop at 1v3 when only 1v1-1v3 exist.
    tif_files = sorted(session_path.joinpath("raw").glob("*.tif"))
    if len(tif_files) < len(behavior_files):
        raise FileNotFoundError(
            f"Found {len(behavior_files)} behavior blocks but only "
            f"{len(tif_files)} TIFF files in {session_path / 'raw'}"
        )
    block_tif_files = tif_files[: len(behavior_files)]

    master_bhv = []
    frame_count = 0
    num_trial = 0

    for behavior_file, tif_path in zip(behavior_files, block_tif_files):
        df = load_master_bhv_file(behavior_file)

        for column in ["lick_frame", "reward_frame", "tone_frame"]:
            is_valid = df[column] != 1000000.0
            df.loc[is_valid, column] += frame_count

        df["trial_number"] += num_trial
        master_bhv.append(df)

        with TiffFile(tif_path) as tif_file:
            frame_count += len(tif_file.pages) // 2
        num_trial += len(df)

    merged = pd.concat(master_bhv, ignore_index=True)
    output_path = session_path / "master_behavior_file.txt"
    merged.to_csv(output_path, header=False, index=False)
    return output_path


if __name__ == "__main__":
    ANIMAL_ID = 318

    # Either form works:
    DAYS = [1, 7]
    # DAYS = 1
    # DAYS = range(1, 4)

    output_path = merge_behavior_text(ANIMAL_ID, DAYS)
    print(output_path)
