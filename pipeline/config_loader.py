import yaml
from pathlib import Path

def load_config(config_path="config.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    cfg["general"]["base_dir"] = Path(cfg["general"]["base_dir"])
    cfg["general"]["previous_antibodies_db"] = cfg["general"]["base_dir"] / cfg["general"].get("previous_antibodies_db", "data/All_mAb_20251106_FACS_BLI.xlsx")
    return cfg