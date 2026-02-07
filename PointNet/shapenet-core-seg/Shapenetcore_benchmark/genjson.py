import os
import json

# ---------------- CONFIG ----------------
ROOT_DIR = "."  # directory containing synset folders
SYNSET_FILE = "synsetoffset2category.txt"
# ----------------------------------------

# Read synset mapping
synset_to_class = {}
synset_to_index = {}

with open(SYNSET_FILE, "r") as f:
    for idx, line in enumerate(f):
        class_name, synset = line.strip().split()
        synset_to_class[synset] = class_name
        synset_to_index[synset] = idx

train_split = []
val_split = []
test_split = []

# Iterate over synset folders
for folder_name in sorted(os.listdir(ROOT_DIR)):
    folder_path = os.path.join(ROOT_DIR, folder_name)

    if not os.path.isdir(folder_path):
        continue

    # folder_name is a synset ID like "03001627"
    if folder_name not in synset_to_class:
        continue

    class_name = synset_to_class[folder_name]
    class_index = synset_to_index[folder_name]

    points_dir = os.path.join(folder_path, "points")
    labels_dir = os.path.join(folder_path, "points_label")

    if not (os.path.isdir(points_dir) and os.path.isdir(labels_dir)):
        continue

    # Deterministic order
    files = sorted(f for f in os.listdir(points_dir) if f.endswith(".npy"))
    assert len(files) == 100, f"{folder_name} has {len(files)} files (expected 100)"

    for i, file in enumerate(files):
        base = file[:-4]
        entry = [
            class_index,
            class_name,
            f"{folder_name}/points/{file}",
            f"{folder_name}/points_label/{base}.seg"
        ]

        if i < 70:
            train_split.append(entry)
        elif i < 80:
            val_split.append(entry)
        else:
            test_split.append(entry)

# Sanity checks
assert len(train_split) == 350
assert len(val_split) == 50
assert len(test_split) == 100

# Write JSON files
with open("train_split.json", "w") as f:
    json.dump(train_split, f, indent=2)

with open("val_split.json", "w") as f:
    json.dump(val_split, f, indent=2)

with open("test_split.json", "w") as f:
    json.dump(test_split, f, indent=2)

print("✅ Per-folder deterministic split complete")
print(f"Train: {len(train_split)}")
print(f"Val:   {len(val_split)}")
print(f"Test:  {len(test_split)}")