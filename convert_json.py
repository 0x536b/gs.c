import json
import csv

file = r"cameras.json"
csv_file = r"cameras.csv"

with open(file, "r") as f:
    data = json.load(f)

entries = [
    ["id", "width", "height", "pos_x", "pos_y", "pos_z", "rot_00", "rot_01", "rot_02", "rot_10", "rot_11", "rot_12", "rot_20", "rot_21", "rot_22", "focal_x", "focal_y", "img_name"]
]

for i in data:
    entry = []
    entry.append(i["id"])
    entry.append(i["width"])
    entry.append(i["height"])
    entry.append(i["position"][0])
    entry.append(i["position"][1])
    entry.append(i["position"][2])
    entry.append(i["rotation"][0][0])
    entry.append(i["rotation"][0][1])
    entry.append(i["rotation"][0][2])
    entry.append(i["rotation"][1][0])
    entry.append(i["rotation"][1][1])
    entry.append(i["rotation"][1][2])
    entry.append(i["rotation"][2][0])
    entry.append(i["rotation"][2][1])
    entry.append(i["rotation"][2][2])
    entry.append(i["fx"])
    entry.append(i["fy"])
    entry.append(i["img_name"])
    entries.append(entry)

with open(csv_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(entries)