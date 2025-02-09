#!/usr/bin/env python3

from pathlib import Path

thesis_mode = False

with open("targets") as f:
    for filename in f:
        filename = filename.strip()
        target = Path("../thesis/") / filename
        s = ""
        with open(target) as target_file:
            for line in target_file:
                line = line.replace(r"\subsubsubsection", r"\paragraph")
                line = line.replace(r"Case 5", r"Case 1")
                line = line.replace(r"Case 6", r"Case 2")

                if line.strip() == r"\ifthesis":
                    thesis_mode = True
                    continue
                elif line.strip() == r"\ifpaper":
                    continue
                elif line.strip() == r"\fi":
                    thesis_mode = False
                    continue
                if thesis_mode is False:
                    s += line
        with open(filename, "w") as out_file:
            print(s.strip(), file=out_file)
