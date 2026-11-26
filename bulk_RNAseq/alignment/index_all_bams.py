#!/usr/bin/env python3

import os
import subprocess

bam_root = "/mnt/s3/FXN/Bulk_RNA/01-aligned"

count = 0
skipped = 0
failed = 0

for root, dirs, files in os.walk(bam_root):
    for file in files:

        if file != "Aligned.sortedByCoord.out.bam":
            continue

        bam = os.path.join(root, file)

        bai1 = bam + ".bai"
        bai2 = bam.replace(".bam", ".bai")

        if os.path.exists(bai1) or os.path.exists(bai2):
            print(f"Already indexed: {bam}")
            skipped += 1
            continue

        print(f"Indexing: {bam}")

        result = subprocess.run(
            ["samtools", "index", bam],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            count += 1
        else:
            failed += 1
            print(result.stderr)

print("\nDone")
print(f"Indexed : {count}")
print(f"Skipped : {skipped}")
print(f"Failed  : {failed}")
