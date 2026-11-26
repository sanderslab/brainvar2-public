#!/bin/bash
set -euo pipefail

# bgzip & tabix
process_file() {
  FILE=$1
  OUT="${FILE}.gz"

  echo "[INFO] Compressing (keeping original): ${FILE} -> ${OUT}"
  bgzip -@4 -c "${FILE}" > "${OUT}"
  tabix -f -p bed "${OUT}"
}

export -f process_file

# 처리할 파일 리스트
FILES=$(ls /home/ec2-user/data/cluster/leafcutter_perind.counts.nochr.gz.qqnorm_* | grep -E '_([0-9]+)$' || true)

# 병렬로 실행
echo "${FILES}" | parallel -j4 process_file