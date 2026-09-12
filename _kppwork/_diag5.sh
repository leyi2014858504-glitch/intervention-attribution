#!/bin/bash
echo "=== disk ==="
df -h / /root /mnt/d /tmp 2>&1
echo ""
echo "=== inodes ==="
df -i / /root 2>&1
echo ""
echo "=== /root/kpp6 identity ==="
ls -ld /root/kpp6
ls -la /root/kpp6/cb05.exe
md5sum /root/kpp6/cb05.exe /root/isotest/cb05.exe 2>&1
echo ""
echo "=== run twice in a fresh dir ==="
rm -rf /root/twice
cp -r /root/isotest /root/twice
cd /root/twice
for i in 1 2 3; do
  rm -f traj.dat run.log state_out.dat
  timeout 60 ./cb05.exe > run.log 2>&1
  echo "run$i rc=$? traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null) logbytes=$(wc -c < run.log 2>/dev/null)"
done
echo ""
echo "=== /root usage ==="
du -sh /root 2>/dev/null
ls -la /root | head -20
