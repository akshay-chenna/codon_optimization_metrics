while read -r l ; do grep -A 1 AUP $l.metrics | awk 'BEGIN {found=0} /AUP/ {found=1; next} found {print; found=0}' ; done < list.txt  >> AUP.txt
while read -r l ; do grep -A 1 MFE $l.metrics | awk 'BEGIN {found=0} /MFE/ {found=1; next} found {print; found=0}' | head -1 ; done < list.txt  >> MFE.txt
while read -r l ; do grep -A 1 CAI $l.metrics | awk 'BEGIN {found=0} /CAI/ {found=1; next} found {print; found=0}' | head -1 ; done < list.txt  >> CAI.txt
while read -r l ; do tail -1 $l.metrics | awk '{print $1}' | cut -b 2- | cut -d , -f 1 ; done < list.txt  >> degscore.txt
while read -r l ; do tail -1 $l.metrics | awk '{print $2}' | cut -d ")" -f 1 ; done < list.txt >> half_life.txt
while read -r l ; do grep -A 1 "Sequence:" $l.metrics | awk 'BEGIN {found=0} /Sequence/ {found=1; next} found {print; found=0}' >> ${l}.txt ; done < list.txt
while read -r l ; do grep -A 1 MFE $l.metrics | awk 'BEGIN {found=0} /MFE/ {found=1; next} found {print; found=0}' | tail -1 > ${l}.structure ; done < list.txt # Gets the structure
