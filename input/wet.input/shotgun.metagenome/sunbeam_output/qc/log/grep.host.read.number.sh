for i in *_1.txt;do base=$(basename ${i} | awk -F_ '{print $1}'); id+="${base}\n" ; reads+="$(cut -f3,4 $i)\n";echo -e ${id}> id; echo -e ${reads}>reads;echo ${i};echo ${base};done
