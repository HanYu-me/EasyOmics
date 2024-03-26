awk '{i=1;while(i <= NF){col[i]=col[i] $i " ";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' inter_result/snpmatrix.txt |sed 's/IID //g' > inter_result/snpmatrix2.txt
