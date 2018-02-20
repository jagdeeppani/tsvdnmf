#set -x
K=$2
K1=$1
n=$3
br=$4
cf=$(echo "8*($n+1)+4"|bc -l)
pg=$(echo "$br*($cf+8)"|bc -l)
L=$(echo "$pg/$cf" | bc -l)
L=${L/\.[0-9]*/}
m=$(echo "$pg*($K*$br-$L)/($L*($br-1))"|bc -l)
b=$(echo "$m * 0.05" | bc -l)
m=${m/\.[0-9]*/}
b=${b/\.[0-9]*/}
awk -f genPara.awk -v n=$n -v pg=$pg -v K=$K1 -v m=$m -v b=$b para.template > sample.para
echo type Generic-Schema ascii 	> sample.scheme
echo comment // 		>> sample.scheme
echo whitespace \' \' >> sample.scheme
i=1
pr=;
while [[ i -le n ]]
	do echo attr a$i double >> sample.scheme; pr="$pr a$i"; i=$((i+1))
	done
echo $n $pr > sample.proj
