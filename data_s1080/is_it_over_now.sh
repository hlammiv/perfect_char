for f in out*;
        do
                echo -n $f" "
                grep "GMES" $f | wc -l
                echo $f | awk -F[a\_] '{printf "%f ",substr($2,3);}'
#               awk '$0~/GMES:/&&NR>5000{s+=$5;t+=1}END{printf "%f\n",(1.0+s/t/3)}' $f
                awk '$0~/GMES:/{s+=$5;t+=1}END{printf "%f\n",(1-s/t/3)}' $f
done
