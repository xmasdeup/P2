for umbral1 in $(seq 4 0.1 6); do 
    echo -ne "$umbral1\t" 
    scripts/run_vad.sh $umbral1 | fgrep TOTAL 
done | sort -t: -k2n