for sensitivity in $(seq 2 0.1 3); do 
    echo -ne "$sensitivity\t" 
    scripts/run_vad.sh $sensitivity | fgrep TOTAL 
done | sort -t: -k2n