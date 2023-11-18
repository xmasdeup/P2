for sensitivity in $(seq 4 0.2 8); do 
    echo -ne "$sensitivity\t" 
    scripts/run_vad.sh $sensitivity | fgrep TOTAL 
done | sort -t: -k2n