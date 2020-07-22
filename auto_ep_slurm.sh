#!/bin/sh

file_location="data/PLOT3/data/"
kill_radius=(0 3)
inhib_radius=(0 4)
grow_radius=(0 5)
mutation=(0.0 0.1 2)
kill_margin=(0.0 0.3 2)
inhib_margin=(0.0 0.3 2)
version=3

num_jobs=4

# will create N 'emabrasingly parallel' jobs. Must be a power of 2.

split_int_array () {
    strt=${1[0]}
    end=${1[1]}
    sze=$(( $end - $start ))
    mid=$(( $strt + $sze / 2 ))
    result=(($strt $mid) ($(( $mid + 1)) $end))
    return result
}

split_float_array () {
    strt=${1[0]}
    end=${1[1]}
    cnt=${1[2]}
    #rt=$(( ($strt - $end) / ($cnt - 1) ))
    rt=eval "bc <<< 'scale=3; ($end - $strt)/($cnt-1)'"
    c_mid=$(( cnt / 2 ))
    #h1=$(( $strt + $rt * $c_mid))
    h1=eval "bc <<< 'scale=3; $strt + $rt * $c_mid'"
    h2=eval "bc <<< 'scale=3; $h1 + $rt'"
    #h2=$(( $h1 + $rt ))
    result=(($str $h1 $c_mid) ($h2 $end $(( $cnt - $c_mid ))))
}

jobs=(($kill_radius $inhib_radius $grow_radius $mutation $kill_margin $inhib_margin $version))

split_point=0
while [ $num_jobs -gt 1 ]; do

    new_jobs=()
    for j in ${jobs[@]}; do
        if [ ${#(j[${split_point}])} -eq 3 ]
            js = split_float_array ${j[$split_point]}
            for s in ${js[@]}; do
                new_jobs+=((${s[@]:0:$(( $split_point - 1))}, $s, ${s[@]:$((split_point + 1)):$((${#j[@]} - $split_point))})
            done
        fi
    done
    jobs=$new_jobs

    num_jobs=$(( $num_jobs / 2 ))
    split_point=$(( $split_point + 1 ))
done

echo $new_jobs


