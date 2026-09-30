#!/bin/bash
# Grid facts for the MDA8 endpoint-counting question.
F=/root/kpp7/traj.dat
awk 'NR==1{print "HDR: " $0}
     NR==2{print "R1:  " $0; a=$1}
     NR==3{print "R2:  " $0; b=$1}
     {t=$1; n=NR-1}
     END{print "last t = " t
         print "nrows = " n
         printf "dt = %.6f h = %.2f min\n", b-a, (b-a)*60
         printf "8h / dt = %.2f  -> w = %d\n", 8.0/(b-a), int(8.0/(b-a)+0.5)
         printf "span of w samples = %.4f h (95 gaps)\n", (int(8.0/(b-a)+0.5)-1)*(b-a)}' "$F"
