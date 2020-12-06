
set multiplot

set size 1,0.5
set origin 0,0.5

p \
'evalContributionsSectors.dat' u 1:2:-1 lc var w lp ti 'Contribution to GF', \
                            '' u 1:2:4 w labels offset 0,1 noti

set origin 0,0
p 'evalContributionsSectors.dat' u 1:3:-1 lc var w lp ti 'Energy'
unset multiplot
