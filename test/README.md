#Test
    python ../src/AllowedFraction.py 0 50 0.01 >frac
    diff frac frac_0_50_0.01 && rm frac

#lj at T=1.3 with coeffcient from doi:10.1135/cccc2009113 
    ##equal to python ../src/AllowedFraction.py 0 50 0.01 -6.33642 9.41046 20.2825 9.28669
    bash lj.sh -3.31774 2.57993 2.9115 0.698 >lj_1.3  
