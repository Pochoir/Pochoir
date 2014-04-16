#declare -a arr=("element1" "element2" "element3")
declare -a CASES=("heat_1D_NP" "heat_2D_NP" "heat_2D_P" "3dfd" "heat_3D_NP" "heat_4D_NP" "life" "rna" "psa_struct" "lcs" "apop" "LBM/BIN/lbm_tang")
## now loop through the above array
for i in "${CASES[@]}"
do
   echo "$i"
   # or do whatever with individual element of the array
   mv $i $1/
done
