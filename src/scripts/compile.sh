suffix="pointer iter macro type"
#file="tb_heat_2D_NP tb_heat_2D_P tb_life tb_3dfd psa rna"
file="tb_3dfd"
set -x
for f in $file; do
    ./pochoir $f".cpp"
    mv $f"_pochoir" $f"_zero_pointer" 
    ./pochoir $f".cpp" -split-macro-shadow
    mv $f"_pochoir" $f"_macro"
    ./pochoir $f".cpp" -split-pointer
    mv $f"_pochoir" $f"_type"
    ./pochoir $f".cpp" -split-opt-pointer
    mv $f"_pochoir" $f"_iter"
done
set +x
