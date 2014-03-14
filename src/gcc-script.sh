
#!/bin/bash
echo "Problem Dimention"
read dim
echo "Problem name"
read name
# @ EKA: replaced the below line with your pochour src directory. All files will be replaced with the modified version so please do this in a test copy of pochoir. 
startdirectory="./"
searchterm0="F"
searchterm1="template <typename F, typename BF>"
searchterm2="template <typename F, typename BF>"
searchterm3="template <typename F>"
searchterm4="template <typename BF>"
searchterm5=", BF>" 
searchterm6=", F>" 
searchterm7=", typename F>"
searchterm8=", typename BF>"



replaceterm1=" "
replaceterm2=", "
replaceterm3="> "

# **********************************************************

echo "***************************************************"
echo "* Search and Replace in Files *"
echo "***************************************************"

i=0; 

mkdir -p "$startdirectory/backup"
  for file in $(grep -l -R $searchterm0 $startdirectory)
    do
#      cp $file $startdirectory/backup/$file.bak
      sed -e "s/$searchterm1/$replaceterm1/ig" -e "s/$searchterm2/$replaceterm1/ig" -e "s/$searchterm3/$replaceterm1/ig"  -e "s/$searchterm4/$replaceterm1/ig"   -e "s/$searchterm5/$replaceterm3/ig"  -e "s/$searchterm6/$replaceterm3/ig"  -e "s/$searchterm7/$replaceterm3/ig" -e "s/$searchterm8/$replaceterm3/ig"   $file > tempfile.tmp
      mv tempfile.tmp $file

    #let i++;

#    echo "Modified: " $file
    done

 
echo "
#define N_RANK2 $dim 
#if N_RANK2==8 
typedef void  (*F)(int,int, grid_info<8> const &) ;  // 
typedef void  (*F1)(int,int,int,int,int,int,int,int,int) ;  // 
typedef void  (*BF)(int,int,int,int,int,int,int,int,int) ;
#elif N_RANK2==7 
typedef void  (*F)(int,int,  grid_info<7> const &) ;  // 
typedef void  (*F1)(int,int,int,int,int,int,int,int, ) ;  // 
typedef void  (*BF)(int,int,int,int,int,int,int,int, ) ;
#elif N_RANK2==6 
typedef void  (*F)(int,int,  grid_info<6> const &) ;  // 
typedef void  (*F1)(int,int,int,int,int,int,int) ;  // 
typedef void  (*BF)(int,int,int,int,int,int,int) ;
#elif N_RANK2==5 
typedef void  (*F)(int,int, grid_info<5> const &) ;  // 
typedef void  (*F1)(int,int,int,int,int,int) ;  // 
typedef void  (*BF)(int,int,int,int,int,int) ;
#elif N_RANK2==4  
typedef void  (*F)(int,int, grid_info<4> const &) ;  // 
typedef void  (*F1)(int, int,int,int,int) ;  // 
typedef void  (*BF)(int,int,int,int,int) ;
#elif N_RANK2==3 
typedef void  (*F)(int,int, grid_info<3> const &) ;  // 
typedef void  (*F1)(int,int,int,int) ;  // 
typedef void  (*BF)(int,int,int,int) ;
#elif N_RANK2==2 
typedef void  (*F)(int,int, grid_info<2> const &) ;  // 
typedef void  (*F1)(int,int,int) ;  // 
typedef void  (*BF)(int,int,int) ; 
#else  
typedef void  (*F)(int,int, grid_info<1> const &) ;  // 
typedef void  (*F1)(int,int) ;  // 
typedef void  (*BF)(int,int) ;
#endif 
" > $startdirectory/GCC-Kerry.hpp
number=1
cd $startdirectory
replaceterm="#if N_RANK2==$number \\
template<> \\
struct meta_grid_boundary  <"
FFile='pochoir_walk.hpp'
FFile2='PMain.hs'
cd $startdirectory
searchterm1="struct meta_grid_boundary <"
searchterm0="struct meta_grid_interior <"
    COO="8 7 6 5 4 3 2 1"
    for number in `echo $COO`
        do
           if [ $number -eq 8 ]; then 
                replaceterm1="#if N_RANK2==$number \\
template<> \\
struct meta_grid_interior  <"
           elif [ $number -ne 1 ]; then 
                replaceterm1="#elif N_RANK2==$number \\
template<> \\
struct meta_grid_boundary  <"
            else
                replaceterm1="#else \\
template<> \\
struct meta_grid_boundary  <"
            fi
           sed -e  "0,/$searchterm1/{s/$searchterm1/$replaceterm1/}" $FFile > tempfile.tmp
           mv tempfile.tmp $FFile
      done        

    COO="8 7 6 5 4 3 2 1"
    for number in `echo $COO`
        do
           if [ $number -eq 8 ]; then
              replaceterm0="#if N_RANK2==$number \\
template<> \\
struct meta_grid_interior  <"
           elif [ $number -ne 1 ]; then
                replaceterm0="#elif N_RANK2==$number \\
template<> \\
struct meta_grid_interior  <"
            else
                replaceterm0="#else \\
template<> \\
struct meta_grid_interior  <"
            fi
           sed -e "0,/$searchterm0/{s/$searchterm0/$replaceterm0/}"  $FFile > tempfile.tmp
           mv tempfile.tmp $FFile
      done

number=1
searchterm="template <int N_RANK>"
replaceterm="#endif \n template <int N_RANK>"

COO="1 2"
    for number in `echo $COO`
        do
           if [ $number -eq 1 ]; then
replaceterm="template < int N_RANK>"
            else
replaceterm="#endif \ntemplate <int N_RANK>"
            fi
     sed -e "0,/$searchterm/{s/$searchterm/$replaceterm/}"  $FFile > tempfile.tmp
     mv tempfile.tmp $FFile
      done

searchterm="using namespace std;"
replaceterm="using namespace std; \\
#include \"GCC-Kerry.hpp\" "
sed -e "s/$searchterm/$replaceterm/"  $FFile > tempfile.tmp
mv tempfile.tmp $FFile

searchterm="static inline void set_worker_count"
replaceterm="#endif \\
static inline void set_worker_count"
sed -e "s/$searchterm/$replaceterm/"  $FFile > tempfile.tmp
mv tempfile.tmp $FFile


searchterm="then iccPPFlags"
replaceterm="then iccPPFlags ++ [\"tb_$name.i\"]"
sed -e "s/$searchterm/$replaceterm/"  $FFile2 > tempfile.tmp
mv tempfile.tmp $FFile2


searchterm="then iccPPDebugFlags" 
replaceterm="then iccPPDebugFlags ++ [\"tb_$name.i\"]"
sed -e "s/$searchterm/$replaceterm/"  $FFile2 > tempfile.tmp
mv tempfile.tmp $FFile2


searchterm="initial_grid, F"
replaceterm="initial_grid, F1"
sed -e "s/$searchterm/$replaceterm/"  $FFile > tempfile.tmp
mv tempfile.tmp $FFile


cd ..








echo " *** All Done! *** Modified files:" 
