#include <cstdio>
#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <cmath>
#include <string.h>

using namespace std;
int main(int argc, char * argv[])
{
	// Current date/time based on current system
	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	char time [100] ;
	strftime(time, 100, "%m%d%Y_%H%M%S", localtm) ;
	cout << "The local date is: " << time << endl;
	char cmd [500], tmp [500], dir[100] ;
	char file [500]; 
	char prefix[100] ;
	char suffix[100] ;

	sprintf(dir, "tc_tune_sub_space_%s", time) ;
	sprintf(cmd, "mkdir %s", dir) ;
	system (cmd) ;

	//cd into the dir
	//sprintf(cmd, "cd %s", dir) ;
	//system (cmd) ;
	int option = 0 ; 
	if (argc > 1) 
	{
		option = atoi(argv[1]) ;
	}
	sprintf(prefix, "CILK_NWORKERS=1 ") ; 
	int NCORES = 1 ;
	//simulate each of the test cases for TRAP
	for (int i = 0 ; i < 1 ; i ++)
	{
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_NP %d %d %d >> %s 2>&1", 1000, 2000, 512, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
          
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_NP %d %d %d >> %s 2>&1", 1000, 2000, 512, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 4000, 4000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 4000, 4000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	/*sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/life %d %d %d >> %s 2>&1", 3000, 2000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	*/
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/3dfd %d %d %d %d >> %s 2>&1", 200, 200, 200, 16, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/3dfd %d %d %d %d >> %s 2>&1", 200, 200, 200, 16, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lbm_tang") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/lbm_tang %d %d %d %d %s %d %d >> %s 2>&1", 64, 100, 100, 130, "result_lbm", 0, 0, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lbm_tang") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/lbm_tang %d %d %d %d %s %d %d >> %s 2>&1", 64, 100, 100, 130, "result_lbm", 0, 0, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/apop") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/apop -s %d -t %d >> %s 2>&1", 2000000, 524288, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/apop") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/apop -s %d -t %d >> %s 2>&1", 2000000, 524288, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
         
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P_10000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 10000, 10000, 4096, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
         
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P_10000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 10000, 10000, 4096, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P2_100_20000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 100, 20000, 8192, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P2_100_20000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_P %d %d %d >> %s 2>&1", 100, 20000, 8192, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/heat_4D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_4D_NP %d %d %d %d %d >> %s 2>&1", 70, 70, 70, 70, 32, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/heat_4D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_4D_NP %d %d %d %d %d >> %s 2>&1", 70, 70, 70, 70, 32, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life_2000_3000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/life %d %d %d >> %s 2>&1", 2000, 3000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life_2000_3000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/life %d %d %d >> %s 2>&1", 2000, 3000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	/*
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life_10000_10000") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/life %d %d %d >> %s 2>&1", 10000, 10000, 4096, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lbm_tang_1000^3") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/lbm_tang %d %d %d %d %s %d %d >> %s 2>&1", 256, 1000, 1000, 1000, "result_lbm", 0, 0, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	*/
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd_1000^3") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/3dfd %d %d %d %d >> %s 2>&1", 1000, 1000, 1000, 64, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd_1000^3") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/3dfd %d %d %d %d >> %s 2>&1", 1000, 1000, 1000, 64, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
        
        sprintf(suffix, "_trap_%d_core", NCORES) ;
        strcpy(file, dir) ;
        strcat(file, "/heat_2D_NP_16000") ;
        strcat(file, suffix) ;
        strcpy(cmd, prefix) ;
        sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_NP %d %d %d >> %s 2>&1", 16000, 17000, 8192, file) ;
        strcat(cmd, tmp) ;
        cout << cmd << endl ;
        system(cmd) ;
            
        sprintf(suffix, "_trap_%d_core", NCORES) ;
        strcpy(file, dir) ;
        strcat(file, "/heat_2D_NP_16000") ;
        strcat(file, suffix) ;
        strcpy(cmd, prefix) ;
        sprintf(tmp, "~/research/git/pochoir/examples/fixed_space/tc_tune_sub_space/trap/heat_2D_NP %d %d %d >> %s 2>&1", 16000, 17000, 8192, file) ;
        strcat(cmd, tmp) ;
        cout << cmd << endl ;
        system(cmd) ;

	}

	sprintf(dir, "tc_tune_sub_space_%s/", time) ;
	//sprintf(cmd, "mv *sawzoid %s", dir) ;
	//system (cmd) ;
	sprintf(cmd, "mv *trap %s", dir) ;
	system (cmd) ;

	return 0 ;
}
