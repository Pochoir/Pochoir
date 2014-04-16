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

	sprintf(dir, "tc_tune_sub_space_sub_time_%s", time) ;
	sprintf(cmd, "mkdir %s", dir) ;
	system (cmd) ;

	//cd into the dir
	//sprintf(cmd, "cd %s", dir) ;
	//system (cmd) ;
	

	sprintf(prefix, "CILK_NWORKERS=1 ") ; 
	int NCORES = 1 ;
	//simulate each of the test cases for TRAP
	//for (int i = 0 ; i < 2 ; i ++)
	for (int i = 0 ; i < 1 ; i ++)
	{
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/heat_2D_NP %d %d %d >> %s 2>&1", 1000, 2000, 512, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/heat_2D_NP %d %d %d >> %s 2>&1", 1000, 2000, 512, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_trap_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/heat_2D_P %d %d %d >> %s 2>&1", 4000, 4000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ; 
	strcpy(file, dir) ;
	strcat(file, "/heat_2D_P") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/heat_2D_P %d %d %d >> %s 2>&1", 4000, 4000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/life %d %d %d >> %s 2>&1", 3000, 2000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/life") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/life %d %d %d >> %s 2>&1", 3000, 2000, 1024, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/3dfd %d %d %d %d >> %s 2>&1", 200, 200, 200, 16, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/3dfd") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/3dfd %d %d %d %d >> %s 2>&1", 200, 200, 200, 16, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lbm_tang") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/lbm_tang %d %d %d %d %s %d %d >> %s 2>&1", 64, 100, 100, 130, "result_lbm", 0, 0, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lbm_tang") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/lbm_tang %d %d %d %d %s %d %d >> %s 2>&1", 64, 100, 100, 130, "result_lbm", 0, 0, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/heat_4D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/heat_4D_NP %d %d >> %s 2>&1", 70, 32, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/heat_4D_NP") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/heat_4D_NP %d %d >> %s 2>&1", 70, 32, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/apop") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/apop -s %d -t %d >> %s 2>&1", 2000000, 524288, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/apop") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/apop -s %d -t %d >> %s 2>&1", 2000000, 524288, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/rna") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/rna -r %d >> %s 2>&1", 300, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/rna") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/rna -r %d >> %s 2>&1", 300, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/psa_struct") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/psa_struct -r %d %d >> %s 2>&1", 100000, 100000, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/psa_struct") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/psa_struct -r %d %d >> %s 2>&1", 100000, 100000, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_trap_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lcs") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/trap/lcs -r %d %d >> %s 2>&1", 200000, 200000, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;

	sprintf(suffix, "_sawzoid_%d_core", NCORES) ;
	strcpy(file, dir) ;
	strcat(file, "/lcs") ;
	strcat(file, suffix) ;
	strcpy(cmd, prefix) ;
	sprintf(tmp, "~/research/git/pochoir/examples/tc_tune_sub_space_sub_time/sawzoid/lcs -r %d %d >> %s 2>&1", 200000, 200000, file) ;
	strcat(cmd, tmp) ;
	cout << cmd << endl ;
	system(cmd) ;
	
	sprintf(prefix, "CILK_NWORKERS=12 ") ; 
	NCORES=12 ;
	}

	sprintf(dir, "tc_tune_sub_space_sub_time_%s/", time) ;
	sprintf(cmd, "mv *sawzoid %s", dir) ;
	system (cmd) ;
	sprintf(cmd, "mv *trap %s", dir) ;
	system (cmd) ;

	return 0 ;
}
