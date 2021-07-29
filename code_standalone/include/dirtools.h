#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>

int prepare_directory(const char * mydir){
	DIR* dir = opendir(mydir);
	if (dir) {
		printf("Directory exists!\n");
    		closedir(dir);
		return 1;
	} else{
		if (ENOENT == errno){
			printf("Directory does not exist! Proceeding to create...\n");
			if(mkdir(mydir, 0757) != 0){
				printf("Error creating directory!\n");
				return -2;
			}
			else{
				printf("Directory created successfully!\n");
				return 0;
			}
		}
		else {
			printf("Error detecting directory.\n");
			return -1;
		}
	}
}
