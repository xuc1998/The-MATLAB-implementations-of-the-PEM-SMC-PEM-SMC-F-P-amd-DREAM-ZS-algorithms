#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h> 
#include <math.h>
#include <libgen.h> // For dirname
#include <limits.h> // For PATH_MAX 


int main(int argc, char **argv)
{
  
   char exe_path[PATH_MAX];//用来存储可执行文件的绝对路径
   char *exe_dir;//用来存储可执行文件所在的目录
   char script_path[PATH_MAX];//用于存储jobclm.csh的完整路径
   char chmod_command[PATH_MAX+10];// 用于存储chmod命令

// 通过/proc/self/exe获取当前可执行文件的路径
  ssize_t len=readlink("/proc/self/exe",exe_path,sizeof(exe_path)-1);

// 获取当前工作目录
  if (len !=-1){
   
   exe_path[len]='\0'; //添加字符串结束符
   exe_dir=dirname(exe_path);//获取可执行文件所在的目录

  // 打印可执行文件所在目录，确认路径是否正确
//  printf("Executable directory:%s\n",exe_dir);

  // 构建jobclm.csh的完整路径
  snprintf(script_path,sizeof(script_path),"%s/jobclm.csh",exe_dir);


  // 构建chmod命令，给绝对路径的jobclm.csh文件添加执行权限
  snprintf(chmod_command,sizeof(chmod_command),"chmod u+x %s",script_path);

  // 给.csh文件添加执行权限
  system(chmod_command);

  // 使用绝对路径运行jobclm.csh
  system(script_path);
  
}
else{
    perror("readlink error");
    return 1;
} 
return 0;   
}

