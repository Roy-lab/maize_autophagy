#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
int 
main(int argc, const char** argv)
{
	if(argc!=2)
	{
		cout <<"Give me a number" << endl;
		exit(0);
	}
	printf("The number is %.2e\n",atof(argv[1]));
	return 0;
}
