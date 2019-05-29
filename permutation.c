#include <stdio.h>
#include <stdlib.h>
 
int perm(int index){	
	int temp;
	int numbers=3;
	int a[numbers+1], upto = 2, temp2;
	int vector[index*3],index = 0,i;
	for( temp2 = 1 ; temp2 <= numbers; temp2++){
		a[temp2]=-2;
	}
	a[numbers]=-2-1;
	temp=numbers;
	while(1){
		if(a[temp]==upto){
			temp--;
			if(temp==0)
				break;
		}
		else{
			a[temp]++;
			while(temp<numbers){
				temp++;
				a[temp]=-2;
			}
 
		printf("(");
			for( temp2 = 1 ; temp2 <= numbers; temp2++){
			    printf("%d", a[temp2]);
	
				//vector[index] = a[temp2];
			//	index++;
			}
			printf(")");
		}
	}
	
//	for(i=0;i<125*3;i++)
//	{
//	    if(vector[i]==1)
//	       vector[i] = -2;
//	  if(vector[i]==2)
	//     vector[i] = -1;
//	if(vector[i]==3)
//	   vector[i] = 0;
//	if(vector[i]==4)
//	    vector[i] = 1;
//	if(vector[i]==5)
//	   vector[i] = 2;
//	}
//	printf("{");
//	for(i=0;i<125*3;i++)
//	{
//	    printf("%d,",vector[i]);
//	}
//	printf("}");
	
	return 0;
}