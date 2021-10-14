#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 1000000
#define cutoff 3.5
#define exradii 100.0
#define Psize 2 /*number of particle*/

const int psize = Psize;
int total,total1;
double daxis[3],maxis[3],axis[3];
double center[psize][3];
double pi;
const double fixc = 94.0;
int massnum=1;
char elementlist[6][BUF] = {
"  1   12.011\n",
"  2    1.008\n",
"  3   15.999\n",
"  4  195.078\n",
"  5   32.060\n",
"  6   18.998\n"};

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void write_kbdump (int option[ATOM][2], double array[ATOM][6]);
void read_center();
void write_check (int option[ATOM][2], double array[ATOM][6]);
void write_inputrd (int option[ATOM][2], double array[ATOM][6],int tot);
void merge_sort (int array[], double value[], int left, int right);

int main(int argc,char **argv)
{
	int i,j,a,b,c,d,e,f,g,flag;
	double	(*pos1v)[6];
	int	(*atom1c)[2],*index,d_num,cont;
	double	mass,*tempv;
	double	x,y,z,r;
	char buf[BUF],buf2[BUF],file1[BUF],file2[BUF],*temp,command[BUF];
	FILE *fp1,*fp2,*fp3;

	sprintf(file1,"%s",argv[1]); /* dump */

	pi=4*atan(1);
	atom1c=malloc(sizeof(int)*ATOM*2);
	pos1v=malloc(sizeof(double)*ATOM*10);
	
	total1=read_dump(file1,atom1c,pos1v);
	index=malloc(sizeof(int)*total1);
    tempv=malloc(sizeof(double)*total1);
	
	for(i = 0; i <total1;i++){
		index[i]=i;
		tempv[i]=pos1v[i][2];
	}
	merge_sort(index,tempv, 0, total1-1);
	d_num=0;
	for(i = 0; i<total1;i++){
		a=index[i];
		for(j = i+1; j<total1;j++){
			b=index[j];
			z=pos1v[b][2]-pos1v[a][2];
			if (z > cutoff) break;
			x=pos1v[b][0]-pos1v[a][0];
			y=pos1v[b][1]-pos1v[a][1];
			r=sqrtl(x*x+y*y+z*z);
			if (r<cutoff) { atom1c[a][1]+=1; atom1c[b][1]+=1;}
		}
		if (atom1c[a][1] >= 5) atom1c[a][1]=5;
		if (atom1c[a][1] <= 3) d_num+=1;
	}
	total=total1-d_num;
	printf("the number of atom %d - %d = %d\n",total1,d_num,total);
	
    read_center();
	write_check(atom1c,pos1v);
	write_inputrd(atom1c,pos1v,total);
	write_kbdump(atom1c,pos1v);
	printf("kb.dump, dump.check, and newinput.rd are generated\n");
}

int read_dump (char *file, int option[ATOM][2], double arr[ATOM][6]) {
	int i, j, k, a, b, c,d;
    double x,y,z,vx,vy,vz,fx,fy,fz,q;
    FILE *fp;
    char buf[BUF];
    
    fp = fopen(file,"r");/*rd file*/	
	if( fp==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file);
		fclose(fp);
		exit(1);
	}
   	    
    for (i = 0; i < 4;i++){fgets( buf, BUF, fp );}
    sscanf( buf, "%d",&a);
    fgets( buf, BUF, fp );
    for (i = 0; i < 3;i++){
        fgets( buf, BUF, fp );
        sscanf( buf, "%lf %lf",&daxis[i], &maxis[i]);
        axis[i]=maxis[i]-daxis[i];
    }
    fgets( buf, BUF, fp );
    for (i = 0; i < a;i++){
        fgets( buf, BUF, fp );
/*        sscanf( buf, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                     &c,&b,&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&q);*/
    	sscanf( buf, "%d%d%lf%lf%lf%lf%lf%lf", &c,&b,&x,&y,&z,&vx,&vy,&vz);
//    	d=i;
    	d=c-1;
        option[d][0]=b;
        option[d][1]=0;
        arr[d][0]=x-daxis[0];
        arr[d][1]=y-daxis[1];
        arr[d][2]=z-daxis[2];
        arr[d][3]=vx;
        arr[d][4]=vy;
        arr[d][5]=vz;
    }
    fclose(fp);
	
	return a;
}

void read_center(){

	int i;
    double x, y, z;
    float center[psize]={};
    FILE *fp;
    char buf[BUF];
        
    fp = fopen("center.txt","r");  //center
    if( fp==NULL ){
   	    printf("Error : These is no %s !! This process will stop.\n","config.rd");
        fclose(fp);
        exit(1);
    }
	while( fgets( buf, BUF, fp )!=NULL ){
    fgets( buf, BUF, fp );
        for (i=0; i<psize; i++){
	        sscanf( buf, "%f %f %f", &x, &y, &z);
            center[i][0]=x;
            center[i][1]=y;
            center[i][2]=z;
        }
    }
}

void write_kbdump (int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen("kb.dump","w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	a=0;
    for (i = 0; i < total1;i++){
    	if (option[i][1]>3){
    		a++;
	    	fprintf(fp,"%d  1",a);
	    	fprintf(fp," %.3f %.3f %.3f\n",arr[i][0],arr[i][1],arr[i][2]);
    	}
	}
    fclose(fp);
}


void write_check (int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen("dump.check","w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total1);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	a=0;
    for (i = 0; i < total1;i++){
    	if (option[i][1]<100){
    		a++;
	    	fprintf(fp,"%d %2d",a,option[i][1]);
	    	fprintf(fp," %.3f %.3f %.3f\n",arr[i][0],arr[i][1],arr[i][2]);
    	}
	}
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], double arr[ATOM][6], int tot) {
    int i, j, k, a, b, c;
    const int ctotal = total ;
    int mask[ctotal];
    double ddx, ddy, ddz;
    double ddr;
    FILE *fp;
    fp = fopen("newinput.rd","w");
    
	fprintf(fp,"#cellx  0.000  %7.3lf\n",axis[0]);
	fprintf(fp,"#celly  0.000  %7.3lf\n",axis[1]);
	fprintf(fp,"#cellz  0.000  %7.3lf\n",axis[2]);
	fprintf(fp,"\n#masses %d\n",massnum);

    fprintf(fp, "\n#fix 1 rigid  x y z\n");
    
    for (i = 0; i < massnum;i++){fprintf(fp,"%s",elementlist[i]);}
//	fprintf(fp,"\n#walls 1\n");
//	fprintf(fp,"  %7.3lf  %7.3lf  %7.3lf  0.0  %7.3lf\n",axis[0]*0.5,axis[1]*0.5,axis[2]*0.5,exradii);
	fprintf(fp,"\n#atoms %d\n",total);
	a=0;
                        
    for (i = 0; i < total1;i++){
    	if (option[i][1]>3){
    		a++;
            mask[i]=0;
            if(option[i][0]==1 && i <= total1){
                for(j=0; j<2; j++){
                    ddx = arr[i][0]-center[j][0];
                    ddy = arr[i][1]-center[j][1];
                    ddz = arr[i][2]-center[j][2];
                    ddr = ddx*ddx+ddy*ddy+ddz*ddz;
                    if(ddr <= fixc*fixc ){
                        mask[i]=1;
                        option[i][1] = 30;
                       //change fixcarbon number  
                    }
                }
            }
        	fprintf(fp,"%8d  %d   %d",a,option[a][0],mask[i]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf",arr[i][0],arr[i][1],arr[i][2]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][3],arr[i][4],arr[i][5]);
    	}
	}
    fclose(fp);
}

void merge_sort (int array[], double value[], int left, int right) {
	int i, j, k, mid;
	int work[right+1]; 
	if (left < right) {
		mid = (left + right)/2; 
		merge_sort(array, value, left, mid);
		merge_sort(array, value, mid+1, right);
		for (i = mid; i >= left; i--) { work[i] = array[i]; }
		for (j = mid+1; j <= right; j++) {
			work[right-(j-(mid+1))] = array[j];
		}
		i = left; j = right;
		for (k = left; k <= right; k++) {
			if (value[work[i]] < value[work[j]]) { array[k] = work[i++]; }
			else                   { array[k] = work[j--]; }
		}
	}
}
