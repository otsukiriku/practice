#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 2000000
#define cutoff 0.3
#define d_ratio 80

char file1[BUF]="dump.pos.0";

int total,total1,total2;
int tstep,interval,dd;
int unit[3] = {1,1,1};
double vdwradii[6]={1.70,1.2,1.52,1.75,1.80,1.47};
double mass[6]={12.011,1.008,15.999,195.078,32.060,18.998};
// 1 C 2 H 3 O 4 Pt 5 S 6 F
int molecules[4][2]={{122296,1},{3,16757},{2,2506},{416,339}};
char elementlist[6][BUF] = {
"  1   12.011\n",
"  2    1.008\n",
"  3   15.999\n",
"  4  195.078\n",
"  5   32.060\n",
"  6   18.998\n"};
int vdw_n = sizeof(vdwradii)/sizeof(vdwradii[0]);
int mol_n = sizeof(molecules)/sizeof(molecules[0]);
int massnum=sizeof(elementlist)/sizeof(elementlist[0]);

double daxis[3],maxis[3],laxis[3],axis[3];
double pi;

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void read_config();
void write_check (int option[ATOM][2], double array[ATOM][6]);
void write_inputrd (int option[ATOM][2], double array[ATOM][6]);
void merge_sort (int array[], double value[], int left, int right);
void minimum (double array[ATOM][6],double constant);

int main(int argc,char **argv)
{
	int i,j,k,l,s,t,u,a,b,c,d,e,f,g,sw,flag,cont;
	double	(*pos1v)[6];
	double	pi,x,y,z,r,bo,rv,rr;
	double  vx,vy,vz;
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;

	int	(*atom1c)[2];
	int	modify,vd[3];
	char buf1[BUF],buf2[BUF],*temp;
	char command[256];
	FILE *fp1,*fp2,*fp3,*fp4;
	
	srand((unsigned int)time(NULL));
	atom1c=malloc(sizeof(int)*ATOM*2);
	pos1v=malloc(sizeof(double)*ATOM*6);

	total1=read_dump(file1,atom1c,pos1v);
	for(i = 0; i < 3;i++){axis[i]=laxis[i];}
	total=total1;
	
	dd=0; k=mol_n-1; c=0;
    for(i=0;i<k;i++){
        c+=molecules[i][0]*molecules[i][1];
    }
	for(i = 0; i <molecules[k][1];i++){
		for(j = 0; j <molecules[k][0];j++){
			a=c+j+i*molecules[k][0];
			if (atom1c[a][0]==5){
				b=rand()%100;
				if (b<d_ratio){
                   	atom1c[a][1]=15;
                   	atom1c[a+4][1]=99;
                   	dd++;
				}
			}
		}
	}
	printf("%d SO3- are generated\n",dd);
	int limit=dd;
	b=total-dd;
	dd=0; 
	write_inputrd(atom1c,pos1v);
	write_check(atom1c,pos1v);
	fp1 = fopen("new_input.rd","a");
	fp2 = fopen("dump.check","a");

	k=1;
	c=molecules[0][0]*molecules[0][1];
	d=molecules[1][1];
	for(i = 0; i <molecules[k][1];i++){
		a=c+i*molecules[k][0];
		if (rand()%d <= limit+100){
			x0=pos1v[a+1][0]-pos1v[a][0];
			x1=pos1v[a+2][0]-pos1v[a][0];
			y0=pos1v[a+1][1]-pos1v[a][1];
			y1=pos1v[a+2][1]-pos1v[a][1];
			z0=pos1v[a+1][2]-pos1v[a][2];
			z1=pos1v[a+2][2]-pos1v[a][2];
			x2=(x0+x1)*0.75; y2=(y0+y1)*0.75; z2=(z0+z1)*0.75;
			x3=pos1v[a][0]-x2;y3=pos1v[a][1]-y2;z3=pos1v[a][2]-z2;
			if (x3+0.0005<0) x3+=axis[0];if (y3+0.0005<0) y3+=axis[1];if (z3+0.0005<0) z3+=axis[2];
			if (x3+0.0005>=axis[0]) x3-=axis[0];if (y3+0.0005>=axis[1]) y3-=axis[1];if (z3+0.0005>=axis[2]) z3-=axis[2];
			dd++;
        	fprintf(fp1,"%8d  %d   0",b+dd,2);
        	fprintf(fp1,"   %7.3lf   %7.3lf   %7.3lf",x3,y3,z3);
        	fprintf(fp1,"     0.000     0.000     0.001\n");
			fprintf(fp2,"%d %2d",b+dd,12);
        	fprintf(fp2,"   %7.3lf   %7.3lf   %7.3lf\n",x3,y3,z3);
		}
		if (dd==limit) break;
	}
	if (limit!=dd) printf("calculation is failed\n");
	fclose(fp1);fclose(fp2);
}

int read_dump (char *file, int option[ATOM][2], double arr[ATOM][6]) {
	int i, j, k, a, b, c;
    double x,y,z,vx,vy,vz,fx,fy,fz,q;
    FILE *fp;
    char buf[BUF];
    
    fp = fopen(file,"r");/*rd file*/
    if( fp!=NULL ){
   	    
    for (i = 0; i < 4;i++){fgets( buf, BUF, fp );}
    sscanf( buf, "%d",&a);
    fgets( buf, BUF, fp );
    for (i = 0; i < 3;i++){
        fgets( buf, BUF, fp );
        sscanf( buf, "%lf %lf",&daxis[i], &maxis[i]);
        laxis[i]=maxis[i]-daxis[i];
    }
    fgets( buf, BUF, fp );
    for (i = 0; i < a;i++){
        fgets( buf, BUF, fp );
/*        sscanf( buf, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                     &c,&b,&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&q);*/
    	sscanf( buf, "%d%d%lf%lf%lf%lf%lf%lf", &c,&b,&x,&y,&z,&vx,&vy,&vz);
        option[c-1][0]=b;
        option[c-1][1]=b;
        arr[c-1][0]=x-daxis[0];
        arr[c-1][1]=y-daxis[1];
        arr[c-1][2]=z-daxis[2];
        arr[c-1][3]=vx;
        arr[c-1][4]=vy;
        arr[c-1][5]=vz;
    }
    fclose(fp);
    }
	
	return a;
}

void read_config () {
	int i, j, k, a, b, c;
    FILE *fp;
    char buf[BUF];
        
    fp = fopen("config.rd","r");/*rd file*/
    if( fp==NULL ){
   	    printf("Error : These is no %s !! This process will stop.\n","config.rd");
        fclose(fp);
        exit(1);
    }
	while( fgets( buf, BUF, fp )!=NULL ){
		if( !strncmp( buf, "TotalStep", 9 ) ){
			sscanf( buf, "%*s%d",&tstep);
		}
		if( !strncmp( buf, "FileStep", 8 ) ){
			sscanf( buf, "%*s%d",&interval);
		}
	}
    fclose(fp);
}

void write_check (int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen("dump.check","w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total-dd);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	a=0;
    for (i = 0; i < total;i++){
    	if (option[i][1]!=99){
    		a++;
	    	fprintf(fp,"%d %2d",a,option[i][1]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][0],arr[i][1],arr[i][2]);
    	}
	}
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen("new_input.rd","w");
    
	fprintf(fp,"#cellx  0.000  %7.6lf\n",axis[0]);
	fprintf(fp,"#celly  0.000  %7.6lf\n",axis[1]);
	fprintf(fp,"#cellz  0.000  %7.6lf\n",axis[2]);
	fprintf(fp,"\n#masses %d\n",massnum);
	for (i = 0; i < massnum;i++){fprintf(fp,"%s",elementlist[i]);}
//	fprintf(fp,"\n#walls 1\n");
//	fprintf(fp,"  %7.3lf  %7.3lf  %7.3lf  0.0  %7.3lf\n",axis[0]*0.5,axis[1]*0.5,axis[2]*0.5,exradii);
	fprintf(fp,"\n#atoms %d\n",total-dd);
	a=0;
    for (i = 0; i < total;i++){
    	if (option[i][1]!=99){
    		a++;
        	fprintf(fp,"%8d  %d   0",a,option[i][0]);
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


void minimum (double array[ATOM][6], double constant){
	int i, j;
	double mini[3], max[3],merge[3];
	
	for (i = 0; i < 3;i++){
		mini[i]=axis[i];
		max[i]=0.0;
		merge[i]=0.0;
	}
	
	for (i = 0; i < total;i++){
		for (j = 0; j < 3;j++){
			if (mini[j]>array[i][j]) mini[j]=array[i][j];
			if (max[j]<array[i][j]) max[j]=array[i][j];
		}
	}
	
	for (i = 0; i < 3;i++){
		merge[i]=mini[i];
		if (axis[i]-max[i] < mini[i]) merge[i]=axis[i]-max[i];
		merge[i]-=constant;
		axis[i]-=2.0*merge[i];
	}
	for (i = 0; i < total;i++){
		for (j = 0; j < 3;j++){
			array[i][j]-=merge[j];
		}
	}	
}
