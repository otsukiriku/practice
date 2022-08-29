#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 3000000
#define cutoff 0.3

char file1[BUF]="completed_pos";
char file2[BUF]="completed_bond";

int total,total1,total2;
int tstep,interval;
int unit[3] = {1,1,1};
double vdwradii[6]={1.70,1.2,1.52,1.75,1.80,1.47};
double mass[6]={12.011,1.008,15.999,195.078,32.060,18.998};
// 1 C 2 H 3 O 4 Pt 5 S 6 F
int molecules[2][2]={{3,24000},{416,500}};
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
	double	(*pos1v)[6],(*bond1v)[20];
	double	pi,x,y,z,r,bo,rv,rr;
	double  vx,vy,vz;
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;

	int	(*atom1c)[2],(*bond1c)[20];
	int	modify,vd[3];
	char buf1[BUF],buf2[BUF],*temp;
	char command[256];
	FILE *fp1,*fp2,*fp3,*fp4;
	
	srand((unsigned int)time(NULL));
	atom1c=malloc(sizeof(int)*ATOM*2);
	pos1v=malloc(sizeof(double)*ATOM*6);

	total1=read_dump(file1,atom1c,pos1v);
	for(i = 0; i < 3;i++){axis[i]=laxis[i];}
	bond1c=malloc(sizeof(int)*total1*20);
	bond1v=malloc(sizeof(double)*total1*20);
	total=total1;
	
	fp1 = fopen(file2,"r");/*dump file*/
	while(fgets(buf1,BUF,fp1) != NULL){
		if(!strncmp(buf1,"Atom",4)){
			sscanf(buf1,"%s %d %d",&temp,&a,&b);
			atom1c[a-1][1]=b;
			for(i = 0; i < b;i++){
				fgets(buf1,BUF,fp1);
				sscanf(buf1,"%d%lf%lf%lf%lf",&c,&x,&y,&z,&bo);
				bond1c[a-1][i]=c-1;
				bond1v[a-1][i]=bo;
			}
		}
	}
	
	c=0;
	for(k = 0; k <mol_n;k++){
		for(i = 0; i <molecules[k][1];i++){
			for(j = 0; j <molecules[k][0];j++){
				a=c+j+i*molecules[k][0];
				b=atom1c[a][1];
				for(l =0 ; l < b;l++){
					d=bond1c[a][l];
					if ((bond1v[a][l]>cutoff)&&(a<d)){
						dx=pos1v[a][0]-pos1v[d][0];
						dy=pos1v[a][1]-pos1v[d][1];
						dz=pos1v[a][2]-pos1v[d][2];
						rr=sqrt(dx*dx+dy*dy+dz*dz);
						if (rr>20.0){
							for(s = -1; s < 2; s++){
							for(t = -1; t < 2; t++){
							for(u = -1; u < 2; u++){
								dx=pos1v[a][0]-pos1v[d][0]-s*axis[0];
								dy=pos1v[a][1]-pos1v[d][1]-t*axis[1];
								dz=pos1v[a][2]-pos1v[d][2]-u*axis[2];
								r=sqrt(dx*dx+dy*dy+dz*dz);
								if (r<rr){
									vd[0]=s; vd[1]=t; vd[2]=u;
									rr=r;
								}
							}}}
							pos1v[d][0]+=(double)vd[0]*axis[0];
							pos1v[d][1]+=(double)vd[1]*axis[1];
							pos1v[d][2]+=(double)vd[2]*axis[2];
						}
//						printf("%d %d %lf\n",a+1,d+1,rr);
					}
				}
			}
		}
		c+=molecules[k][0]*molecules[k][1];
	}	
	
	write_check(atom1c,pos1v);
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
    fp = fopen("dump.modified","w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z vx vy vz\n");
	a=0;
    for (i = 0; i < total;i++){
//    	if (option[i][1]<10){
    		a++;
	    	fprintf(fp,"%d %2d",a,option[i][0]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf",arr[i][0],arr[i][1],arr[i][2]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][3],arr[i][4],arr[i][5]);
//    	}
	}
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen("new_input.rd","w");
    
	fprintf(fp,"#cellx  0.000  %7.3lf\n",axis[0]);
	fprintf(fp,"#celly  0.000  %7.3lf\n",axis[1]);
	fprintf(fp,"#cellz  0.000  %7.3lf\n",axis[2]);
	fprintf(fp,"\n#masses %d\n",massnum);
	for (i = 0; i < massnum;i++){fprintf(fp,"%s",elementlist[i]);}
//	fprintf(fp,"\n#walls 1\n");
//	fprintf(fp,"  %7.3lf  %7.3lf  %7.3lf  0.0  %7.3lf\n",axis[0]*0.5,axis[1]*0.5,axis[2]*0.5,exradii);
	fprintf(fp,"\n#atoms %d\n",total);
	a=0;
    for (i = 0; i < total;i++){
    		a++;
        	fprintf(fp,"%8d  %d   0",a,option[i][0]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf",arr[i][0],arr[i][1],arr[i][2]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][3],arr[i][4],arr[i][5]);
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
