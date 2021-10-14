/* produced by Nobuki Ozawa ver. 1.00  25/10/2020  
   Prepare dump file for initial structures  */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



#define BUF 256
#define ATOM 2000000
#define	pi 3.141592653589

/* Input */
char header[BUF] = {"colored.pos."};
int molecules[4][2]={{353283,1},{326,14},{3,16827},{416,276}};
int mol_k = sizeof(molecules)/sizeof(molecules[0]);
double	mgrid=3.0;
double mmgrid[3];

int target1=13; //F: 6 O of SO3H: 3 O of H2O: 13

// 1 Zr 2 Y 3 O 4 Ni 5 H

//Groval
int grid[3];
int vtotal;
int total,tstep,interval;
double daxis[3],maxis[3],laxis[3];

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void merge_sort (int array[], double value[], int left, int right);
void read_config ();
void write_voxel  (char* file, double value[]);
void write_voxel2  (char* file, double value[], double array[ATOM][6]);

int main(int argc,char **argv){
	
	int i,j,k,a,b,c,d,e,f,g,flag;
	int is0,is1,is2,js0,js1,js2;
	int (*atomc)[2];
	short *prop,*pprop;
	double *voxev,*vvoxev;
	int nn[3];
	
	double	(*posv)[6];
	double	r,x0,y0,z0,dx,dy,dz;

	
	char buf1[BUF],file1[BUF],file2[BUF],*temp;
	
	read_config ();	

	for(k = 0; k <= tstep;k+=interval){
		atomc=malloc(sizeof(int)*ATOM*2);
		posv=malloc(sizeof(double)*ATOM*6);
	    sprintf(file1,"%s%d",header,k);
	    sprintf(file2,"densitycar%d",k);
		total=read_dump(file1,atomc,posv);
		vtotal=1;
		for (j = 0; j < 3;j++){
			grid[j]=(int)(laxis[j]/mgrid);
			mmgrid[j]=laxis[j]/(double)grid[j];
			vtotal*=grid[j];
		}		
		printf("step %d number of grid %d %d %d\n",k,grid[0],grid[1],grid[2]);
		printf("%f %f %f\n",mmgrid[0],mmgrid[1],mmgrid[2]);
        double deltav=mmgrid[0]*mmgrid[1]*mmgrid[2];
		
		
		voxev=malloc(sizeof(double)*vtotal);
		prop=malloc(sizeof(short)*vtotal);
		
		for(i = 0; i<vtotal;i++){voxev[i]=0.0;}
		

		int Pt_start=molecules[0][0];//Start of Pt
		int w_start=molecules[0][0]+molecules[1][0]*molecules[1][1];//Start of Water
		for(i = w_start; i < total; i++){
			if (atomc[i][1]==target1){
				for(j = 0; j < 3; j++){nn[j]=(int)(posv[i][j]/mmgrid[j]);}
				a=nn[0]+nn[1]*grid[0]+nn[2]*grid[0]*grid[1];
				voxev[a]+=1.0/deltav;
//				printf("check %d %d %d %d\n",a,nn[0],nn[1],nn[2]);
			}
		}
        
//			for(i = 0; i<vtotal;i++){printf("%lf\n",voxev[i]);}
		
		 write_voxel2(file2,voxev,posv);
		
		free(atomc);
		free(posv);
		free(voxev);
		free(prop);
		free(vvoxev);
		free(pprop);
	}

}

int read_dump (char *file, int option[ATOM][2], double arr[ATOM][6]) {
	int i, j, k, a, b, c;
    double x,y,z,vx,vy,vz,fx,fy,fz,q;
    FILE *fp;
    char buf[BUF];
    
    fp = fopen(file,"r");/*rd file*/
	if( fp==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file);
        fclose(fp);
        exit(1);
	}
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


void write_voxel (char *file, double value[]) {
    int i, j, k, a, b, c,d;
    double x,y,z;
    FILE *fp;
    fp = fopen(file,"w");
    
	d=0;
	for (i = 0; i < grid[0];i++){
    for (j = 0; j < grid[1];j++){
    for (k = 0; k < grid[2];k++){
    	b=k+j*grid[2]+i*grid[1]*grid[2];
    	if (value[b] > 1) d+=1;
    }}}
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",d);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",maxis[j]-daxis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	
	a=1;
    for (i = 0; i < grid[0];i++){
    for (j = 0; j < grid[1];j++){
    for (k = 0; k < grid[2];k++){
    	b=k+j*grid[2]+i*grid[1]*grid[2];
    	if (value[b] > 0) {
    		x=(double)i*mmgrid[0]+0.5*mmgrid[0];
    		y=(double)j*mmgrid[1]+0.5*mmgrid[1];
    		z=(double)k*mmgrid[2]+0.5*mmgrid[2];
	    	fprintf(fp,"%d %2d",a,value[b]);
	    	fprintf(fp," %.2f %.2f %.2f\n",x,y,z);
    		a++;
    	}
	}}}
    fclose(fp);
}

void write_voxel2 (char *file, double value[], double arr[ATOM][6]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen(file,"w");
    
/*    fprintf(fp," CPMD CUBE FILE.\n");
	fprintf(fp," OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
	fprintf(fp," %4d    0.000000    0.000000    0.000000\n",1);
	fprintf(fp," %4d    %.6f    0.000000    0.000000\n",grid[0],mmgrid[0]);
	fprintf(fp," %4d    0.000000    %.6f    0.000000\n",grid[1],mmgrid[1]);
	fprintf(fp," %4d    0.000000    0.000000    %.6f\n ",grid[2],mmgrid[2]);
	fprintf(fp,"24   24.000000    0.000000    0.000000    0.000000\n");*/
	
    fprintf(fp,"voxel data\n");
	fprintf(fp,"   1.00000000000000     \n");
	fprintf(fp,"   %.6f    0.000000    0.000000\n",grid[0]*mmgrid[0]);
	fprintf(fp,"   0.000000    %.6f    0.000000\n",grid[1]*mmgrid[1]);
	fprintf(fp,"   0.000000    0.000000    %.6f\n",grid[2]*mmgrid[2]);
    fprintf(fp,"     C  Pt\n");
	a=molecules[0][0]*molecules[0][1];
	b=molecules[1][0]*molecules[1][1];
	c=a+b;
    fprintf(fp,"   %d   %d\n",a,b);
    fprintf(fp,"Direct\n");
	for (i = 0; i < c;i++){
		fprintf(fp,"  %8.6lf  %8.6lf  %8.6lf\n",(arr[i][0]-daxis[0])/laxis[0],(arr[i][1]-daxis[1])/laxis[1],(arr[i][2]-daxis[2])/laxis[2]);
	}
	
    fprintf(fp,"\n");
    fprintf(fp," %d %d %d\n",grid[0],grid[1],grid[2]);

    for (i = 1; i <= vtotal;i++){
	    fprintf(fp," %lf",value[i-1]);
    	if (i%6==0) fprintf(fp,"\n");;
	}
    fclose(fp);
}
