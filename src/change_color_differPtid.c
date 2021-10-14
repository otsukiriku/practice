/* produced by Nobuki Ozawa ver. 1.00  25/10/2020  
   Prepare dump file for initial structures  */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#define BUF 256
#define ATOM 3000000
#define ABORT_INIT 1
#define Ptatom 356 //number of Pt atoms in each cluster
int na = 241615;
int nb = 494551;

int pid, nprocs;
MPI_Status status;

int total,tstep,interval,ka;
double daxis[3],maxis[3],laxis[3];
double pi;

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][10]);
void read_config();
void write_dump (char *file, int option[][2], double array[][10]);

int main(int argc,char **argv){
	int i,j,k,a,b,c,d,e,f,g,count=0;
	int (*atomc)[2];
	double	(*posv)[10];
	double dx,dy,dz,r;
	int p,pp;
    FILE *fp1,*fp2;

/*start MPI -------------------------------------------------------------*/
	MPI_Init(&argc,&argv); /* Initialize the MPI environment */
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);  /* My processor ID */
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
/* Screen out Settings --------------------------------------------------*/
	
	atomc=malloc(sizeof(int)*ATOM*2);
	posv=malloc(sizeof(double)*ATOM*10);
	
	char buf1[BUF],file1[BUF],file2[BUF],*temp;

	pi=4*atan(1);
	read_config ();
	int kt=(tstep/interval)/nprocs;
	for(k = 0; k <= kt;k++){
		p=0;
		ka=(k*nprocs+pid)*interval;
	    sprintf(file1,"dump.pos.%d",ka);
		fp1 = fopen(file1,"r");
		sprintf(file2,"colored.pos.%d",ka);
		fp2 = fopen(file2,"r");
		if((fp1!=NULL)){
			total=read_dump(file1,atomc,posv);
//     change color for ovito 
            int idx = 40;
			for(i = 0; i < total;i++){
/*      	posv[i][0]: x(A) posv[i][1]: y(A) posv[i][2]: z(A)
        	atomc[i][0]:index for element atomc[i][1]:index for color in ovito
        	example        
        	if ((posv[i][2]<10)&&(posv[i][2]<20)) atomc[i][1]=21;
        
        	a=67842; b=3;
        	if (i+1>a) atomc[i][1]=21;
        	if ((i+1>a)&&(atomc[i][0]==b)) atomc[i][1]=22;*/
				if ((i+1 <= na) &&(atomc[i][0]==1)) {
					atomc[i][1]=21;
					dx=posv[i][0]-laxis[0]*0.5;
					dy=posv[i][1]-laxis[1]*0.5;
					dz=posv[i][2]-laxis[2]*0.5;
					r=dx*dx+dy*dy+dz*dz;
					if (r > 98.0*98.0) atomc[i][1]=22;
				}
				if ((i+1 > na)&&(i+1 < nb)&&(atomc[i][0]==2)) atomc[i][1]=12;
				if ((i+1 > na)&&(i+1 < nb)&&(atomc[i][0]==3)) atomc[i][1]=13;

// This part is for OHadd model
//                if ((i+1 > nb)&&(atomc[i][0]==2)) atomc[i][1]=14;
//                if ((i+1 > nb)&&(atomc[i][0]==3)) atomc[i][1]=15;
                
                if(atomc[i][0]==4){
                    atomc[i][1] = idx;
                    count++;
                    if(count==326){
                        count = 0;
                        idx++;
                    }
                }
			}
//     change color for ovito
			write_dump(file2,atomc,posv);
			printf("%s is generated.\n",file2);
		}
	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Allreduce(&p,&pp,nprocs,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
//		if ((pid==0)&&(pp == 0)) break;
	}
    MPI_Finalize();
}

int read_dump (char *file, int option[ATOM][2], double arr[ATOM][10]) {
	int i, j, k, a, b, c;
    double x,y,z,vx,vy,vz,fx,fy,fz,q;
    FILE *fp;
    char buf[BUF];
    
    fp = fopen(file,"r");/*rd file*/
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
        sscanf( buf, "%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
                     &c,&b,&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&q);
        option[c-1][0]=b;
        option[c-1][1]=b;
        arr[c-1][0]=x-daxis[0];
        arr[c-1][1]=y-daxis[1];
        arr[c-1][2]=z-daxis[2];
        arr[c-1][3]=vx;
        arr[c-1][4]=vy;
        arr[c-1][5]=vz;
        arr[c-1][6]=fx;
        arr[c-1][7]=fy;
        arr[c-1][8]=fz;
        arr[c-1][9]=q;
    }
    fclose(fp);
	
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

void write_dump (char *file, int option[][2], double arr[][10]) {
    int i, j, k, a, b, c;
    FILE *fp;
    fp = fopen(file,"w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",maxis[j]-daxis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z vx vy vz fx fy fz q\n");
    for (i = 0; i < total;i++){
	    fprintf(fp,"%d %2d",i+1,option[i][1]);
	    fprintf(fp," %.3f %.3f %.3f",arr[i][0],arr[i][1],arr[i][2]);
	    fprintf(fp," %.3f %.3f %.3f",arr[i][3],arr[i][4],arr[i][5]);
	    fprintf(fp," %.3f %.3f %.3f %.3f\n",arr[i][6],arr[i][7],arr[i][8],arr[i][9]);
	}
    fclose(fp);
}

