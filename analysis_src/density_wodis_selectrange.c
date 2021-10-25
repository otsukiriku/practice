/* Density of Nafion */
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
int molecules[4][2]={{226566,1},{326,14},{3,16827},{416,276}};
int mol_k = sizeof(molecules)/sizeof(molecules[0]);
int w_start;
const int yrange = 400;
const int zrange = 430;
float range1 = 124.000;//smaller center x coordinate
float range2 = 304.000;//larger center x coordinate
double	mgrid=3.0;
double mmgrid[3];

const int target1=1; //F: 6 O of SO3H: 3 O of H2O: 13 C of Nafion:1

// 1 Zr 2 Y 3 O 4 Ni 5 H

//Groval
int grid[3];
int vtotal;
int total,tstep,interval;
double daxis[3],maxis[3],laxis[3];


int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void merge_sort (int array[], double value[], int left, int right);
void read_config ();
//float read_center(float center[]);
void write_voxel2  (char* file, int option[ATOM][2], double array[ATOM][6]);

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
	    sprintf(file2,"wodisdensity%d",k);
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
        
//		for(i = 0; i<vtotal;i++){printf("%lf\n",voxev[i]);}

        write_voxel2(file2,atomc,posv);
		
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
        option[c-1][0]=c;
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

/*float read_center(float center[]){

	int i, x, b, c;
    float center[particle_num]={};
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
        for (i=0; i<particle_num; i++){
	        sscanf( buf, "%f %f %f", &x, &b, &c);
            center[i]=x;
        }
    }
}
*/

void write_voxel2 (char *file, int option[ATOM][2], double arr[ATOM][6]) {
    int i, j, k, l, x, y, z;
    const int xrange=grid[0];
    const int yrange=grid[1];
    const int zrange=grid[2];
    int number[1000][1000][1000]={};
    int hantei[100]={};
    int count[100]={};
    FILE *fp;
    fp = fopen(file,"w");
	fprintf(fp, "#Cell range X=%f, Y=%f, Z=%f\n", grid[0]*mmgrid[0], grid[1]*mmgrid[1], grid[2]*mmgrid[2]);
    fprintf(fp, "\n#density data\n");

    //count number of particle in each cell
    for(i=0; i<total; i++){
        if(option[i][1] == target1){
            if(range1 < arr[i][0] && arr[i][0] < range2){
                x = round(arr[i][0])/mmgrid[0];
                y = round(arr[i][1])/mmgrid[1];
                z = round(arr[i][2])/mmgrid[2];
                number[x][y][z] += 1;
            }
            else{
                continue;
            }
        }
    }
    
    //count number of cell
    for(i=0; i<=xrange; i++)
        for(j=0; j<=yrange; j++)
           for(k=0; k<=zrange; k++){
                count[number[i][j][k]-1]++;
      //          printf("%d\n",number[x][y][z]);
            }

    for(l=0;l<100;l++){
    fprintf(fp, "%d %d\n",l,l*count[l]/abs(range2-range1));
    //fprintf(fp, "%.3f %d\n",hantei[l]*54.531/6.02*10/(mgrid*mgrid),count[l]);
    //94.342 is average mass of Nafion per one carbon
    //unit is g/cm-3
    }

    fclose(fp);

}



