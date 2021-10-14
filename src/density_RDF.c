/* Density of Nafion */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>



#define BUF 256
#define ATOM 2000000
#define	pi 3.141592653589
#define particle_num 2
#define	AVOGADRO    6.02205e23	    /* [/mol] */

/* Input */
char header[BUF] = {"colored.pos."};
int molecules[4][2]={{226566,1},{326,14},{3,16827},{416,276}};
int mol_k = sizeof(molecules)/sizeof(molecules[0]);
int w_start;
const int yrange = 400;
const int zrange = 430;
double	mgrid=3.0;
double mmgrid[3];

const int target1=15;
//Center of RDF
//F: 6 / SofSO3H: 5 / OofH2O:13 HofH2O:12 / CofNafion:1 OofNafion15/
const int target2=12;
// HofH2O:12
const int target3=2;
//H of terminalH on Carbon:2
const int totalHterminal=11384;
//number of total H terminal of catalyst


// 1 Zr 2 Y 3 O 4 Ni 5 H

//Groval
const double cutoff =10*10; 
int grid[3];
int vtotal;
int total,tstep,interval;
double daxis[3],maxis[3],laxis[3];


int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void merge_sort (int array[], double value[], int left, int right);
void read_config ();
//void read_center();
void write_voxel  (char* file, int option[ATOM][2], double array[ATOM][6]);

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
	    sprintf(file2,"RDF%d",k);
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
//      read_center();
        write_voxel(file2,atomc,posv);
		
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

/*void read_center(){

	int i, x, b, c;
    FILE *fp;
    char buf[BUF];
        
    fp = fopen("center.txt","r");  //center
    if( fp==NULL ){
   	    printf("Error : These is no %s !! This process will stop.\n","center.txt");
        fclose(fp);
        exit(1);
    }
	while( fgets( buf, BUF, fp )!=NULL ){
    fgets( buf, BUF, fp );
        for (i=0; i<particle_num; i++){
	        sscanf( buf, "%lf %lf %lf", &x, &b, &c);
            center[i][0]=x;
            center[i][1]=b;
            center[i][2]=c;
        }
    }
}
*/

void write_voxel (char *file, int option[ATOM][2], double arr[ATOM][6]){
    int i, j, l, r, r100;
    const int xrange=grid[0];
    const int yrange=grid[1];
    const int zrange=grid[2];
    int HofH2O[1000]={}, Hofterminal[1000]={}, distance[1000]={}, sum=0;
    int count=0;
    double x,y,z,k,rr;
    FILE *fp;
    fp = fopen(file,"w");
	fprintf(fp, "#Cell range X=%f, Y=%f, Z=%f\n", grid[0]*mmgrid[0], grid[1]*mmgrid[1], grid[2]*mmgrid[2]);
    fprintf(fp, "\n#density data\n");
    fprintf(fp, "rfromOofNafion / HofH2O(num/angs^3) / HofterminalH(num/angs^3) /number_of_terminalH \n");

    for(i=0; i<total; i++){
        if(option[i][1] == target1 && option[i-1][1]!=1){
        count++;
            for(j=0; j<total; j++){
                if(option[j][1] == target2 ||option[j][1] == target3){
                    x = arr[j][0]-arr[i][0];
                    y = arr[j][1]-arr[i][1];
                    z = arr[j][2]-arr[i][2];
                    rr = x*x+y*y+z*z;
                    if(rr <= cutoff){
                        r100 = (int)100*sqrt(rr); 
                        distance[r100]++;
                        if(option[j][1] == target2) HofH2O[r100]++;
                        else if(option[j][1] == target3) Hofterminal[r100]++;
                    }
                }
            }
        }
    }

    //fprintf(fp,"number of O %lf \n", count);

    for(i=1; i<1000; i++){
        fprintf(fp, "%.2f %.3lf %.3lf %d %d %.1lf\n" ,(float)i/100,(double)10000*HofH2O[i]/(4.0/3.0*pi*(pow(i+1,3)-pow(i,3))),(double)10000*Hofterminal[i]/(4.0/3.0*pi*(pow(i+1,3)-pow(i,3))), Hofterminal[i], sum, (double)100*sum/totalHterminal);
        sum += Hofterminal[i]; 
        
    }
    
    fclose(fp);
}
