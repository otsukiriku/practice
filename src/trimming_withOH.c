#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 2000000
#define cutoff 0.3
#define exradii 107.0
#define tcutoff 3.2


int total,total1;
int tstep,interval;
int massnum;
char elementlist[6][BUF] = {
"  1   12.011\n",
"  2    1.008\n",
"  3   15.999\n",
"  4  195.078\n",
"  5   32.060\n",
"  6   18.998\n"};
double center[2][3]={{124.920,124.957,124.867},{364.920,124.957,124.867}};
int centnum = sizeof(center)/sizeof(center[0]);

double daxis[3],maxis[3],laxis[3],axis[3];
double pi;
double fixc = 94.0; //radius of fixcarbon

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void read_config();
void write_check (int option[ATOM][2], double array[ATOM][6]);
void write_inputrd (int option[ATOM][2], double array[ATOM][6], int tot);
void merge_sort (int array[], double value[], int left, int right);

int main(int argc,char **argv)
{
	int i,j,k,s,a,b,c,d,e,f,g,sw,flag,cont;
	double	(*pos1v)[6];
	double	pi,x,y,z,r,bo,sumbo,r0,r1;
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;
	int	(*atom1c)[2],*atomi;
	char buf1[BUF],buf2[BUF],file1[BUF],file2[BUF],*temp;
	char command[256];
	FILE *fp1,*fp2,*fp3,*fp4;
	int sp[9],db[2];
	read_config ();
	massnum=sizeof(elementlist)/sizeof(elementlist[0]);
	
	for(k = 0; k <= tstep;k+=interval){
		printf("step %d\n",k);
		
		atom1c=malloc(sizeof(int)*ATOM*2);
		pos1v=malloc(sizeof(double)*ATOM*6);
		
		sprintf(file1,"dump.pos.%d",k);
		sprintf(file2,"dump.bond.%d",k);
		total1=read_dump(file1,atom1c,pos1v);
		
		total=total1;
		for(i = 0; i <3;i++){axis[i]=laxis[i];}
		for(i = 0; i <9;i++){sp[i]=0;}
		for(i = 0; i <2;i++){db[i]=0;}
		for(i = 0; i <total;i++){
			atom1c[i][1]=0;
			if (atom1c[i][0]==4)atom1c[i][1]=4;
		}
		int COnum=0;
		int CHnum=0;
		
		fp3 = fopen(file2,"r");/*dump file*/
		while(fgets(buf2,BUF,fp3) != NULL){
			if(!strncmp(buf2,"Atom",4)){
				cont=0;
				sumbo=0.0;
				sscanf(buf2,"%s %d %d",&temp,&a,&b);
				if (atom1c[a-1][0]==1){
					for(i = 0; i < b;i++){
						fgets(buf2,BUF,fp3);
						sscanf(buf2,"%d%lf%lf%lf%lf",&c,&x,&y,&z,&bo);
						if ((bo>cutoff))cont+=1;
					    if ((bo>cutoff))sumbo+=bo;
						if ((atom1c[c-1][0]==2)&&(bo>cutoff)) atom1c[c-1][1]=2;
						if ((atom1c[c-1][0]==3)&&(bo>cutoff)) atom1c[c-1][1]=3;
						if ((atom1c[c-1][0]==2)&&(bo>cutoff)) CHnum+=1;
						if ((atom1c[c-1][0]==3)&&(bo>cutoff)) COnum+=1;
					}
//          40:sp3 30:sp2 20:sp 50:radical
					if (cont==4) atom1c[a-1][1]=40;
					if (cont==3) atom1c[a-1][1]=30;
					if (cont<=2) atom1c[a-1][1]=20;
					if ((cont==3)&&(sumbo>tcutoff)) atom1c[a-1][1]=30;
					if ((cont==3)&&(sumbo<tcutoff)) atom1c[a-1][1]=41;
					if ((cont==2)&&(sumbo>tcutoff)) atom1c[a-1][1]=20;
					if ((cont==2)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) atom1c[a-1][1]=31;
					if ((cont==2)&&(sumbo<tcutoff-1.0)) atom1c[a-1][1]=42;
					if ((cont<1)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) atom1c[a-1][1]=21;
					if ((cont<1)&&(sumbo<tcutoff-1.0)) atom1c[a-1][1]=32;
					if ((cont<1)&&(sumbo>tcutoff)) atom1c[a-1][1]=50;
				}
				
				if (atom1c[a-1][0]==4){
					for(i = 0; i < b;i++){
						fgets(buf2,BUF,fp3);
						sscanf(buf2,"%d%lf%lf%lf%lf",&c,&x,&y,&z,&bo);
						if ((atom1c[c-1][0]==2)&&(bo>cutoff)) atom1c[c-1][1]=2;
						if ((atom1c[c-1][0]==3)&&(bo>cutoff)) atom1c[c-1][1]=3;
					}
					
				}
			}		
		}
		fclose(fp3);
		
		fp3 = fopen(file2,"r");/*dump file*/
		while(fgets(buf2,BUF,fp3) != NULL){
			if(!strncmp(buf2,"Atom",4)){
				cont=0;
				sumbo=0.0;
				sscanf(buf2,"%s %d %d",&temp,&a,&b);
				if (atom1c[a-1][1]==3){
					for(i = 0; i < b;i++){
						fgets(buf2,BUF,fp3);
						sscanf(buf2,"%d%lf%lf%lf%lf",&c,&x,&y,&z,&bo);
						if ((atom1c[c-1][0]==2)&&(bo>cutoff)) {
							atom1c[c-1][1]=2;
						}
					}
				}
			}		
		}
		fclose(fp3);
		
		for(i = 0; i <total;i++){
//				if (atom1c[i][1]==40) sp[3]+=1;
				if (atom1c[i][1]==40) sp[0]+=1;
				if (atom1c[i][1]==40) db[0]+=1;
//				if ((cont==3)&&(sumbo>tcutoff)) atom1c[a-1][1]=30;
//				if (atom1c[i][1]==30) sp[2]+=1;
				if (atom1c[i][1]==30) sp[1]+=1;
				if (atom1c[i][1]==30) db[0]+=1;
//				if ((cont==2)&&(sumbo>tcutoff)) atom1c[a-1][1]=20;
//				if (atom1c[i][1]==20) sp[1]+=1;
				if (atom1c[i][1]==20) sp[2]+=1;
				if (atom1c[i][1]==20) db[0]+=1;
//				if ((cont<1)&&(sumbo>tcutoff)) atom1c[a-1][1]=50;
//				if (atom1c[i][1]==50) sp[0]+=1;
				if (atom1c[i][1]==50) sp[3]+=1;
//				if ((cont==3)&&(sumbo<tcutoff)) atom1c[a-1][1]=41;
//				if (atom1c[i][1]==41) sp[3]+=1;
				if (atom1c[i][1]==41) sp[4]+=1;
				if (atom1c[i][1]==41) db[1]+=1;
//				if ((cont==2)&&(sumbo<tcutoff-1.0)) atom1c[a-1][1]=42;
//				if (atom1c[i][1]==42) sp[3]+=1;
				if (atom1c[i][1]==42) sp[5]+=1;
				if (atom1c[i][1]==42) db[1]+=1;
//				if ((cont==2)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) atom1c[a-1][1]=31;
//				if (atom1c[i][1]==31) sp[2]+=1;
				if (atom1c[i][1]==31) sp[6]+=1;
				if (atom1c[i][1]==31) db[1]+=1;
//				if ((cont<1)&&(sumbo<tcutoff-1.0)) atom1c[a-1][1]=32;
//				if (atom1c[i][1]==32) sp[2]+=1;
				if (atom1c[i][1]==32) sp[7]+=1;
				if (atom1c[i][1]==32) db[1]+=1;
//				if ((cont<1)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) atom1c[a-1][1]=21;
//				if (atom1c[i][1]==21) sp[1]+=1;
				if (atom1c[i][1]==21) sp[8]+=1;
				if (atom1c[i][1]==21) db[1]+=1;
		}
//		a=sp[0]+sp[1]+sp[2]+sp[3];
		a=sp[0]+sp[1]+sp[2]+sp[3]+sp[4]+sp[5]+sp[6]+sp[7]+sp[8];
		db[0]=sp[0]+sp[1]+sp[2]+sp[3];
		db[1]=sp[4]+sp[5]+sp[6]+sp[7]+sp[8];
//		printf("number of carbon %d  total atom %d\n",a,total);
		printf("%d %d %d %d %d %d %d %d %d!\n",sp[0],sp[1],sp[2],sp[3],sp[4],sp[5],sp[6],sp[7],sp[8]);
		printf("noDB %d wDB %d\n",db[0],db[1]);
		printf("C-O %d C-H %d\n",COnum,CHnum);
//		printf("%.2lf %.2lf %.2lf %.2lf\n",sp[0]/(double)a*100.0,sp[1]/(double)a*100.0,sp[2]/(double)a*100.0,sp[3]/(double)a*100.0);
	}
	
	total=0;
	for(i = 0; i <total1;i++){
		x0=pos1v[i][0]-center[0][0];
		y0=pos1v[i][1]-center[0][1];
		z0=pos1v[i][2]-center[0][2];
		r0=sqrtl(x0*x0+y0*y0+z0*z0);
		x0=pos1v[i][0]-center[1][0];
		y0=pos1v[i][1]-center[1][1];
		z0=pos1v[i][2]-center[1][2];
		r1=sqrtl(x0*x0+y0*y0+z0*z0);
		if ((r0>exradii)&&(r1>exradii)&&(atom1c[i][0]!=4)) atom1c[i][1]=0;
		if (atom1c[i][1]!=0)total+=1;
	}
	write_inputrd(atom1c,pos1v, total);
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
    fp = fopen("dump.check","w");
    
    fprintf(fp,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n");
	fprintf(fp,"%d\n",total);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	a=0;
    for (i = 0; i < total1;i++){
    	if (option[i][1]!=0){
	    	fprintf(fp,"%d %2d",a+1,option[i][1]);
	    	fprintf(fp," %.3f %.3f %.3f\n",arr[i][0],arr[i][1],arr[i][2]);
    		a++;
    	}
	}
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], double arr[ATOM][6],int tot) {
    int i, j, k, a, b, c;
    double ddx, ddy, ddz, ddr;
    const int ctotal = total1 ;
    int mask[ctotal];
    FILE *fp;
    fp = fopen("newinput.rd","w");
    
	fprintf(fp,"#cellx  0.000  %7.3lf\n",axis[0]);
	fprintf(fp,"#celly  0.000  %7.3lf\n",axis[1]);
	fprintf(fp,"#cellz  0.000  %7.3lf\n",axis[2]);
	fprintf(fp,"\n#masses %d\n",massnum);
	for (i = 0; i < massnum;i++){fprintf(fp,"%s",elementlist[i]);}
//	fprintf(fp,"\n#walls 1\n");
//	fprintf(fp,"  %7.3lf  %7.3lf  %7.3lf  0.0  %7.3lf\n",axis[0]*0.5,axis[1]*0.5,axis[2]*0.5,exradii);
	fprintf(fp,"\n#atoms %d\n",total);
	a=0;
    for (i = 0; i < total1;i++){
            mask[i]=0;
            if(option[i][0]==1 && i <= total1){
                for(j=0; j<2; j++){
                    ddx = arr[i][0]-center[j][0];
                    ddy = arr[i][1]-center[j][1];
                    ddz = arr[i][2]-center[j][2];
                    ddr = ddx*ddx+ddy*ddy+ddz*ddz;
                    if(ddr <= fixc*fixc ){
                        mask[i]=1;
                        option[i][1] = 99;
                        //fixed carbon's option is changed to 99 in dump.check
                    } 
                }
            }
	    	if (option[i][1]!=0){
	        	fprintf(fp,"%8d  %d   %d",a+1,option[i][0],mask[i]);
	        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf",arr[i][0],arr[i][1],arr[i][2]);
	        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][3],arr[i][4],arr[i][5]);
	    		a++;
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
