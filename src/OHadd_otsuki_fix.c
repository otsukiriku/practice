#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 1000000
#define exradii 97.0
#define cutoff 0.3  // for bond population
#define tcutoff 3.2
#define dcutoff 2.2 // minimum distance between OH
#define ittr_limit 1000
#define particle_num 6

char file1[BUF] = {"dump.pos.0"};
char file2[BUF] = {"dump.bond.0"};

int total,total1,target;
double daxis[3],maxis[3],laxis[3],axis[3];
double pi;
double fixc = 94.0; //radius of fixcarbon
double center[particle_num][3]={};
double vdwradii[6]={1.70,1.2,1.52,1.75,1.80,1.47};
double mass[6]={12.011,1.008,15.999,195.078,32.060,18.998};
// 1 C 2 H 3 O 4 Pt 5 S 6 F
char elementlist[6][BUF] = {
"  1   12.011\n",
"  2    1.008\n",
"  3   15.999\n",
"  4  195.078\n",
"  5   32.060\n",
"  6   18.998\n"};
int centnum = sizeof(center)/sizeof(center[0]);
int vdw_n = sizeof(vdwradii)/sizeof(vdwradii[0]);
int massnum=sizeof(elementlist)/sizeof(elementlist[0]);

int rcut=23;
double	mgrid=3.5 ;//vdwradii[0]*2+0.1;
double OHratio=0.75; // ratio of OH group

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void read_config();
void write_check (int option[ATOM][2], double array[ATOM][6]);
void write_inputrd (int option[ATOM][2], double array[ATOM][6], int tot);
void write_kbdump (int option[ATOM][2], double array[ATOM][6]);
void merge_sort (int array[], double value[], int left, int right);
void surfacearea (int option[ATOM][2], double array[ATOM][6]);
int DBcheck (char* file, int option[ATOM][2]);
int addOH(int num, int t, double array[ATOM][6], double out[2][6]);
int addH(int num, int t, double array[ATOM][6], double out[2][6]);
void read_center();

int main(int argc,char **argv)
{
	int	i,j,k,a,b,c,d,e,f,g,s,t,u,flag;
	double	(*pos1v)[6],(*posv)[6],v;
	int	(*atom1c)[2],(*atomc)[2],*index,d_num,cont;
	double	mass,temp1v[2][6];
	double	x,y,z,r,x0,y0,z0,x1,y1,z1;
	char buf1[BUF],buf2[BUF],*temp,command[BUF];
	FILE *fp1,*fp2,*fp3;

    read_center();
    for(t=0;t<particle_num;t++)printf("%3.6lf %3.6lf %3.6lf\n",center[t][0],center[t][1],center[t][2]);

	srand((unsigned int)time(NULL));
//	srand(11);
	pi=4*atan(1);
	atom1c=malloc(sizeof(int)*ATOM*2);
	pos1v=malloc(sizeof(double)*ATOM*6);
	
	total1=read_dump(file1,atom1c,pos1v);
	total=total1;
	surfacearea(atom1c,pos1v);
	
	target=10;
	int db_total;
	db_total=DBcheck(file2, atom1c);
	
	int (*dblist)[2],(*ddblist)[2];
	int OHnum;
	OHnum=(int)(db_total*OHratio/2);
	total=total1+OHnum*3;
	atomc=malloc(sizeof(int)*total*2);
	posv=malloc(sizeof(double)*total*6);
	dblist=malloc(sizeof(int)*db_total*2);
	ddblist=malloc(sizeof(int)*2*OHnum*2);
	a=0;
	for (i = 0; i < total1;i++){
		if ((atom1c[i][1]%10!=0)&&(atom1c[i][1]>=20)){
			dblist[a][0]=i;
			dblist[a][1]=0;
			a+=1;
		}
		atomc[i][0]=atom1c[i][0];
		atomc[i][1]=atom1c[i][1];
		posv[i][0]=pos1v[i][0];
		posv[i][1]=pos1v[i][1];
		posv[i][2]=pos1v[i][2];
		posv[i][3]=pos1v[i][3];
		posv[i][4]=pos1v[i][4];
		posv[i][5]=pos1v[i][5];
	}
	t=0;
	c=0;
	while(t != 2*OHnum){
		flag=0;
		while (flag!=1){
			a=rand()%db_total;
			if (dblist[a][1]==0){
				ddblist[t][0]=dblist[a][0];
				ddblist[t][1]=1;
				flag=1;
				for (j = 0; j < t;j++){
					b=ddblist[j][0];
					
					x0=posv[dblist[a][0]][0]-posv[b][0];
					y0=posv[dblist[a][0]][1]-posv[b][1];
					z0=posv[dblist[a][0]][2]-posv[b][2];
					r=sqrtl(x0*x0+y0*y0+z0*z0);
					if ((r<dcutoff)&&(ddblist[j][1]==1)) {
//						printf("%d %d %lf %lf\n",a,b,r,dcutoff);
						flag=0;
						break;
					}
				}
								
				if (flag==1){
					d=c+total1;
					e=addOH(ddblist[t][0], d, posv, temp1v);
					if (e==ittr_limit){
//						printf("replace %d \n",a);
						flag=0;
					}
					atomc[d][0]=3;
					atomc[d][1]=3;
					posv[d][0]=temp1v[0][0];
					posv[d][1]=temp1v[0][1];
					posv[d][2]=temp1v[0][2];
					posv[d][3]=temp1v[0][3];
					posv[d][4]=temp1v[0][4];
					posv[d][5]=temp1v[0][5];
					if (flag==1) c++;
					d=c+total1;
					atomc[d][0]=2;
					atomc[d][1]=2;
					posv[d][0]=temp1v[1][0];
					posv[d][1]=temp1v[1][1];
					posv[d][2]=temp1v[1][2];
					posv[d][3]=temp1v[1][3];
					posv[d][4]=temp1v[1][4];
					posv[d][5]=temp1v[1][5];
					if (flag==1) c++;
				}
			}
			if (flag==1) dblist[a][1]=1;
			if (flag==1) t++;
		}		
		flag=0;
		while (flag!=1){
			a=rand()%db_total;
			if (dblist[a][1]==0){
				ddblist[t][0]=dblist[a][0];
				ddblist[t][1]=2;
				flag=1;
								
				if (flag==1){
					d=c+total1;
					e=addH(ddblist[t][0], d, posv, temp1v);
					if (e==ittr_limit){
//						printf("replace %d \n",a);
						flag=0;
					}
					atomc[d][0]=2;
					atomc[d][1]=2;
					posv[d][0]=temp1v[0][0];
					posv[d][1]=temp1v[0][1];
					posv[d][2]=temp1v[0][2];
					posv[d][3]=temp1v[0][3];
					posv[d][4]=temp1v[0][4];
					posv[d][5]=temp1v[0][5];
					if (flag==1) c++;
					
				}
			}
			
			if (flag==1) dblist[a][1]=2;
			if (flag==1) t++;
		}
		
		if ((t%100==0)&&(flag==1)) printf("checkpoint %d %d\n",t,2*OHnum);

	}

	write_inputrd(atomc,posv, total);
	write_check(atomc,posv);
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
    for (i = 0; i < total;i++){
    	if (option[i][1]>2){
    		a++;
	    	fprintf(fp,"%d  %d",a,option[i][0] );
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
	fprintf(fp,"%d\n",total);
	fprintf(fp,"ITEM: BOX BOUNDS pp pp pp\n");
	for(j = 0; j < 3;j++){fprintf(fp,"0.000 %.3f\n",axis[j]);}
	fprintf(fp,"ITEM: ATOMS id type x y z\n");
	a=0;
    for (i = 0; i < total;i++){
    	if (option[i][1]<100){
    		a++;
	    	fprintf(fp,"%d %2d",a,option[i][1]);
	    	fprintf(fp," %.3f %.3f %.3f\n",arr[i][0],arr[i][1],arr[i][2]);
    	}
	}
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], double arr[ATOM][6],int tot) {
    int i, j, k, a, b, c;
    double ddx, ddy, ddz, ddr;
    const int ctotal = total ;
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
    for (i = 0; i < total;i++){
    		a++;
            mask[i]=0;
            if(option[i][0]==1 && i <= total1){
                for(j=0; j<particle_num; j++){
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
        	fprintf(fp,"%8d  %d   %d",i+1,option[i][0],mask[i]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf",arr[i][0],arr[i][1],arr[i][2]);
        	fprintf(fp,"   %7.3lf   %7.3lf   %7.3lf\n",arr[i][3],arr[i][4],arr[i][5]);
    }
    fclose(fp);
}

/* produced by Nobuki Ozawa ver. 1.00  16/6/2021  */
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

void surfacearea (int option[ATOM][2], double arr[ATOM][6]) {
	int i,j,k,a,b,c,d,e,f,g,s,t,u,flag;
	double x0,y0,z0,r0,r[particle_num];
	double mmgrid[3];
	int grid[3];
	short *voxev,*prop;
	int nn[3];
    

	int vtotal=1;
	for (j = 0; j < 3;j++){
		grid[j]=(int)(axis[j]/mgrid);
		mmgrid[j]=axis[j]/(double)grid[j];
		vtotal*=grid[j];
	}		
	printf("number of grid %d %d %d\n",grid[0],grid[1],grid[2]);
	printf("%f %f %f\n",mmgrid[0],mmgrid[1],mmgrid[2]);
	double deltav=mmgrid[0]*mmgrid[1]*mmgrid[2];
	voxev=malloc(sizeof(short)*vtotal);
	prop=malloc(sizeof(short)*vtotal);
	for(i = 0; i<vtotal;i++){
		voxev[i]=0; prop[i]=0;
	}
	for(i = 0; i < total; i++){
		for(j = 0; j < 3; j++){nn[j]=(int)(arr[i][j]/mmgrid[j]);}
		a=nn[0]+nn[1]*grid[0]+nn[2]*grid[0]*grid[1];
//		printf("check %d %d %d %d\n",a,nn[0],nn[1],nn[2]);
		voxev[a]=1;
	}
	
	for(k = 1; k < grid[2]-1; k++){
	for(j = 1; j < grid[1]-1; j++){
	for(i = 1; i < grid[0]-1; i++){
		flag=0;
		for(s = -1; s < 2; s++){
		for(t = -1; t < 2; t++){
		for(u = -1; u < 2; u++){
			b=(i+s)+(j+t)*grid[0]+(k+u)*grid[0]*grid[1];
			if (voxev[b]>0) flag+=1;
		}}}
		a=i+j*grid[0]+k*grid[0]*grid[1];
		prop[a]=flag;
	}}}
	
//	for(i = 0; i<vtotal;i++){if(prop[i]!=0)printf("%d %d\n",i,prop[i]);}
	
	for(i = 0; i < total; i++){
	    flag = 0;
        if (option[i][0]==4)option[i][1]=4;
		if (option[i][0]==1){
			option[i][1]=11;
			
            for(k=0; k<particle_num; k++){
			    x0=arr[i][0]-center[k][0];
			    y0=arr[i][1]-center[k][1];
			    z0=arr[i][2]-center[k][2];
		    	r[k]=sqrtl(x0*x0+y0*y0+z0*z0);
			
		    	if (r[k]>exradii){
	    			for(j = 0; j < 3; j++){nn[j]=(int)(arr[i][j]/mmgrid[j]);}
    				a=nn[0]+nn[1]*grid[0]+nn[2]*grid[0]*grid[1];
			        flag +=1 ;
                }
			if (prop[a]<rcut && flag == 6) option[i][1]=10;
            }
		}
	}
}

int DBcheck (char *file, int option[ATOM][2]) {
	int i, j, k, a, b, c, d, cont;
	int sp[9],db[2];
    double x,y,z,bo,sumbo;

    FILE *fp,*temp;
    char buf[BUF];
    
	for(i = 0; i <9;i++){sp[i]=0;}
	for(i = 0; i <2;i++){db[i]=0;}
	
    fp = fopen(file,"r");/*rd file*/	
	if( fp==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file);
		fclose(fp);
		exit(1);
	}
	while(fgets(buf,BUF,fp) != NULL){
		if(!strncmp(buf,"Atom",4)){
			cont=0;
			sumbo=0.0;
			sscanf(buf,"%s %d %d",&temp,&a,&b);
			if (option[a-1][1]==target){
				for(i = 0; i < b;i++){
					fgets(buf,BUF,fp);
					sscanf(buf,"%d%lf%lf%lf%lf",&c,&x,&y,&z,&bo);
					if ((bo>cutoff))cont+=1;
				    if ((bo>cutoff))sumbo+=bo;
				}

				//          40:sp3 30:sp2 20:sp 50:radical
				if (cont==4) option[a-1][1]=40;
				if (cont==3) option[a-1][1]=30;
				if (cont<=2) option[a-1][1]=20;
				if ((cont==3)&&(sumbo>tcutoff)) option[a-1][1]=30;
				if ((cont==3)&&(sumbo<tcutoff)) option[a-1][1]=41;
				if ((cont==2)&&(sumbo>tcutoff)) option[a-1][1]=20;
				if ((cont==2)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) option[a-1][1]=31;
				if ((cont==2)&&(sumbo<tcutoff-1.0)) option[a-1][1]=42;
				if ((cont<1)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) option[a-1][1]=21;
				if ((cont<1)&&(sumbo<tcutoff-1.0)) option[a-1][1]=32;
				if ((cont<1)&&(sumbo>tcutoff)) option[a-1][1]=50;
			}
		}		
	}
	for(i = 0; i <total;i++){
//		if (option[i][1]==40) sp[3]+=1;
		if (option[i][1]==40) sp[0]+=1;
		if (option[i][1]==40) db[0]+=1;
//		if ((cont==3)&&(sumbo>tcutoff)) option[a-1][1]=30;
//		if (option[i][1]==30) sp[2]+=1;
		if (option[i][1]==30) sp[1]+=1;
		if (option[i][1]==30) db[0]+=1;
//		if ((cont==2)&&(sumbo>tcutoff)) option[a-1][1]=20;
//		if (option[i][1]==20) sp[1]+=1;
		if (option[i][1]==20) sp[2]+=1;
		if (option[i][1]==20) db[0]+=1;
//		if ((cont<1)&&(sumbo>tcutoff)) option[a-1][1]=50;
//		if (option[i][1]==50) sp[0]+=1;
		if (option[i][1]==50) sp[3]+=1;
//		if ((cont==3)&&(sumbo<tcutoff)) option[a-1][1]=41;
//		if (option[i][1]==41) sp[3]+=1;
		if (option[i][1]==41) sp[4]+=1;
		if (option[i][1]==41) db[1]+=1;
//		if ((cont==2)&&(sumbo<tcutoff-1.0)) option[a-1][1]=42;
//		if (option[i][1]==42) sp[3]+=1;
		if (option[i][1]==42) sp[5]+=1;
		if (option[i][1]==42) db[1]+=1;
//		if ((cont==2)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) option[a-1][1]=31;
//		if (option[i][1]==31) sp[2]+=1;
		if (option[i][1]==31) sp[6]+=1;
		if (option[i][1]==31) db[1]+=1;
//		if ((cont<1)&&(sumbo<tcutoff-1.0)) option[a-1][1]=32;
//		if (option[i][1]==32) sp[2]+=1;
		if (option[i][1]==32) sp[7]+=1;
		if (option[i][1]==32) db[1]+=1;
//		if ((cont<1)&&(sumbo<tcutoff)&&(sumbo>tcutoff-1.0)) option[a-1][1]=21;
//		if (option[i][1]==21) sp[1]+=1;
		if (option[i][1]==21) sp[8]+=1;
		if (option[i][1]==21) db[1]+=1;
	}
//		a=sp[0]+sp[1]+sp[2]+sp[3];
	a=sp[0]+sp[1]+sp[2]+sp[3]+sp[4]+sp[5]+sp[6]+sp[7]+sp[8];
	db[0]=sp[0]+sp[1]+sp[2]+sp[3];
	db[1]=sp[4]+sp[5]+sp[6]+sp[7]+sp[8];
//		printf("number of carbon %d  total atom %d\n",a,total);
	printf("%d %d %d %d %d %d %d %d %d!\n",sp[0],sp[1],sp[2],sp[3],sp[4],sp[5],sp[6],sp[7],sp[8]);
	printf("noDB %d wDB %d\n",db[0],db[1]);
	
	return db[1];
}


int addOH(int num, int t, double arr[ATOM][6], double out[2][6]) {
	int i, j, k, a, b, c, d,flag;
	int ittr;
	double	r,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
	double rcut=6.0;
	
	a=num;
	x1=0.0; y1=0.0; z1=0.0;
	for (i = 0; i < t;i++){
		x0=arr[a][0]-arr[i][0];
		y0=arr[a][1]-arr[i][1];
		z0=arr[a][2]-arr[i][2];
		r=sqrtl(x0*x0+y0*y0+z0*z0);
		if ((r<rcut)&&(i!=a)) {x1+=x0/r; y1+=y0/r; z1+=z0/r;}
	}
	
	flag=1;
	ittr=0;
	while (flag!=0){
		flag=0;
		r=sqrtl(x1*x1+y1*y1+z1*z1);
		x2=x1/r*1.6+arr[a][0]; 
		y2=y1/r*1.6+arr[a][1]; 
		z2=z1/r*1.6+arr[a][2];
		x3=x1/r*2.5+arr[a][0]; 
		y3=y1/r*2.5+arr[a][1];
		z3=z1/r*2.5+arr[a][2];
	    for (i = 0; i < t;i++){
			x0=x2-arr[i][0];
			y0=y2-arr[i][1];
			z0=z2-arr[i][2];
			r=sqrtl(x0*x0+y0*y0+z0*z0);
	    	if ((r<1.0)&&(i!=a)) {
//	    		printf("warning O in OH %d %d %.3lf\n",ittr,a,r);
	    		flag=1;
	    		break;
	    	}
	    	x0=x3-arr[i][0];
			y0=y3-arr[i][1];
			z0=z3-arr[i][2];
			r=sqrtl(x0*x0+y0*y0+z0*z0);
	    	if ((r<1.5)&&(i!=a)) {
//	    		printf("warning H in OH %d %d %.3lf \n",ittr,a,r);
	    		flag=1;
	    		break;
	    	}
	    }
		if (flag==1) {
			d=-500+rand()%1000;
			x1+=d/1000.0;
			d=-500+rand()%1000;
			y1+=d/1000.0;
			d=-500+rand()%1000;
			z1+=d/1000.0;
		}
		
		ittr+=1;
		if(ittr==ittr_limit) flag=0;
	}
	
	r=sqrtl(x1*x1+y1*y1+z1*z1);
	out[0][0]=x2; out[0][1]=y2; out[0][2]=z2;
	out[0][3]=1.0*x1/r/100.0; out[0][4]=1.0*y1/r/100.0; out[0][5]=1.0*z1/r/100.0;
	out[1][0]=x3; out[1][1]=y3; out[1][2]=z3;
	out[1][3]=1.0*x1/r/100.0; out[1][4]=1.0*y1/r/100.0; out[1][5]=1.0*z1/r/100.0;
	
	return ittr;
}

int addH(int num, int t, double arr[ATOM][6], double out[2][6]) {
	int i, j, k, a, b, c, d,flag;
	int ittr;
	double	r,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
	double rcut=4.0;
		
	a=num;
	x1=0.0; y1=0.0; z1=0.0;
	for (i = 0; i < t;i++){
		x0=arr[a][0]-arr[i][0];
		y0=arr[a][1]-arr[i][1];
		z0=arr[a][2]-arr[i][2];
		r=sqrtl(x0*x0+y0*y0+z0*z0);
		if ((r<rcut)&&(i!=a)) {x1+=x0/r; y1+=y0/r; z1+=z0/r;}
	}
	
	flag=1;
	ittr=0;
	while (flag!=0){
		flag=0;
		r=sqrtl(x1*x1+y1*y1+z1*z1);
		x2=x1/r*1.2+arr[a][0]; 
		y2=y1/r*1.2+arr[a][1]; 
		z2=z1/r*1.2+arr[a][2];
	    for (i = 0; i < t;i++){
			x0=x2-arr[i][0];
			y0=y2-arr[i][1];
			z0=z2-arr[i][2];
			r=sqrtl(x0*x0+y0*y0+z0*z0);
	    	if ((r<1.5)&&(i!=a)) {
//	    		printf("warning O in OH %d %d %.3lf\n",ittr,a,r);
	    		flag=1;
	    		break;
	    	}
	    }
		if (flag==1) {
			d=-500+rand()%1000;
			x1+=d/1000.0;
			d=-500+rand()%1000;
			y1+=d/1000.0;
			d=-500+rand()%1000;
			z1+=d/1000.0;
		}
		
		ittr+=1;
		if(ittr==ittr_limit) flag=0;
	}
	
	r=sqrtl(x1*x1+y1*y1+z1*z1);
	out[0][0]=x2; out[0][1]=y2; out[0][2]=z2;
	out[0][3]=1.0*x1/r/100.0; out[0][4]=1.0*y1/r/100.0; out[0][5]=1.0*z1/r/100.0;

	return ittr;
	
}

void read_center(){

	int i, x;
    double a, b, c;
    FILE *fp;
    char buf[BUF];
        
    fp = fopen("center.txt","r");  //center
    if( fp==NULL ){
   	    printf("Error : These is no %s !! This process will stop.\n","center.txt");
        fclose(fp);
        exit(1);
    }
    while( fgets( buf, BUF, fp )!=NULL ){
        for (i=0; i<particle_num; i++){
	        sscanf(buf, "%d %lf %lf %lf",&x, &a, &b, &c);
            center[x][0]=a;
            center[x][1]=b;
            center[x][2]=c;
        }
    }
    fclose(fp);
}
