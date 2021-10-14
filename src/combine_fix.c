#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define ATOM 2000000
#define exradius 98.0
#define temper 300.0	/* [K] */
#define	AVOGADRO    6.02205e23	    /* [/mol] */
#define	BOLTZMANN   1.380658e-23    /* [J/K] */
#define	velocity   0.5    /* [km/s] */
#define Psize 2 /*number of paori*/

char file1[BUF]="ptcb2_120.dump";
char file2[BUF]="nafionw.dump";


const int psize = Psize;
int total,total1,total2;
int tstep,interval;
int unit[3] = {4,4,4};
double fixc = 94.0;
double rmax[2]={155.0,160.0};
double rmin[2]={98.0,98.0};
double vdwradii[6]={1.70,1.2,1.52,1.75,1.80,1.47};
double mass[6]={12.011,1.008,15.999,195.078,32.060,18.998};
// 1 C 2 H 3 O 4 Pt 5 S 6 F
int molecules[2][2]={{3,41472},{416,864}};
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
double  center[psize][3];

int  read_dump (char* file, int option[ATOM][2], double array[ATOM][6]);
void read_config();
void write_check (int option[ATOM][2], double array[ATOM][6]);
void write_inputrd (int option[ATOM][2], int tot, double array[ATOM][6], double lcenter[Psize][3]);
void write_center (double center[psize][3]);
void merge_sort (int array[], double value[], int left, int right);
void minimum (double array[ATOM][6],double constant);

int main(int argc,char **argv)
{
	int i,j,k,s,t,u,a,b,c,d,e,f,g,sw,flag,cont;
	int m,n,flag1,flag2,panum,pa;
    double	(*pos0v)[6],(*pos1v)[6],(*pos2v)[6],(*pos3v)[6];
	double	pi,x,y,z,r,rv,tmass,phi,theta,rr;
	double  vx,vy,vz;
	double	axis1[3],axis2[3],(*paori)[3];
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;
	double	xx,yy,zz;

	int	(*atom0c)[2],(*atom1c)[2],(*atom2c)[2],(*atom3c)[2];
	int *del_F, modify,vd[3];
	char buf1[BUF],buf2[BUF],*temp;
	char command[256];
	FILE *fp1,*fp2,*fp3,*fp4;
	
	srand((unsigned int)time(NULL));
	atom1c=malloc(sizeof(int)*ATOM*2);
	pos1v=malloc(sizeof(double)*ATOM*6);
	atom2c=malloc(sizeof(int)*ATOM*2);
	pos2v=malloc(sizeof(double)*ATOM*6);
	
    total1=read_dump(file1,atom1c,pos1v);
	for(i = 0; i < 3;i++){axis1[i]=laxis[i];}
    total2=read_dump(file2,atom2c,pos2v);
	for(i = 0; i < 3;i++){axis2[i]=laxis[i];}
	for(i = 0; i < 3;i++){axis[i]=axis2[i]*unit[i];}
	for(i = 0; i < total1;i++){
		pos1v[i][0]=pos1v[i][0]-0.5*axis1[0]+0.5*axis[0];
		pos1v[i][1]=pos1v[i][1]-0.5*axis1[1]+0.5*axis[1];
		pos1v[i][2]=pos1v[i][2]-0.5*axis1[2]+0.5*axis[2];
		dx=pos1v[i][0]-0.5*axis[0];
		dy=pos1v[i][1]-0.5*axis[1];
		dz=pos1v[i][2]-0.5*axis[2];
		r=sqrt(dx*dx+dy*dy+dz*dz);
		if (atom1c[i][0]==1) atom1c[i][1]=21;
		if ((atom1c[i][0]==1)&&(r>exradius)) atom1c[i][1]=22;
	}
	
	fp4 = fopen("center.txt","r");
	if( fp4==NULL ){
		printf("Warnning : These is no center.txt \n");
		panum=1;
		paori=malloc(sizeof(double)*panum * 3);
		paori[0][0]=0.5*axis[0]; paori[0][1]=0.5*axis[1]; paori[0][2]=0.5*axis[2];
	}else{
		fp4 = fopen("center.txt","r");
		i=0;
		while( fgets( buf1, BUF, fp4 )!=NULL ){
			sscanf( buf1, "%d%*s",&a);
			panum=a;
			paori=malloc(sizeof(double)*panum * 3);
			for(i = 0; i < panum;i++){
				fgets( buf1, BUF, fp4 );
				sscanf( buf1, "%d%lf%lf%lf",&a,&x,&y,&z);
		    	paori[i][0]=x-0.5*axis1[0]+0.5*axis[0];
				paori[i][1]=y-0.5*axis1[1]+0.5*axis[1];
				paori[i][2]=z-0.5*axis1[2]+0.5*axis[2];
                center[i][0]=paori[i][0];
                center[i][1]=paori[i][1];
                center[i][2]=paori[i][2];
				printf("%d %f %f %f \n",panum,paori[i][0],paori[i][1],paori[i][2]);
            }
		}
	}
	fclose(fp4);
	
	c=0;
//  correct coordination
	for(k = 0; k <mol_n;k++){
		for(i = 0; i <molecules[k][1];i++){			
			tmass=0.0;
			x0=0.0; y0=0.0; z0=0.0;
			for(j = 0; j <molecules[k][0];j++){
				a=c+j+i*molecules[k][0];
				d=atom2c[a][0];
				tmass+=mass[d];
				x0+=mass[d]*pos2v[a][0];
				y0+=mass[d]*pos2v[a][1];
				z0+=mass[d]*pos2v[a][2];
			}
			x0=x0/tmass; y0=y0/tmass; z0=z0/tmass;
			
			for(j = 1; j <molecules[k][0];j++){
				rr=9999.0;
				a=c+j+i*molecules[k][0];
				for(s = -1; s < 2; s++){
				for(t = -1; t < 2; t++){
				for(u = -1; u < 2; u++){
					dx=pos2v[a][0]+s*axis2[0]-x0;
					dy=pos2v[a][1]+t*axis2[1]-y0;
					dz=pos2v[a][2]+u*axis2[2]-z0;
					r=sqrt(dx*dx+dy*dy+dz*dz);
					if (r<rr){
						vd[0]=s; vd[1]=t; vd[2]=u;
						rr=r;
					}
				}}}
				pos2v[a][0]+=(double)vd[0]*axis2[0];
				pos2v[a][1]+=(double)vd[1]*axis2[1];
				pos2v[a][2]+=(double)vd[2]*axis2[2];
//				printf("%d\n",a);
			}
		}
		c+=molecules[k][0]*molecules[k][1];
	}
	
	xx=axis2[0]*(rand()%1000)/1000.0;
	yy=axis2[1]*(rand()%1000)/1000.0;
	zz=axis2[2]*(rand()%1000)/1000.0;
//	xx=0.0; yy=0.0; zz=0.0;
	
	fp1 = fopen("temp.rd","w");
	c=0;
	int mol_total=0;
	for(k = 0; k <mol_n;k++){
		modify=0;
		for(i = 0; i <molecules[k][1];i++){
			
			for(s = 0; s < unit[0]; s++){
			for(t = 0; t < unit[1]; t++){
			for(u = 0; u < unit[2]; u++){
			    tmass=0.0;
			    x0=0.0;
			    y0=0.0;
			    z0=0.0;
				for(j = 0; j <molecules[k][0];j++){
					a=c+j+i*molecules[k][0];
					d=atom2c[a][0];
					x1=pos2v[a][0]+s*axis2[0]-xx;
					y1=pos2v[a][1]+t*axis2[1]-yy;
					z1=pos2v[a][2]+u*axis2[2]-zz;
					tmass+=mass[d];
					x0+=mass[d]*x1;
					y0+=mass[d]*y1;
					z0+=mass[d]*z1;
			    }
				
				x0=x0/tmass;
				y0=y0/tmass;
				z0=z0/tmass;
				
				pa=0;
				dx=x0-paori[pa][0];
				dy=y0-paori[pa][1];
				dz=z0-paori[pa][2];
				r=sqrt(dx*dx+dy*dy+dz*dz);
				rr=r; 
				flag=0; flag1=0; flag2=0;
// hantei
				for(int j = 0; j < panum; j++) {
					dx=x0-paori[j][0];
					dy=y0-paori[j][1];
					dz=z0-paori[j][2];
					r=sqrt(dx*dx+dy*dy+dz*dz);
					if (r >= rmin[k]) flag1+=1;
					if (r <= rmax[k]) flag2+=1;
					if (r < rr) pa=j;
                    if (r < rr) rr=r;
				}
				if ((flag1==panum)&&(flag2>0)) flag=1;
				
				if (flag==1){
					modify+=1;
					for(j = 0; j <molecules[k][0];j++){
						a=c+j+i*molecules[k][0];
						x0=pos2v[a][0]+s*axis2[0]-xx;
						y0=pos2v[a][1]+t*axis2[1]-yy;
						z0=pos2v[a][2]+u*axis2[2]-zz;
						fprintf(fp1,"%d %.16lf %.16lf %.16lf\n",atom2c[a][0],x0,y0,z0);
					}
				}
// hantei
			}}}
		}
		c+=molecules[k][0]*molecules[k][1];
		molecules[k][1]=modify;
		mol_total+=molecules[k][1];
	}
	fclose(fp1);
	
	int total3=0;
	for(k = 0; k <mol_n;k++){
		total3+=molecules[k][0]*molecules[k][1];
	}
	del_F=malloc(sizeof(int)*mol_total);
	atom3c=malloc(sizeof(int)*total3*2);
	pos3v=malloc(sizeof(double)*total3*6);
	
	fp1 = fopen("temp.rd","r");
	for (i = 0; i < total3;i++){
		fgets( buf1, BUF, fp1 );
    	sscanf( buf1, "%d%lf%lf%lf", &a,&x,&y,&z);
        atom3c[i][0]=a;
        atom3c[i][1]=a;
        pos3v[i][0]=x;
        pos3v[i][1]=y;
        pos3v[i][2]=z;
    }
	fclose(fp1);
	
	total=total1;
	c=0;
	e=0;
	for(k = 0; k <mol_n;k++){
		for(i = 0; i <molecules[k][1];i++){
			for(j = 0; j <molecules[k][0];j++){
				a=c+j+i*molecules[k][0];
				
				for (s = 0; s <total1;s++){
					rv=vdwradii[atom1c[s][0]-1]+vdwradii[atom3c[a][0]-1];
					dx=pos3v[a][0]-pos1v[s][0];
					dy=pos3v[a][1]-pos1v[s][1];
					dz=pos3v[a][2]-pos1v[s][2];
					r=sqrt(dx*dx+dy*dy+dz*dz);
					if (r <rv)break;
				}
				if (r <rv){
					del_F[e]=0;
					break;
				}
				del_F[e]=1;
			}
			e+=1;
		}
		c+=molecules[k][0]*molecules[k][1];
	}	
	
	total=total1+total3;
	atom0c=malloc(sizeof(int)*total*2);
	pos0v=malloc(sizeof(double)*total*6);
	c=0;
	for(i = 0; i < total1;i++){
		pos0v[c][0]=pos1v[c][0];
		pos0v[c][1]=pos1v[c][1];
		pos0v[c][2]=pos1v[c][2];
		pos0v[c][3]=pos1v[c][3];
		pos0v[c][4]=pos1v[c][4];
		pos0v[c][5]=pos1v[c][5];
		atom0c[c][0]=atom1c[c][0];
		atom0c[c][1]=atom1c[c][1];
		c++;
	}
	b=0;
	e=0;
	total=total1;
	for(k = 0; k <mol_n;k++){
		modify=0;
		for(i = 0; i <molecules[k][1];i++){
			if (del_F[e]==1) {
			    tmass=0.0;
			    x0=0.0;
			    y0=0.0;
			    z0=0.0;
				for(j = 0; j <molecules[k][0];j++){
					a=b+j+i*molecules[k][0];
					d=atom3c[a][0];
					tmass+=mass[d];
					x0+=mass[d]*pos3v[a][0];
					y0+=mass[d]*pos3v[a][1];
					z0+=mass[d]*pos3v[a][2];
				}
				x0=x0/tmass-0.5*axis[0];
				y0=y0/tmass-0.5*axis[1];
				z0=z0/tmass-0.5*axis[2];
				r=sqrt(x0*x0+y0*y0+z0*z0);
				theta=asin(z0/r);
				phi=asin(y0/r/cos(theta));
				vx=velocity*cos(phi)*cos(theta);
				if (vx*x0>0) vx*=-1.0;
				vy=-1.0*velocity*sin(phi)*cos(theta);
				vz=-1.0*velocity*sin(theta);
				for(j = 0; j <molecules[k][0];j++){
					a=b+j+i*molecules[k][0];
					pos0v[c][0]=pos3v[a][0];
					pos0v[c][1]=pos3v[a][1];
					pos0v[c][2]=pos3v[a][2];
					pos0v[c][3]=vx/100.0;
					pos0v[c][4]=vy/100.0;
					pos0v[c][5]=vz/100.0;
					atom0c[c][0]=atom3c[a][0];
					atom0c[c][1]=atom3c[a][1];
//					printf("%d %.3lf %.3lf %.3lf \n",atom0c[c][0],pos0v[c][0],pos0v[c][1],pos0v[c][2]);
					total++;
					c++;
		    	}
				modify+=1;
			}
			e+=1;
		}
		printf("number of molecule %d %d\n",k+1,modify); 
        b+=molecules[k][0]*molecules[k][1];
	}
	
	minimum(pos0v, 5.0);
    write_center(center);
	write_inputrd(atom0c,total1,pos0v,center);
	write_check(atom0c,pos0v);
	
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
    for (i = 0; i < total;i++){
//    	if (option[i][1]<10){
    		a++;
            fprintf(fp,"%d %2d ",a,option[i][1]);
		    fprintf(fp,"%.3f %.3f %.3f\n",arr[i][0],arr[i][1],arr[i][2]);
//    	}
	}
    fclose(fp);
}

void write_center (double center[psize][3]) {
    int i;
    FILE *fp;
    fp = fopen("new_center.txt","w");
    
	fprintf(fp,"center\n");
	
    for (i = 0; i < psize;i++){
		fprintf(fp,"%d %.3lf %.3lf %.3lf\n",i,center[i][0],center[i][1],center[i][2]);	
    }
    fclose(fp);
}

void write_inputrd (int option[ATOM][2], int tot, double arr[ATOM][6], double center[Psize][3]){
    int i, j, k, a, b, c;
    const int ctotal = total ;
    int mask[ctotal];
    double ddx, ddy, ddz;
    double ddr;
    FILE *fp;
    fp = fopen("new_input.rd","w");
    
	fprintf(fp,"#cellx  0.000  %7.3lf\n",axis[0]);
	fprintf(fp,"#celly  0.000  %7.3lf\n",axis[1]);
	fprintf(fp,"#cellz  0.000  %7.3lf\n",axis[2]);
    
    fprintf(fp,"\n#masses %d\n",massnum);
	for (i = 0; i < massnum;i++){fprintf(fp,"%s",elementlist[i]);}
//fix inside if carbon
    fprintf(fp, "\n#fix 1 rigid  x y z\n");

//	fprintf(fp,"\n#walls 1\n");
//	fprintf(fp,"  %7.3lf  %7.3lf  %7.3lf  0.0  %7.3lf\n",axis[0]*0.5,axis[1]*0.5,axis[2]*0.5,exradii);
	fprintf(fp,"\n#atoms %d\n",total);
	a=0;
    for (i = 0; i < total;i++){
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
                        //fixed carbon's option is changed to 30 in dump.check
                    } 
                }
            }
        	fprintf(fp,"%8d  %d   %d",a,option[i][0],mask[i]);
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
    for (i =0; i<psize; i++){
        for (j=0; j < 3; j++){
            center[i][j]-=merge[j];
        }
    }
}
