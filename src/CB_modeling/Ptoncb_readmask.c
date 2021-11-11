/* produced by Nobuki Ozawa ver. 1.00  16/6/2020  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define sradii 10
#define radii 115
#define ATOM 3000000
#define	AVOGADRO    6.02205e23	    /* [/mol] */
#define ptnum 84 // the number of Pt particles
#define cutoff 70 // distance between Pt particles
#define velocity 0.1 // km/s
#define Ptclustnum 309
#define particle_num 6 //carbon_num
#define file1 "./dump.pos.0"
#define file2 "./Pt309.rd"

void merge_sort (int array[], double value[], int left, int right) {
	int i, j, k, mid;
	int work[right+1];  // 作業用配列
	if (left < right) {
		mid = (left + right)/2; // 真ん中
		merge_sort(array, value, left, mid);  // 左を整列
		merge_sort(array, value, mid+1, right);  // 右を整列
		for (i = mid; i >= left; i--) { work[i] = array[i]; } // 左半分
		for (j = mid+1; j <= right; j++) {
			work[right-(j-(mid+1))] = array[j]; // 右半分を逆順
		}
		i = left; j = right;
		for (k = left; k <= right; k++) {
			if (value[work[i]] < value[work[j]]) { array[k] = work[i++]; }
			else                   { array[k] = work[j--]; }
		}
	}
}

int main(int argc,char **argv)
{
	int i,j,k,s,a,b,c,d,e,f,g,sw,flag,cont,q,m;
	int total,ctotal,ptotal,panum;
    int mask[ATOM];
	double	(*pos1v)[6],(*pos2v)[3],axis[3],daxis[3],maxis[3],laxis[3];
	double	x,y,z,r,vx,vy,vz;
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;
	double	pi,theta,phi,ramda,(*paori)[3],(*ptmid)[3],(*ptr)[4];
	int	(*atom1c)[2],(*ptint)[4];
	char buf1[BUF],buf2[BUF],*temp;
	char command[256];
	FILE *fp1,*fp2,*fp3,*fp4;

	pi=4*atan(1);
	ptint=malloc(sizeof(int)*ptnum * 4);
	ptr=malloc(sizeof(double)*ptnum * 4);
	ptmid=malloc(sizeof(double)*ptnum * 3);
	srand((unsigned int)time(NULL));


	fp1 = fopen(file1,"r");/*rd file*/
	if( fp1==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file1);
		fclose(fp1);
		exit(1);
	}
	for (i = 0; i < 4;i++){fgets( buf1, BUF, fp1 );}
	sscanf( buf1, "%d",&ctotal);
	atom1c=malloc(sizeof(int)*ctotal*2);
	pos1v=malloc(sizeof(double)*ctotal * 6);
	fgets( buf1, BUF, fp1 );
	for (i = 0; i < 3;i++){
	    fgets( buf1, BUF, fp1 );
        sscanf( buf1, "%lf %lf",&daxis[i], &maxis[i]);
        axis[i]=maxis[i]-daxis[i];
	}
	fgets( buf1, BUF, fp1 );
//	printf("axis %lf %lf %lf\n",axis[0],axis[1],axis[2]);

    for (i = 0; i < ctotal;i++){
    	fgets( buf1, BUF, fp1 );
    	sscanf( buf1, "%d%d%lf%lf%lf%lf%lf%lf%lf%d",&a,&b,&x,&y,&z,&vx,&vy,&vz,&q,&m);
    	atom1c[i][1]=b;
    	pos1v[i][0]=x;
    	pos1v[i][1]=y;
    	pos1v[i][2]=z;
    	pos1v[i][3]=vx*100.0;
    	pos1v[i][4]=vy*100.0;
    	pos1v[i][5]=vz*100.0;
        mask[i]=m;
//    	printf("%d %lf %lf %lf\n",b,x,y,z);
    	for (j = 0; j < 3;j++){pos1v[i][j]-=daxis[j];}
//	printf("check %d %lf %lf %lf\n",a,x,y,z);
    }
	fclose(fp1);

	
	fp1 = fopen(file2,"r");
	if( fp1==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file1);
		fclose(fp1);
		exit(1);
	}
	pos2v=malloc(sizeof(double)*ptotal * 3);
	sw=0; i=0;
	while( fgets( buf1, BUF, fp1 )!=NULL ){
		if( !strncmp( buf1, "axis", 4 ) ){
			a = sscanf( buf1, "%*s%lf%lf%lf",&laxis[0], &laxis[1], &laxis[2]);
			laxis[0] = laxis[0] * 1e10;
			laxis[1] = laxis[1] * 1e10;
			laxis[2] = laxis[2] * 1e10;
		}
		if( !strncmp( buf1, "atoms", 5 ) ){
			sscanf( buf1, "%*s%d",&a);
			ptotal=a;
			pos2v=malloc(sizeof(double)*ptotal * 3);
		}
		if( !strncmp( buf1, "end", 3 ) )sw = 0;
		if(sw==2) {
			sscanf( buf1, "%d%d%d%d%d%lf%lf%lf",&e,&a,&b,&c,&d,&x,&y,&z);
			pos2v[i][0]=(x-0.5)*laxis[0];
			pos2v[i][1]=(y-0.5)*laxis[1];
			pos2v[i][2]=(z-0.5)*laxis[2]+radii;
//			printf("%d %f %f %f\n",i,pos2v[i][0],pos2v[i][1],pos2v[i][2]);
			i+=1;
		}		
		if( !strncmp( buf1, "elem", 4 ) )sw = 1;
		if( !strncmp( buf1, "part", 4 ) )sw = 2;
	}
	fclose(fp1);
	total=ctotal+ptotal*ptnum;
		
	fp1 = fopen("center.txt","r");
	if( fp1==NULL ){
		printf("Error : These is no center.txt !! This process will stop.\n");
		fclose(fp1);
		exit(1);
	}
	if((fp1 = fopen("center.txt","r"))==NULL){
		paori[i][0]=0.5*axis[0]; paori[i][1]=0.5*axis[1]; paori[i][2]=0.5*axis[2];
	}else{
		fp1 = fopen("center.txt","r");
		i=0;
		while( fgets( buf1, BUF, fp1 )!=NULL ){
			sscanf( buf1, "%d%*s",&a);
			panum=a;
			paori=malloc(sizeof(double)*panum * 3);
			for(i = 0; i < panum;i++){
				fgets( buf1, BUF, fp1 );
				sscanf( buf1, "%d%lf%lf%lf",&a,&x,&y,&z);
		    	paori[i][0]=x; paori[i][1]=y; paori[i][2]=z;
//				printf("%d %f %f %f \n",panum,paori[i][0],paori[i][1],paori[i][2]);
			}
		}
	fclose(fp1);
	}
	
//  Pt cluster positions	
	if((fp4 = fopen("ptcluster_cb.txt","r"))==NULL){
	for(i = 0; i < 4;i++){ptint[0][i]=0;}
	flag=0; cont=0;
	for(i = 0; i < ptnum;i++){
		ptr[i][0]=radii;
	}
	for(i = 1; i < ptnum;i++){
		while (flag != i) {
			ptint[i][0]=rand()%panum;
			ptint[i][1]=rand()%360;
			ptint[i][2]=rand()%360;
			ptint[i][3]=rand()%360;
			a=ptint[i][0];
			theta=pi*ptint[i][1]/180.0;
			phi=pi*ptint[i][2]/180.0;
			x1=paori[a][0]+radii*sin(theta)*sin(phi);
			y1=paori[a][1]+radii*cos(theta)*sin(phi);
			z1=paori[a][2]+radii*cos(phi);
			j=0; flag=0;
			while (j<i) {
				b=ptint[j][0];
				theta=pi*ptint[j][1]/180.0;
				phi=pi*ptint[j][2]/180.0;
				x2=paori[b][0]+radii*sin(theta)*sin(phi);
				y2=paori[b][1]+radii*cos(theta)*sin(phi);
				z2=paori[b][2]+radii*cos(phi);
				dx=fabs(x2-x1);
				if (dx>0.5*axis[0]){dx-=axis[0];}
				dy=fabs(y2-y1);
				if (dy>0.5*axis[1]){dy-=axis[1];}
				dz=fabs(z2-z1);
				if (dz>0.5*axis[2]){dz-=axis[2];}
				r=sqrt(dx*dx+dy*dy+dz*dz);
				printf("Number %3d %3d  %.3f %.3f\n", i+50,j+50,r,cutoff);
				if (r>cutoff) flag+=1;
				if (r<cutoff) break;
				j+=1;
			}
			
			a=ptint[i][0];
			theta=pi*ptint[i][1]/180.0;
			phi=pi*ptint[i][2]/180.0;
			ramda=pi*ptint[i][3]/180.0;
			for(j = 0; j < ptotal;j++){
				x0=pos2v[j][0]*cos(ramda)+pos2v[j][1]*sin(ramda);
				y0=-1.0*pos2v[j][0]*sin(ramda)+pos2v[j][1]*cos(ramda);
				x1=x0;
				y1=y0*cos(phi)+pos2v[j][2]*sin(phi);
				z1=-1.0*y0*sin(phi)+pos2v[j][2]*cos(phi);
				x2=x1*cos(theta)+y1*sin(theta);
				y2=-1.0*x1*sin(theta)+y1*cos(theta);
				z2=z1;
				x=x2+paori[a][0];
				y=y2+paori[a][1];
				z=z2+paori[a][2];
				if (x<0){x+=axis[0];}
				if (y<0){y+=axis[1];}
				if (z<0){z+=axis[2];}
				if (x>axis[0]){x-=axis[0];}
				if (y>axis[1]){y-=axis[1];}
				if (z>axis[2]){z-=axis[2];}
				for(k = 0; k < panum;k++){
					if (a!=k){
						dx=fabs(x-paori[k][0]);
						dy=fabs(y-paori[k][1]);
						dz=fabs(z-paori[k][2]);
						if (dx>0.5*axis[0]){dx-=axis[0];}
						if (dy>0.5*axis[1]){dy-=axis[1];}
						if (dz>0.5*axis[2]){dz-=axis[2];}
						r=sqrt(dx*dx+dy*dy+dz*dz);
						if (r < radii) break;
					}
				}
				if (r < radii-2) {
					printf("%d is contained in particles %d %lf A\n",i+50,k,r);
					flag=0;
					break;
				}
			}
			
			if (flag!=i){cont+=1;}
			if (cont==2000){i=1; cont=0;}
//		printf("%d %d %d %d\n",i,ptint[i][0],ptint[i][1],ptint[i][2]);
		}
		cont=0;
	}
		
	fp4 = fopen("ptcluster_cb.txt","w");
		for(i = 0; i < ptnum;i++){fprintf(fp4,"%3d %3d %3d %3d %3d %7.2f 0.0 0.0 0.0\n",50+i,ptint[i][0],ptint[i][1],ptint[i][2],ptint[i][3],ptr[i][0]);}	
	fclose(fp4);
	
	a=0;
	fp3 = fopen("dump.check","w");
	fprintf(fp3,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n",total);
	fprintf(fp3,"ITEM: BOX BOUNDS pp pp pp\n");
	for(i = 0; i < 3;i++){fprintf(fp3,"0.0  %10.7f\n",axis[i]);}
	fprintf(fp3,"ITEM: ATOMS id type x y z\n");
	for(i = 0; i < ctotal;i++){
		fprintf(fp3,"%d %d %.3f %.3f %.3f\n",a,atom1c[i][1],pos1v[i][0],pos1v[i][1],pos1v[i][2]);
		a+=1;
	}
	for(i = 0; i < ptnum;i++){
		b=ptint[i][0];
		c=ptint[i][1];
		d=ptint[i][2];
		e=ptint[i][3];
		theta=pi*(double)c/180.0;
		phi=pi*(double)d/180.0;
		ramda=pi*(double)e/180.0;
		for(j = 0; j < ptotal;j++){
			x0=pos2v[j][0]*cos(ramda)+pos2v[j][1]*sin(ramda);
			y0=-1.0*pos2v[j][0]*sin(ramda)+pos2v[j][1]*cos(ramda);
			z0=pos2v[j][2]-radii+ptr[i][0];
			x1=x0;
			y1=y0*cos(phi)+z0*sin(phi);
			z1=-1.0*y0*sin(phi)+z0*cos(phi);
			x2=x1*cos(theta)+y1*sin(theta);
			y2=-1.0*x1*sin(theta)+y1*cos(theta);
			z2=z1;
			x=x2+paori[b][0];
			y=y2+paori[b][1];
			z=z2+paori[b][2];
			if (x<0){x+=axis[0];}
			if (y<0){y+=axis[1];}
			if (z<0){z+=axis[2];}
			if (x>axis[0]){x-=axis[0];}
			if (y>axis[1]){y-=axis[1];}
			if (z>axis[2]){z-=axis[2];}

			fprintf(fp3,"%d %d %.3f %.3f %.3f\n",a,50+i,x,y,z);
			a+=1;
		}
	}

	}else{
		printf("Reading ptcluster_cb.txt\n");
		fp4 = fopen("ptcluster_cb.txt","r");
		i=0;
		while( fgets( buf1, BUF, fp1 )!=NULL ){
			sscanf( buf1, "%d%d%d%d%d%lf%lf%lf%lf",&a,&b,&c,&d,&e,&r,&x,&y,&z);
			ptint[i][0]=b; ptint[i][1]=c; ptint[i][2]=d; ptint[i][3]=e; ptr[i][0]=r;
			ptr[i][1]=x; ptr[i][2]=y; ptr[i][3]=z;
			i+=1;
		}
		fclose(fp4);
		a=0;
		fp3 = fopen("dump.check","w");
		fprintf(fp3,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\n",total);
		fprintf(fp3,"ITEM: BOX BOUNDS pp pp pp\n");
		for(i = 0; i < 3;i++){fprintf(fp3,"0.0  %10.7f\n",axis[i]);}
		fprintf(fp3,"ITEM: ATOMS id type x y z\n");
		for(i = 0; i < ctotal;i++){
			fprintf(fp3,"%d %d %.3f %.3f %.3f\n",a,atom1c[i][1],pos1v[i][0],pos1v[i][1],pos1v[i][2]);
			a+=1;
		}
		for(i = 0; i < ptnum;i++){
		b=ptint[i][0];
		c=ptint[i][1];
		d=ptint[i][2];
		e=ptint[i][3];
		theta=pi*(double)c/180.0;
		phi=pi*(double)d/180.0;
		ramda=pi*(double)e/180.0;
		for(j = 0; j < ptotal;j++){
			x0=pos2v[j][0]*cos(ramda)+pos2v[j][1]*sin(ramda);
			y0=-1.0*pos2v[j][0]*sin(ramda)+pos2v[j][1]*cos(ramda);
			z0=pos2v[j][2]-radii+ptr[i][0];
			x1=x0;
			y1=y0*cos(phi)+z0*sin(phi);
			z1=-1.0*y0*sin(phi)+z0*cos(phi);
			x2=x1*cos(theta)+y1*sin(theta);
			y2=-1.0*x1*sin(theta)+y1*cos(theta);
			z2=z1;
			x=x2+paori[b][0]-ptr[i][1];
			y=y2+paori[b][1]-ptr[i][2];
			z=z2+paori[b][2]-ptr[i][3];
			if (x<0){x+=axis[0];}
			if (y<0){y+=axis[1];}
			if (z<0){z+=axis[2];}
			if (x>axis[0]){x-=axis[0];}
			if (y>axis[1]){y-=axis[1];}
			if (z>axis[2]){z-=axis[2];}

			fprintf(fp3,"%d %d %.3f %.3f %.3f\n",a,50+i,x,y,z);
			a+=1;
		}
	}
		
	}
	
	fclose(fp3);
	
	printf("Checking the generated structure\n");
	for(i = 0; i < ptnum;i++){
		for(j = 0; j < 3;j++){ptmid[i][j]=0.0;}
		a=ptint[i][0];
		c=ptint[i][1];
		d=ptint[i][2];
		e=ptint[i][3];
		theta=pi*(double)c/180.0;
		phi=pi*(double)d/180.0;
		ramda=pi*(double)e/180.0;
		for(j = 0; j < ptotal;j++){
			x0=pos2v[j][0]*cos(ramda)+pos2v[j][1]*sin(ramda);
			y0=-1.0*pos2v[j][0]*sin(ramda)+pos2v[j][1]*cos(ramda);
			z0=pos2v[j][2]-radii+ptr[i][0];
			x1=x0;
			y1=y0*cos(phi)+z0*sin(phi);
			z1=-1.0*y0*sin(phi)+z0*cos(phi);
			x2=x1*cos(theta)+y1*sin(theta);
			y2=-1.0*x1*sin(theta)+y1*cos(theta);
			z2=z1;
			x=x2+paori[a][0]-ptr[i][1];
			y=y2+paori[a][1]-ptr[i][2];
			z=z2+paori[a][2]-ptr[i][3];
			ptmid[i][0]+=x;
			ptmid[i][1]+=y;
			ptmid[i][2]+=z;
		}
		for(j = 0; j < 3;j++){ptmid[i][j]/=(double)ptotal;}
//		printf("check %d %f %f %f\n",i+50,ptmid[i][0],ptmid[i][1],ptmid[i][2]);
	}

	for(i = 0; i < ptnum-1;i++){
		for(j = i+1 ; j < ptnum;j++){
			dx=fabs(ptmid[i][0]-ptmid[j][0]);
			dy=fabs(ptmid[i][1]-ptmid[j][1]);
			dz=fabs(ptmid[i][2]-ptmid[j][2]);
			if (dx>0.5*axis[0]){dx-=axis[0];}
			if (dy>0.5*axis[1]){dy-=axis[1];}
			if (dz>0.5*axis[2]){dz-=axis[2];}
			r=sqrt(dx*dx+dy*dy+dz*dz);
			if (r < cutoff*2){
				printf("Distance %d(%d) %d(%d) %.2f\n",i+50,ptint[i][0],j+50,ptint[j][0],r);
			}
		}
	}
    
    
	fp4 = fopen("ptoncb.rd","w");
	i=0; sw=0;
	fprintf(fp4,"#cellx  0.000  %7.3lf\n",axis[0]);
	fprintf(fp4,"#celly  0.000  %7.3lf\n",axis[1]);
	fprintf(fp4,"#cellz  0.000  %7.3lf\n",axis[2]);
	fprintf(fp4,"\n#masses 2\n");
	fprintf(fp4,"  1   12.011\n");
	fprintf(fp4,"  4  195.078\n");
	fprintf(fp4,"\n#thermofree 5\n");
	fprintf(fp4,"\n#atoms %d\n",total);
	for(i = 0; i < ctotal;i++){
		fprintf(fp4,"%8d  1   %d",i+1, mask[i]);
            fprintf(fp4,"   %7.3lf   %7.3lf   %7.3lf",pos1v[i][0],pos1v[i][1],pos1v[i][2]);
            fprintf(fp4,"   %7.3lf   %7.3lf   %7.3lf\n",pos1v[i][3]/100,pos1v[i][4]/100,pos1v[i][5]/100);
	}
	k=1;
	for(i = 0; i < ptnum;i++){
		a=ptint[i][0];
		c=ptint[i][1];
		d=ptint[i][2];
		e=ptint[i][3];
		theta=pi*(double)c/180.0;
		phi=pi*(double)d/180.0;
		ramda=pi*(double)e/180.0;
		x1=0.0;
		y1=0.0*cos(phi)+velocity*sin(phi);
		z1=-1.0*0.0*sin(phi)+velocity*cos(phi);
		x2=x1*cos(theta)+y1*sin(theta);
		y2=-1.0*x1*sin(theta)+y1*cos(theta);
		z2=z1;
		vx=-1.0*x2;
		vy=-1.0*y2;
		vz=-1.0*z2;
		for(j = 0; j < ptotal;j++){
			x0=pos2v[j][0]*cos(ramda)+pos2v[j][1]*sin(ramda);
			y0=-1.0*pos2v[j][0]*sin(ramda)+pos2v[j][1]*cos(ramda);
			z0=pos2v[j][2]-radii+ptr[i][0];
			x1=x0;
			y1=y0*cos(phi)+z0*sin(phi);
			z1=-1.0*y0*sin(phi)+z0*cos(phi);
			x2=x1*cos(theta)+y1*sin(theta);
			y2=-1.0*x1*sin(theta)+y1*cos(theta);
			z2=z1;
			x=x2+paori[a][0]-ptr[i][1];
			y=y2+paori[a][1]-ptr[i][2];
			z=z2+paori[a][2]-ptr[i][3];
			if (x<0){x+=axis[0];}
			if (y<0){y+=axis[1];}
			if (z<0){z+=axis[2];}
			if (x>axis[0]){x-=axis[0];}
			if (y>axis[1]){y-=axis[1];}
			if (z>axis[2]){z-=axis[2];}
			fprintf(fp4,"%8d  4   10",ctotal+k);
            fprintf(fp4,"   %7.3lf   %7.3lf   %7.3lf",x,y,z);
            fprintf(fp4,"   %7.6lf   %7.6lf   %7.6lf\n",vx/100,vy/100,vz/100);
			k+=1;
		}
	}
	fclose(fp4);
    double Ptmass,Cmass;
    Ptmass = ptnum*Ptclustnum*195.078/AVOGADRO;
    Cmass = 400000*12.001/AVOGADRO*particle_num;
    printf("Pt/C wt_ratio = %lf\n",Ptmass/Cmass);
	printf("ptoncb.rd is generated !\n");
}
