#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUF 256
#define cutoff 9.0

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

int main(int argc,char **argv)
{
	int i,j,k,s,a,b,c,d,e,f,g,sw,flag,cont;
	int total,ctotal,ptotal,panum,patom;
	int in_total,out_total,in_num,out_num;
	double	(*pos1v)[10],(*pos2v)[6],maxis[3],daxis[3],axis[3];
	double	pi,x,y,z,r,vx,vy,vz,fx,fy,fz,q,*tempv,vol,bo;
	double	x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dx,dy,dz;
	int	(*atom1c)[3],*index;
	char buf1[BUF],buf2[BUF],file1[BUF],file2[BUF],file3[BUF],*temp;
	char command[256],ele[2],tt[8];
	FILE *fp1,*fp2,*fp3,*fp4;
	int coor[111]={1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 48, 49, 50, 65, 66, 86, 87, 88, 90, 102, 103, 104, 106, 107, 108, 110, 112, 114, 115, 116, 117, 120, 121, 122, 129, 130, 142, 143, 159, 160, 180, 181, 182, 184, 186, 188, 190, 192, 194, 196, 197, 198, 202, 206, 207, 215, 216, 228, 229, 245, 246, 248, 250, 252, 254, 256, 258, 260, 262, 268, 272, 273, 281, 282, 284, 294, 295, 296, 298, 300, 302, 304, 310, 314, 316, 320, 322};
	
    pi=4*atan(1);
    flag=0;
	s=0;
	for(k = 20; k <= 20;k++){
                
	printf("%d %d ",s,k*10000);
	sprintf(file1,"new_dump.pos.%d",k*10000);
	sprintf(file2,"distri_dump.pos.%d",k*10000);
	
	fp1 = fopen(file1,"r");/*dump file*/
	if( fp1==NULL ){
		printf("Error : These is no %s !! This process will stop.\n",file1);
		fclose(fp1);
		exit(1);
	}
	fp2 = fopen(file2,"w");/*dump file*/
	for(i = 0; i < 9;i++){
		fgets(buf1,BUF,fp1);
		if(i==3)sscanf(buf1,"%d\n",&total);
		if(i>4)sscanf(buf1,"%lf%lf\n",&daxis[i-5], &maxis[i-5]);
		if(i>4)axis[i-5]=maxis[i-5]-daxis[i-5];
	}
		
	fprintf(fp2,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\n",k*10000,total);
	fprintf(fp2,"ITEM: BOX BOUNDS pp pp pp\n");
	fprintf(fp2,"0.0 %.6lf\n",axis[0]);
	fprintf(fp2,"0.0 %.6lf\n",axis[1]);
	fprintf(fp2,"0.0 %.6lf\n",axis[2]);
	fprintf(fp2,"ITEM: ATOMS id type x y z coverage\n");		
	
	atom1c=malloc(sizeof(int)*total*3);
	pos1v=malloc(sizeof(double)*total*10);
    index=malloc(sizeof(int)*total);
	tempv=malloc(sizeof(double)*total);
	cont=0;
	for(i = 0; i <total;i++){
        fgets(buf1,BUF,fp1);
		sscanf(buf1,"%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&a,&b,&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&q);
		index[a-1]=a;
		atom1c[a-1][0]=b;
		atom1c[a-1][1]=b;
		atom1c[a-1][2]=0;
		if (b==4) {
			for(j = 0; j <111;j++){
				if (cont==coor[j])atom1c[a-1][2]=1;
			}
			cont++;
		}
			
		if (cont==326) cont=0;
		if (b==21) atom1c[a-1][0]=1;
		if (b==22) atom1c[a-1][0]=1;
		pos1v[a-1][0]=x;
		pos1v[a-1][1]=y;
		pos1v[a-1][2]=z;
		pos1v[a-1][3]=vx;
		pos1v[a-1][4]=vy;
		pos1v[a-1][5]=vx;
		pos1v[a-1][6]=fx;
		pos1v[a-1][7]=fy;
		pos1v[a-1][8]=0.0;
		pos1v[a-1][9]=0.0;
	}
	
	
	
	in_num=0;
	in_total=0;
	out_num=0;
	out_total=0;
	c=4; //C 1 Pt 4
    d=6; //F 6 SO3H 5
	g=0; //Start of C
	e=289934 ;//Start of Pt
	f=294498 ;//Start of Nafion
	for(i = e; i < f;i++){
//		if (i%10000==0) printf("%d\n",i);
		if (atom1c[i][0]==c){
		for(j = f; j<total;j++){
			dz=pos1v[j][2]-pos1v[i][2];
			dx=pos1v[j][0]-pos1v[i][0];
			dy=pos1v[j][1]-pos1v[i][1];
			r=sqrtl(dx*dx+dy*dy+dz*dz);
			if ((r<cutoff)&&(atom1c[j][0]==d)) pos1v[i][9]-=1.0;
		}
			
		dx=pos1v[i][0]-0.5*axis[0];
		dy=pos1v[i][1]-0.5*axis[1];
		dz=pos1v[i][2]-0.5*axis[2];
		r=sqrtl(dx*dx+dy*dy+dz*dz);
		if ((r<101.0)&&(atom1c[i][2]==1)&&(pos1v[i][9]<0)) in_num+=1;
		if ((r<101.0)&&(atom1c[i][2]==1)) in_total+=1;
		if ((r>=101.0)&&(atom1c[i][2]==1)&&(pos1v[i][9]<0)) out_num+=1;
		if ((r>=101.0)&&(atom1c[i][2]==1)) out_total+=1;
		}
		
	}
	
	printf("%d %d %d %d ",in_num, in_total,out_num,out_total);
	printf("Time(ps) IN/OUT %.2f ",k*10*0.25);
	printf("%.3lf  %.3lf \n",(double)in_num/(double)in_total,(double)out_num/(double)out_total);

	
	for(i = 0; i <total;i++){
		a=atom1c[i][1];
		fprintf(fp2,"%d %d ",i+1,a);
		fprintf(fp2,"%.3f %.3f %.6f %.1f\n",pos1v[i][0],pos1v[i][1],pos1v[i][2],pos1v[i][9]);
	}
	
//	printf("%s is generated !\n",file2);
        s+=1;
	fclose(fp1);
	fclose(fp2);
	}
}
