#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <iomanip>
#include <cstdlib>
#include <math.h>
using namespace std;

#define pi 3.14159265359
#define particle_num 2
#define BUF 256

//-----------------------Parameter & Infomation-----------------------------

string input = "sorteddump.pos.0";
int initial = 451918;//nafionの行の始まり
int finish = 951089;//nafionの行の終わり //原子の最後に値する
int nafion_atom = 416; //nafionの原子数

//-----------------------Parameter & Infomation-----------------------------

double center[particle_num][3] ={};

class Atom{
public:
    int id_ori;
    int id,type;
    double x,y,z;
    double vx,vy,vz;
    double q;
};

double mass[6] = {12.011,1.008,15.999,195.078,32.060,18.988};
                  /*  C      H     O      Pt      S      F  */


////////////////////////////
void Read_dump(vector<Atom> &atoms);
void read_newcenter();
////////////////////////////

int main(){
    vector<Atom> atoms;
    vector<Atom> nafion_g;
    Read_dump(atoms);
    vector<Atom> nafion_tip;
    vector<Atom> cog2center;
    vector<Atom> cog2tip;
    //nafionの本数
    int nafion_number = (finish - initial + 1) / nafion_atom;

    //centerを読み込む
    read_newcenter();
    for(int i=0; i<particle_num; i++)printf("#center x=%lf y=%lf z=%lf\n", center[i][0], center[i][1], center[i][2]);

    //nafionの配列を変更(nafionの本数と同じだけにする)
    nafion_g.resize(nafion_number);
    
    //nafion一本の質量を求める
    double nafion_mass = 0; 
    for(int i = 0; i < nafion_atom; i ++){
        if(atoms[initial + i].type == 1){nafion_mass += 12.011;}  // C
        if(atoms[initial + i].type == 2){nafion_mass += 1.008;}   // H
        if(atoms[initial + i].type == 3){nafion_mass += 15.999;}  // O
        if(atoms[initial + i].type == 4){nafion_mass += 195.078;} // Pt
        if(atoms[initial + i].type == 5){nafion_mass += 32.060;}  // S
        if(atoms[initial + i].type == 6){nafion_mass += 18.988;}  // F
    }

    //nafionの重心を求める
    for(int i =0; i < nafion_number; i++){
        for(int j =0; j < nafion_atom; j++){
            if(atoms[(initial + i*nafion_atom) + j].type == 1 ){
                nafion_g[i].x += mass[0] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[0] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[0] *atoms[(initial + i * nafion_atom) + j].z;
            }
            else if(atoms[(initial + i*nafion_atom) + j].type == 2 ){
                nafion_g[i].x += mass[1] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[1] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[1] *atoms[(initial + i * nafion_atom) + j].z;
            }

            else if(atoms[(initial + i*nafion_atom) + j].type == 3 ){
                nafion_g[i].x += mass[2] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[2] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[2] *atoms[(initial + i * nafion_atom) + j].z;
            }
            
            else if(atoms[(initial + i*nafion_atom) + j].type == 4 ){
                nafion_g[i].x += mass[3] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[3] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[3] *atoms[(initial + i * nafion_atom) + j].z;
            }
            
            else if(atoms[(initial + i*nafion_atom) + j].type == 5 ){
                nafion_g[i].x += mass[4] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[4] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[4] *atoms[(initial + i * nafion_atom) + j].z;
            }

            else if(atoms[(initial + i*nafion_atom) + j].type == 6 ){
                nafion_g[i].x += mass[5] *atoms[(initial + i * nafion_atom) + j].x;
                nafion_g[i].y += mass[5] *atoms[(initial + i * nafion_atom) + j].y;
                nafion_g[i].z += mass[5] *atoms[(initial + i * nafion_atom) + j].z;
            }
        }
    }

    double cog_nafion[nafion_g.size()][3];
    double sum[nafion_g.size()][3];
    double dis_from_center[particle_num]={};
    double xyz[particle_num][3]={};
    int flag[nafion_g.size()]={};
    //nafion_gの合計を求める
    for(int i = 0; i < nafion_g.size(); i++){
        sum[i][0] += nafion_g[i].x;
        sum[i][1] += nafion_g[i].y;
        sum[i][2] += nafion_g[i].z;
        
        //分母÷分子
        cog_nafion[i][0] = sum[i][0] / nafion_mass;
        cog_nafion[i][1] = sum[i][1] / nafion_mass;
        cog_nafion[i][2] = sum[i][2] / nafion_mass;

        //最も近いcenterのflagを求める
        //まずは，中心とNafionの重心の距離を求める．

        for(int j=0;j<particle_num; j++){
            for(int k=0;k<3;k++){
                xyz[j][k]+=cog_nafion[i][k]-center[j][k];
            }
            dis_from_center[j]=xyz[j][0]*xyz[j][0]+xyz[j][1]*xyz[j][1]+xyz[j][2]*xyz[j][2];
        }
        
        int min_dis=10000000;//minimum distance
        for(int j=0;j<particle_num; j++){
            if(min_dis > dis_from_center[j]){
            min_dis = dis_from_center[j];
            flag[i] = j;
            }
        }

    }

 //cout << "nafionの質量:" << nafion_mass << endl;
 //for(int i = 0; i < nafion_g.size(); i ++){
 //    cout << "nafionの重心" << "("<< i+1 << "):" << cog_nafion[i][0] <<","<<  cog_nafion[i][1] <<"," << cog_nafion[i][2]<< endl;} 

 //nafionの先端の座標
 nafion_tip.resize(nafion_number);
 for(int i = 0; i < nafion_number; i++){ 
     nafion_tip[i].x = atoms[(initial + i * nafion_atom)].x;
     nafion_tip[i].y = atoms[(initial + i * nafion_atom)].y;
     nafion_tip[i].z = atoms[(initial + i * nafion_atom)].z;
  //   cout << "nafionの先端の座標" << "("<< i+1 << "):" << nafion_tip[i].x <<","<< nafion_tip[i].y <<"," << nafion_tip[i].z<< endl;
 } 

 //gravity of center と centerで作るベクトル　cog to center
 cog2center.resize(nafion_number);
 int fl;
 for(int i =0;i < nafion_number; i++){
    fl = flag[i];
    cog2center[i].x = center[fl][0] - cog_nafion[i][0];
    cog2center[i].y = center[fl][1] - cog_nafion[i][1];
    cog2center[i].z = center[fl][2] - cog_nafion[i][2];

  //   cout << "nafionの重心と系の中心で作るベクトル" << "("<< i+1 << "):" << cog2center[i].x <<","<< cog2center[i].y <<"," << cog2center[i].z<< endl;
 }

 //gravity of center とnafion_tipで作るベクトル　cog2tip
 cog2tip.resize(nafion_number);
 for(int i =0;i < nafion_number; i++){
    cog2tip[i].x = nafion_tip[i].x - cog_nafion[i][0];
    cog2tip[i].y = nafion_tip[i].y - cog_nafion[i][1];
    cog2tip[i].z = nafion_tip[i].z - cog_nafion[i][2];
  //  cout << "nafionの先端と重心で作るベクトル" << "("<< i+1 << "):" << cog2tip[i].x <<","<< cog2tip[i].y <<"," << cog2tip[i].z<< endl;
 }


 //ベクトルから角度を求める
 double costheta[nafion_number];
 double arccos[nafion_number];
 double part1[nafion_number];
 double part2[nafion_number];
 double part3[nafion_number];
 for(int i =0; i < nafion_number; i++){
      part1[i] = cog2center[i].x * cog2tip[i].x + cog2center[i].y * cog2tip[i].y + cog2center[i].z * cog2tip[i].z; 
      part2[i] = sqrt(pow(cog2center[i].x,2) +pow(cog2center[i].y,2) + pow(cog2center[i].z,2));
      part3[i] = sqrt(pow(cog2tip[i].x,2) +pow(cog2tip[i].y,2) + pow(cog2tip[i].z,2));
      costheta[i] = part1[i]/(part2[i] * part3[i]);
      arccos[i] = acos(costheta[i]);
 //     cout << "cos(θ)" << "("<< i+1 << "):" << costheta[i]<< endl;
 }


 //中心からの距離を求める
 double distance[nafion_number];
 for(int i = 0; i < nafion_number; i++){
     distance[i] = sqrt(pow(cog2center[i].x,2) + pow(cog2center[i].y,2) + pow(cog2center[i].z,2));
 }
 //nafionの数と角度を表示する
 //中心からの距離を記載する
 
 /*
 cout <<"nafionの数:" <<  nafion_g.size() << endl;
 for(int i = 0; i < nafion_number; i++){
     cout << "θ(degree)" << "("<< i+1 << "):" << arccos[i] * 360 /(2*pi)<< endl;
     
 }
 */

 double angle[nafion_number];
 for(int i = 0;i < nafion_number; i++){
     angle[i] = arccos[i] * 360 /(2*pi);
     if(angle[i] >= 90){
         angle[i] = angle[i] - 90;
     }
 }

 for(int i = 0; i < nafion_number;i++){
     cout << distance[i] << "  " << angle[i] << endl;
 }


 /*Hermanの配向度
 double Herman[nafion_number];
 for(int i = 0; i <nafion_number; i ++){
     Herman[i] =1.5*pow(costheta[i],2)-0.5;
     cout << Herman[i] << endl;
 }
 */

}


void Read_dump(vector<Atom> &atoms){
    //cout << "Set sorted dump file:";
    //cin >> input;
    ifstream ifs(input.c_str());
    if(!ifs){
        cout << "Error: No " << input << endl;
        exit(1);
    }
    string buff;
    while(getline(ifs,buff)){
        stringstream ss;
        ss << buff;
        string tmp;
        ss >> tmp; //first word
        ss >> tmp; //second word
        if(tmp =="NUMBER"){
            getline(ifs,buff);
            int numatoms = std::stoi(buff);
            atoms.resize(numatoms);
        }
        if(tmp =="ATOMS"){
            int numatoms =atoms.size();
            for(int i = 0; i < numatoms; i++){
                getline(ifs,buff);
                stringstream sss;
                sss << buff;
                int id, type;
                double x,y,z,vx,vy,vz,q;
                sss >> id >> type >> x >> y >> z >> vx >> vy >> vz >> q;
                atoms[i].id_ori = id;
                atoms[i].id = i+1;
                atoms[i].type = type;
                atoms[i].x = x; atoms[i].y = y; atoms[i].z = z;
                atoms[i].vx = vx; atoms[i].vy = vy; atoms[i].vz = vz;
            }
        }
    }
   // cout << "---" << input << "  was read !" << endl;
}

void read_newcenter(){

	int i, x;
    double a, b, c;
    FILE *fp;
    char buf[BUF];
        
    fp = fopen("new_center.txt","r");  //center
    if( fp==NULL ){
   	    printf("Error : These is no %s !! This process will stop.\n","new_center.txt");
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

