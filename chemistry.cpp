#include<bits/stdc++.h>
using namespace std;
typedef long long Int;
const Int MAXSIZE = 200;
inline Int gcd(Int x,Int y){return x%y==0?y:gcd(y,x%y);}
struct frac{
	Int a,b;//a/b
	void reduce(){Int x=gcd(a,b);a/=x,b/=x;};
	frac operator = (Int x){a=x,b=1;return *this;};
	frac operator = (const frac x){a=x.a,b=x.b;reduce();return *this;};
	frac operator + (const frac x){return (frac){b*x.a+a*x.b,b*x.b};};
	frac operator - (const frac x){return (frac){a*x.b-b*x.a,b*x.b};};
	frac operator * (const frac x){return (frac){a*x.a,b*x.b};};
	frac operator / (const frac x){return (frac){a*x.b,b*x.a};};
	bool operator < (const frac x){return a*x.b<b*x.a;};
	bool operator == (const frac x){return a*x.b==b*x.a;};
};//分数结构体
struct matrix{Int Size_X,Size_Y;vector< vector<frac> > Matrix;} r;//矩阵
Int i,j,k;frac X[MAXSIZE]/*答案*/;
string chemistry/*化学方程式*/,Lchemistry/*等号左边的化学方程式*/,Rchemistry/*等号右边的化学方程式*/;
string lch[MAXSIZE],rch[MAXSIZE];Int lenl=1,lenr=1;//等号左右边物质集合
map<string,Int>chemToline;//每种元素对应矩阵的哪一行
string nxt[MAXSIZE];//(H)*+(OH)-=H2O   As2S3+H2O+(NO3)-=(AsO4)3-+(SO4)2-+NO+(H)*
Int foundNot0(vector<frac>h){ for(Int d=1;d<h.size();d++){ if(h[d].a!=0)return d; }return 998244353; }
//寻找第一个不为零的矩阵项
void processCheml(string chem,Int ret,Int doubl){
	Int chemLength = chem.length();
	if(chem[chemLength-1]=='-'||chem[chemLength-1]=='*'){//电荷守恒
		if(chemToline.count("e")==0){
			Int nxtLine=chemToline.size()+1;
			chemToline["e"]=nxtLine;
		}
		int endnum=chemLength;
		for(Int d=chemLength-1;d>=0;d--)
			if(chem[d]==')'){endnum = d+1;break;}
		int tot = 0;
		for(Int d=endnum;d<chemLength-1;d++)
			tot = tot*10+ ( chem[d]-'0' );
		if(tot==0)tot=1;
		if(chem[chemLength-1]=='-')
			r.Matrix[chemToline["e"]][ret]= r.Matrix[chemToline["e"]][ret] + (frac){-tot,1};
		else r.Matrix[chemToline["e"]][ret]= r.Matrix[chemToline["e"]][ret] + (frac){tot,1};
		string tmp="";
		for(Int d=0;d<endnum;d++)tmp+=chem[d];
		chem=tmp;
	}
	/*元素守恒*/Int d=0;Int nxtnum=0,doub[100];
	while(d<chemLength){
		if( (chem[d]>='A'&&chem[d]<='Z') && 
			(d+1<chem.length()&&chem[d+1]>='a'&&chem[d+1]<='z') ) {
			//有一大一小两个字母的元素
			string chemName;chemName+=chem[d];chemName+=chem[d+1];
			Int e,tot = 0;
			for(e=d+2;e<chemLength;e++){
				if(!(chem[e]>='0'&&chem[e]<='9') )break;
				tot=tot*10+(chem[e]-'0');
			}d = e-1;
			if(chemToline.count(chemName)==0){
				Int nxtLine=chemToline.size()+1;
				chemToline[chemName]=nxtLine;
			}
			if(tot==0)tot=1;
			r.Matrix[chemToline[chemName]][ret]= r.Matrix[chemToline[chemName]][ret] + (frac){tot,1};
		}
		else if(chem[d]>='A'&&chem[d]<='Z') {
			//只有一个大写字母的元素
			string ch;ch+=chem[d];
			Int e,tot = 0;
			if(d+1==chem.length())tot=1;
			else{
				for(e=d+1;e<chem.length();e++){
					if(!(chem[e]>='0'&&chem[e]<='9') )
						break;
					tot=tot*10+(chem[e]-'0');
				}d=e-1;
			}
			if(chemToline.count(ch)==0){
				Int si=chemToline.size()+1;
				chemToline[ch]=si;
			}
			if(tot==0)tot=1;
			r.Matrix[chemToline[ch]][ret]= r.Matrix[chemToline[ch]][ret] + (frac){tot,1};
		}
		else if(chem[d]=='(') {
			//遇到括号，将括号内取出
			Int e;nxtnum++;
			for(e=d+1;e<chemLength;e++){
				if(chem[e]==')')
					break;
				nxt[nxtnum]+=chem[e];
			}
			Int f,tot=0;
			for(f=e+1;f<chemLength;f++){
				if(!(chem[f]>='0'&&chem[f]<='9') )
					break;
				tot=tot*10+(chem[e]-'0');
			}
			if(tot==0)tot=1;
			doub[nxtnum]=tot;
			d = f-1;
		}
		d++;
	}
	for(Int ii=1;ii<=nxtnum;ii++)processCheml(nxt[ii],ret,doubl*doub[ii]);//继续统计
}
void processChemr(string chem,Int ret,Int doubl){
	Int chemLength = chem.length();
	if(chem[chemLength-1]=='-'||chem[chemLength-1]=='*'){//电荷守恒
		if(chemToline.count("e")==0){
			Int nxtLine=chemToline.size()+1;
			chemToline["e"]=nxtLine;
		}
		int endnum=chemLength;
		for(Int d=chemLength-1;d>=0;d--)
			if(chem[d]==')'){endnum = d+1;break;}
		int tot = 0;
		for(Int d=endnum;d<chemLength-1;d++)
			tot = tot*10+ ( chem[d]-'0' );
		if(tot==0)tot=1;
		if(chem[chemLength-1]=='-')
			r.Matrix[chemToline["e"]][ret]= r.Matrix[chemToline["e"]][ret] + (frac){tot,1};
		else r.Matrix[chemToline["e"]][ret]= r.Matrix[chemToline["e"]][ret] + (frac){-tot,1};
		string tmp="";
		for(Int d=0;d<endnum;d++)tmp+=chem[d];
		chem=tmp;
	}/*元素守恒*/Int d=0,nxtnum=0,doub[100]={0};
	while(d<chemLength){
		if( (chem[d]>='A'&&chem[d]<='Z') && 
			(d+1<chem.length()&&chem[d+1]>='a'&&chem[d+1]<='z') ) {
			//有一大一小两个字母的元素
			string ch;ch+=chem[d];ch+=chem[d+1];
			Int e = 0,tot = 0;
			for(e=d+2;e<chem.length();e++){
				if(!(chem[e]>='0'&&chem[e]<='9') )
					break;
				tot=tot*10+(chem[e]-'0');
			}d=e-1;
			if(chemToline.count(ch)==0){
				Int si=chemToline.size()+1;
				chemToline[ch]=si;
			}
			if(tot==0)tot=1;
			r.Matrix[chemToline[ch]][ret] = r.Matrix[chemToline[ch]][ret] - (frac){tot*doubl,1};
		}
		else if(chem[d]>='A'&&chem[d]<='Z') {
			//只有一个大写字母的元素
			string chemName;chemName+=chem[d];
			Int e,tot = 0;
			if(d+1 == chem.length())
				tot=1;
			else{
				for(e=d+1;e<chemLength;e++){
					if(!(chem[e]>='0'&&chem[e]<='9') )
						break;
					tot=tot*10+(chem[e]-'0');
				}d=e-1;
			}
			if(chemToline.count(chemName)==0){Int nxtLine = chemToline.size()+1;chemToline[chemName] = nxtLine;}
			if(tot==0) tot=1;
			r.Matrix[chemToline[chemName]][ret] = r.Matrix[chemToline[chemName]][ret] - (frac){tot*doubl,1};
		}
		else if(chem[d]=='(') {
			//遇到括号，将括号内取出
			Int e;nxtnum++;
			for(e=d+1;e<chemLength;e++){
				if(chem[e]==')')break;
				nxt[nxtnum]+=chem[e];
			}
			Int f,tot=0;
			for(f=e+1;f<chemLength;f++){
				if(!(chem[f]>='0'&&chem[f]<='9') )break;
				tot=tot*10+(chem[f]-'0');
			}
			if(tot==0)tot=1;
			doub[nxtnum]=tot;
			d = f-1;
		}
		d++;
	}
	for(Int ii=1;ii<=nxtnum;ii++)
		processChemr(nxt[ii],ret,doubl*doub[ii]);
}
signed main(){
	getline(cin,chemistry);
	bool found_same = false;
	for(i=0;i<chemistry.size();i++){
		if(chemistry[i]=='='){found_same = true;continue;}
		if(!found_same)Lchemistry+=chemistry[i];
		else Rchemistry+=chemistry[i];
	}
	for(i=0;i<Lchemistry.size();i++){if(Lchemistry[i]=='+'){ lenl++;continue; }lch[lenl]+=Lchemistry[i];}
	for(i=0;i<Rchemistry.size();i++){if(Rchemistry[i]=='+'){ lenr++;continue; }rch[lenr]+=Rchemistry[i];}
	vector<frac>emp;r.Matrix.push_back(emp);
	for(i=1;i<=100;i++){
		r.Matrix.push_back(emp); r.Matrix[i].push_back( (frac){0,1} );
		for(j=1;j<=100;j++) r.Matrix[i].push_back( (frac){0,1} );
	}processChemr(lch[1],lenl+lenr,1);//第一个特殊处理
	for(i=2;i<=lenl;i++)fill(nxt,nxt+MAXSIZE-5,""),processCheml(lch[i],i-1,1);
	for(i=1;i<=lenr;i++)fill(nxt,nxt+MAXSIZE-5,""),processChemr(rch[i],lenl+i-1,1);
	r.Size_X = chemToline.size();
	r.Size_Y=lenl+lenr;
	for(i=1;i<=r.Size_X;i++){
		//手动选择排序
		for(j=i;j<=r.Size_X;j++)
			for(k=j+1;k<=r.Size_X;k++){
				bool canSwap = false,front0=true;
				for(Int f=i;f<=r.Size_Y;f++){
					if(r.Matrix[j][f].a!=0) front0=false;
					//判断当前矩阵元素是否为前导零的一部分，如果是那么交换
					if(r.Matrix[j][f]<r.Matrix[k][f]||front0){ canSwap = true;break;}//如果可交换，那么交换
					else if(r.Matrix[k][f]<r.Matrix[j][f]) break;//如果后边比前边小，放弃交换并退出
				}
				if(canSwap==true) for(Int f=1;f<=r.Size_Y;f++) swap(r.Matrix[j][f],r.Matrix[k][f]);//交换
			}
		Int not0 = foundNot0(r.Matrix[i]);
		if(not0>r.Size_Y)break;
		if( !( r.Matrix[i][not0] == frac{1,1} ) ) {
			frac front = r.Matrix[i][not0];
			for(j=not0;j<=r.Size_Y;j++)
				r.Matrix[i][j]=r.Matrix[i][j] / front,r.Matrix[i][j].reduce();
		}
		for(j=i+1;j<=r.Size_X;j++) {
			frac twice = r.Matrix[j][not0] / r.Matrix[i][not0];
			for(k=not0;k<=r.Size_Y;k++)
				r.Matrix[j][k] = r.Matrix[j][k] - r.Matrix[i][k]*twice,r.Matrix[j][k].reduce();
		}
	}
	for(j=1;j<=r.Size_X;j++)
		for(k=j+1;k<=r.Size_X;k++){
			Int not0_j = foundNot0(r.Matrix[j]),not0_k=foundNot0(r.Matrix[k]);
			if(not0_j>not0_k) for(Int f=1;f<=r.Size_Y;f++) swap(r.Matrix[j][f],r.Matrix[k][f]);//交换
		}
	Int nowx=r.Size_Y-1;
	for(i=r.Size_X;i>=1;i--){
		if(r.Matrix[i][nowx].a==0)continue;X[nowx]=r.Matrix[i][r.Size_Y];
		for(j=nowx+1;j<=r.Size_Y-1;j++){ X[nowx] = X[nowx] - (r.Matrix[i][j] * X[j]);}nowx--;
	}X[0]=(frac){1,1};Int minidouble = X[0].b;
	for(i=1;i<r.Size_Y;i++)minidouble=X[i].b * minidouble/( gcd(X[i].b, minidouble) );
	for(i=0;i<r.Size_Y;i++) X[i] = X[i] * (frac){minidouble,1};
	for(i=0;i<lenl;i++){
		if(X[i].a!=1) cout<<X[i].a;
		if(i<lenl-1)cout<<lch[i+1]<<" + ";
		else cout<<lch[i+1];
	}cout<<" = ";
	for(i=lenl;i<lenl+lenr;i++){
		if(X[i].a!=1) cout<<X[i].a;
		if(i<lenl+lenr-1)cout<<rch[i-lenl+1]<<" + ";
		else cout<<rch[i-lenl+1];
	}cout<<endl;system("pause");return 0;
}
/*以下为hack数据
H2O2+H2CrO4=O2+Cr(OH)3+H2O
KMnO4+H2O2+H2SO4=K2SO4+MnSO4+O2+H2O
*/
