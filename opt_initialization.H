/********************************************************************************/
/*****************必要变量声明及初始化,如MMA***********************************/
/********************************************************************************/


Info<< "Clear output files of last solve \n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("Compliance.txt");

int i = 0, opt = 0;                     // opt用于优化过程参数的放松，在"update.H"里opt++,这里先进行初始化
double del = 0;                       // del就是heviside过滤里的惩罚强度参数betta
scalar n_dof = xh.size();               //变量个数，不考虑并行
scalar m_con = 1;                           // 约束个数
scalar eta0, eta5, eta1, y0, y5;//eta5就是heviside过滤的截断值eta_h，小于该值的伪密度会倾向于罚为0；大于该值的倾向于罚为1

// n = xh.size(); // n:length of x in the i-th core
// N = xh.size(); // N:design freedom;
// reduce(N, sumOp<scalar>());// 把所有核中的值累加在一起，reduce()函数中sumOp描述累加，还有minOp，maxOp等

scalar C(0), V(0), area = gSum(mesh.V());                            //定义了柔度C，体积V，area是单元体积和，gSum() sums over all the processors in a parallel run
double len = std::pow(area / N, 1.0 / 3);                            // 这个len计算出所有单元的平均尺寸
b.value() = 1.0 / (filterR * len / 3.464) / (filterR * len / 3.464); // b就是PDE求导数中的1/( r^2 ),r=filterR/(2*3^1/3) in paper,3.464=2*3^1/3

/*****************************MMA相关*******************************************/

Info<< "Initialization of MMA \n" << endl;
//------------------------------------------------------------------------------------------------------------//