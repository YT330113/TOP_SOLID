/**************************************************************************/
/*****************目标函数和约束函数值计算*******************************/
/**************************************************************************/


Info <<"\nCaculate cost function\n"<<endl;
/****************************计算柔度和体积约束******************************/

gradD = fvc::grad(D);
sigmaD = mu * twoSymm(gradD) + (lambda * I) * tr(gradD);
volScalarField Energy(0.25 * (gradD + gradD.T()) && sigmaD);

C = fvc::domainIntegrate(Energy).value();//f(x)
V = fvc::domainIntegrate(xh).value() / area - voluse; // g(x)<=0;     area和voluse都是无量刚标量
//-----------------------------------------------------------------------------------------------------------//

/****************************输出每一迭代步的目标函数值*********************/

if (Pstream::master())
{
    ofstream outfile1("Voluse.txt", std::ios::app);
    outfile1 << V << "\n";
    outfile1.close();

    ofstream outfile2("Compliance.txt", std::ios::app);
    outfile2 << C << "\n";
    outfile2.close();

    ofstream outfile3("Time.txt", std::ios::app);
    outfile3 << runTime.clockTimeIncrement() << "\n";
    outfile3.close();    
}
if(runTime.writeTime())
{
  // D.write();     //控制输出位移场
}