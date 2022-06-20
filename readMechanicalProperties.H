/********************************************************************************/
/*****************读取材料属性并创建位移场*************************************/
/********************************************************************************/

/************************读取物理属性*******************************************/

Info<< "Reading mechanical properties\n" << endl;
IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
int iCorr = 0, nCorr = 300; // 这两行是solidDisplacement的fvSolution里stressAnalysis子字典里的内容
                            // nCorr就是“nCorrectors”（稳态为1），convergenceTolerance就是“D”，iCorr和nCorr共同控制位移平衡方程在每一个时间步长的求解迭代次数
scalar convergenceTolerance(1e-5), initialResidual(0);

Switch planeStress(mechanicalProperties.lookup("planeStress"));                                                // 是否平面应力问题
dimensionedScalar rhoE("rhoE", dimMass/dimLength/sqr(dimTime),readScalar(mechanicalProperties.lookup("rhoE")));// 定义弹性模量及单位
dimensionedScalar rho("rho", dimMass / dimVolume, readScalar(mechanicalProperties.lookup("rho")));             // 定义密度及单位
scalar Po(readScalar(mechanicalProperties.lookup("Po")));                                                      //定义无量刚数泊松比
dimensionedScalar E(rhoE / rho);                         // 弹性模量正则化，因为解方程时方程两边是同除以rho的 ; 这里的E就是拓扑列式里的Emax
dimensionedScalar Emin(E * 1e-9);                      // Emin=Emax*10^-9
volScalarField mu(xh *xh *xh *E / (2.0 * (1.0 + Po))); //下面两行定义并初始化拉梅常数，openfoam求解位移方程所需
volScalarField lambda(xh *xh *xh *Po *E / ((1.0 + Po) * (1.0 - 2.0 * Po)));
//----------------------------------------------------------------------------------------------------------//

/************************创建位移场等*******************************************/

Info<< "Reading field D\n" << endl;
volVectorField D                                            //// 定义位移向量场
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD                                    //// 定义位移梯度张量场
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD                  // 定义应力张量场，对陈张量场 
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradD) + (lambda*I)*tr(gradD)//应力应变与位移关系显式求得
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
//这与openfoam位移方程求解的分离式解法有关，显示与隐士离散迭代求解，具体见PDF
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

Switch compactNormalStress("yes");//是否 compact Normal Stress，需要修改divSigmaExp，不懂

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());         //？？？
//----------------------------------------------------------------------------------------------------------//

/***********************************定义灵敏度场********************************/

volScalarField fsens //目标函数灵敏度，对x的导数,用于MMA求解
    (
        IOobject(
            "fsens",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
            ),
        -gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD))),////////////////////////////////////////////////////////////////////这句是我改过的，不知有没有问题
        zeroGradientFvPatchScalarField::typeName
        );
        volScalarField fsensh(fsens); //目标函数灵敏度，对xh的导数

volScalarField gsens          //约束函数灵敏度，对x的导数                              
(
    IOobject
    (
        "gsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsensh(gsens); //约束函数灵敏度，对xh的导数
//----------------------------------------------------------------------------------------------------------//
