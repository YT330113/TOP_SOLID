/***********************************************************/
/**********************创建设计变量场**********************/
/***********************************************************/

Info << "Reading optimization parameters\n" << endl;

/********************从optProperties字典文件中读取参数*****************/

IOdictionary optProperties                      //创建字典，用于输入输出
(
   IOobject
   (
       "optProperties",                                     // 字典名称
       runTime.constant(),                              // 字典路径
       mesh,
       IOobject::MUST_READ_IF_MODIFIED, //字典读控制（从字典文件读取到程序）
       IOobject::NO_WRITE                                   // 字典写控制（程序写入到字典文件）
   )
);
scalar voluse(readScalar(optProperties.lookup("voluse")));      // 体积约束
Switch solid_area(optProperties.lookup("solid_area"));              //bool型变量，判断是否存在不可设计域
scalar filterR(readScalar(optProperties.lookup("filterR")));
dimensionedScalar b("b", dimensionSet(0, -2, 0, 0, 0, 0, 0), 1);// b就是PDE求导数中的1/( r^2 ),后续会定义具体值  

//-------------------------------------------------------------------------------//

/**************************定义伪密度场gamma**************************/

// xp是PDE过滤的；xh是Heviside过滤的，xh是最终更新的设计变量，
//用于下一次迭代时控制方程的求解 ; x是优化变量，即MMA优化求解器中的当前点
Info << "Reading pseudo-density field\n"<< endl;
volScalarField xp                           //创建体标量场
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    // 如果在算例0文件夹下存在gamma的初边值文件，则读入文件中的值
    //作为设计变量的初值；否则采用voluse作为初值。
    zeroGradientFvPatchScalarField::typeName   
    // 指定场边界信息，所有边界均为zeroGradient 边界
    // 体心标量场类型的数据不仅包含标量场在控制体积中心存储的离散值，
    // 还包含了该场的边界信息
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField x(xh);

volScalarField drho(xh);                            //这个drho是xh对xp的导数
//--------------------------------------------------------------------------------------------//


/*****************************定义非设计域****************************/

labelList cells;                                        //定义标签列表，存放非设计域单元编号

if (solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName); //找到solid_area cellzone的序号，
                                                                            //注意是整个solid_area的序号，不是每个单元的序号
  cells=mesh.cellZones()[zoneID];       //将非设计域单元编号存在cells列表里

  forAll(cells, i)
  {
     x[cells[i]]=1.0;                                           //对非设计域单元将设计变量置1
  }
}
//---------------------------------------------------------------------------------------------//

