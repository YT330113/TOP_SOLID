/****************************************************************************/
/*****************线弹性平衡方程求解****************************************/
/****************************************************************************/

Info<< "\nslove linear elasticity equation\n" << endl;

initialResidual = 0;  //初始误差初始化
iCorr = 0;//iCorr和nCorr共同控制位移平衡方程在每个时间步长的求解迭代次数 ，nCorr越大每一步收敛越好

do
{
  {
      fvVectorMatrix DEqn
      (
          fvm::d2dt2(D)
       ==
          fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
        + divSigmaExp
        //- fvc::grad(threeKalpha*T)
      );

      initialResidual = DEqn.solve().max().initialResidual();//更新initialResidual

      if (!compactNormalStress)
      {
          divSigmaExp = fvc::div(DEqn.flux());
      }
  }

  {
      volTensorField gradD(fvc::grad(D));
      sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

      if (compactNormalStress)
      {
          divSigmaExp = fvc::div
          (
              sigmaD - (2*mu + lambda)*gradD,
              "div(sigmaD)"
          );
      }
      else
      {
          divSigmaExp += fvc::div(sigmaD);
      }
  }

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);