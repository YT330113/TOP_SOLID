/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    TOP_SOLID

Description

    这是我决定从0开始写的第一个拓扑优化程序，用于固体结构拓扑优化！
    Topology optimization of linear elasticity problem.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"                      // 下面的 creatControl.H 所需

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "********Topology optimization of linear elasticity problem.********"  <<nl;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "readMechanicalProperties.H" 
    #include "opt_initialization.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

     while (simple.loop(runTime))
    {
        #include "update.H"
        #include "LinearElasticity.H"
        #include "costfunction.H"              
        #include "sensitivity.H"
        #include "topSolver.H"
    }

    #include "finalize.H"

    Foam::Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Foam:: Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
