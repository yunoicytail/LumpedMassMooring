/*------------------------------------------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
----------------------------------------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    于含 系泊载液浮体水动力特性的数值及试验研究
    提供的基于集中质量法的锚链力算法

\*-----------------------------------------------------------------------------------------------------------*/
#include "LumpedMassMethodMooring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "DynamicList.H"
// * * * * * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug (LumpedMassMethodMooring, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        LumpedMassMethodMooring,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::LumpedMassMethodMooring
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict),
    anchor_(),
    refAttachmentPt_(),
    cableLength_(),
    totalMass_(),
    EA_(),
    n_(),
    g_()
{
    read(sDoFRBMRDict);
    catenaryCableLength = cableLength_;
    catenaryTotalMass = totalMass_ ;
    supportForcePercent= 0;
    firstPieceSeita = 0.0001;
    Ttotal = 0;
    Th =0;
    Tv =0;
    calculatedZ=0;
    calculatedX=0;
    mooringstate=-1;
}

// * * * * * * * * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::~LumpedMassMethodMooring()
{}

// * * * * * * * * * * * * * * * * * * * * * Member Functions * * * * * * ** * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);//获取当前时间步下浮体锚固点坐标
    restraintMoment = vector::zero;  //初始化力与力矩
    restraintForce = vector::zero;
    scalar direct1=(restraintPosition[0]==anchor_[0]?0:(restraintPosition[0]-anchor_[0]>0?-1:1));
    scalar direct2=(restraintPosition[1]==anchor_[1]?0:(restraintPosition[1]-anchor_[1]>0?-1:1));
    //上面这两个direct指的是力的方向即确定正负 *没有z是因为悬链z方向力肯定是负的*

    scalar x = mag(anchor_[0] - restraintPosition[0]); //进入悬链线平面进行二维计算
    scalar h = mag(anchor_[1] - restraintPosition[1]);
    scalar a1=999999; //以下是求当前高度下**恰好anchor无垂直力**的悬链线长度 l_p
    scalar a2=0.0001; //迭代求解，悬链线方程 y=a(coshx/a -1)
    scalar toleranceA=0.001;
    while (Foam::mag(((a1+a2)/2)*(Foam::cosh(x/((a1+a2)/2))-1)-h) > toleranceA)
    {
        scalar y3 = ((a1+a2)/2) *(Foam::cosh ( x/((a1+a2)/2)) -1);
        if(y3<h)
        {
            a1= (a1+a2)/2;
        }
        if(y3>h)
        {
            a2 = (a1+a2)/2;
        }
    }
    scalar a=(a1+a2)/2;
    scalar maxCatenarycableLength = a * Foam::sinh (x/a); //得到 l_p
    vector r = restraintPosition - anchor_; //两个锚固点连线，定义为矢量
    scalar strightLineLength = mag(r);

    if (cableLength_ >= maxCatenarycableLength) //缆绳长度大于 l_p时，部分拖地
    {
        mooringstate=1;//The first state of the mooring chain, which is mentioned in section 4.1
        mooringstate1loopingZ:

            calculatingZX();
            while (calculatedZ < h) //增加seita继续迭代h
            {
                addSeita();
                calculatingZX();
            }
            while (calculatedX > x)
            {
                calculatingFromNextPiece();
                goto mooringstate1loopingZ;
            }

        restraintForce[0]=direct1*Th;
        restraintForce[1]=direct2*Tv;
        restraintForce[2]=0;

    }

    else if (maxCatenarycableLength > cableLength_ && cableLength_ >= strightLineLength)
    {
        mooringstate=2; //The second state of the mooring chain, which is mentioned in section 4.1
        mooringstate2loopingZ:

            calculatingZX();
            while (calculatedZ < h)
            {
                addSeita();
                calculatingZX();
            }
            while (calculatedX < x)
            {
                changeSupportForcePercent();
                goto mooringstate2loopingZ;
            }

        restraintForce[0]=direct1*Th;
        restraintForce[1]=direct2*Tv;
        restraintForce[2]=0;
    }

    else if(strightLineLength > cableLength_)
    {
        mooringstate=3; //The third state of the mooring chain, which is mentioned in section 4.1
        r/= (strightLineLength + VSMALL);
        scalar elongation = strightLineLength - cableLength_;
        scalar stiffness = EA_ / cableLength_;
        restraintForce = -stiffness * elongation * r;
        restraintForce[2]=0;
    }

    if (motion.report())
    {
        Info<< " magForce : "<<mag(restraintForce)
            //<< " mooringstate: " << mooringstate
            //<< " restraintPosition: " <<restraintPosition
            //<< " anchor: " <<anchor_
            //<< " seita " <<firstPieceSeita
            //<< " calculatedZ " <<calculatedZ
            //<< " calculatedX " <<calculatedX
            << endl;
    }
    catenaryCableLength = cableLength_;
    catenaryTotalMass = totalMass_ ;
    supportForcePercent= 0;
    firstPieceSeita = 0.0001;
    Ttotal = 0;
    Th =0;
    Tv =0;
    calculatedZ=0;
    calculatedX=0;
    mooringstate=-1;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::calculatingZX()const
{
    scalar l = catenaryCableLength / n_;  //每小段长度
    scalar verticalForce = (catenaryTotalMass/n_/2 * g_)*(1-supportForcePercent) ;//节点垂直力
    //the horizontalforce in the first piece which also in every piece
    scalar horizontalForce = verticalForce/Foam::tan(firstPieceSeita);
    scalar seita = firstPieceSeita;
    scalar totalForce = verticalForce/Foam::sin(seita); //该段总力
    scalar finalZ=0;
    scalar finalX=cableLength_ - catenaryCableLength;

    finalZ = finalZ + l*(1+totalForce/EA_)*Foam::sin(seita);
    finalX = finalX + l*(1+totalForce/EA_)*Foam::cos(seita);

    for (int i=1;i<=n_ - 1; i++)
    {
        //calculate the vertical force in each piece
        verticalForce = verticalForce + catenaryTotalMass * g_ / n_ ; //垂直力一段一段的加起来
        //calculate the angle between each cable piece and horizontal
        seita = Foam::atan(verticalForce/horizontalForce) ;
        //calculate total force in each piece
        totalForce = verticalForce/Foam::sin(seita);
        finalZ = finalZ + l*(1+totalForce/EA_)*Foam::sin(seita);
        finalX = finalX + l*(1+totalForce/EA_)*Foam::cos(seita);
    }
    calculatedZ=finalZ;
    calculatedX=finalX;
    Ttotal=totalForce;
    Tv=verticalForce; //垂直力
    Th=horizontalForce; //水平力
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::addSeita()const
{
    firstPieceSeita += 0.0001;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::calculatingFromNextPiece()const
{
    catenaryCableLength = catenaryCableLength - catenaryCableLength/n_ ;
    catenaryTotalMass = catenaryTotalMass - catenaryTotalMass/n_ ;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::changeSupportForcePercent()const
{
    supportForcePercent = supportForcePercent - 1;
}

bool Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);
    sDoFRBMRCoeffs_.lookup("anchor") >> anchor_;
    sDoFRBMRCoeffs_.lookup("refAttachmentPt") >> refAttachmentPt_;
    sDoFRBMRCoeffs_.lookup("cableLength") >> cableLength_;
    sDoFRBMRCoeffs_.lookup("totalMass") >> totalMass_;
    sDoFRBMRCoeffs_.lookup("EA") >> EA_;
    sDoFRBMRCoeffs_.lookup("n") >> n_;
    sDoFRBMRCoeffs_.lookup("g") >> g_;
    return true;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring::write
(
    Ostream& os
) const
{
    os.writeKeyword("anchor")
        << anchor_ << token::END_STATEMENT << nl;
    os.writeKeyword("refAttachmentPt")
        << refAttachmentPt_ << token::END_STATEMENT << nl;
    os.writeKeyword("cableLength")
        << cableLength_ << token::END_STATEMENT << nl;
    os.writeKeyword("totalMass")
        << totalMass_ << token::END_STATEMENT << nl;
    os.writeKeyword("EA")
        << EA_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("g")
        << g_ << token::END_STATEMENT << nl;
}
// ******************************************************************************* //
