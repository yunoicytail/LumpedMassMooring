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
    
    Caculate the 3D;
    Modified by Yuno.

\*-----------------------------------------------------------------------------------------------------------*/
#include "LumpedMassMethodMooring3D.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "DynamicList.H"
// * * * * * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug (LumpedMassMethodMooring3D, 0);
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        LumpedMassMethodMooring3D,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::LumpedMassMethodMooring3D
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
    gMag_ = Foam::mag(g_);
    g_ /= Foam::mag(g_);
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

Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::~LumpedMassMethodMooring3D()
{}

// * * * * * * * * * * * * * * * * * * * * * Member Functions * * * * * * ** * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::restrain
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
    point pos0 = restraintPosition;
    point pos1 = anchor_;
    // Set the local coordinate system
    scalar y0 = (pos0 & g_);
    scalar y1 = (pos1 & g_ );

    scalar x0 = 0.0;

    //系泊点与锚固点直线长度
    vector r = restraintPosition - anchor_; 
    scalar strightLineLength = mag(r);
    scalar x1 = Foam::mag(r - (r & g_)*g_);
    scalar span_ = x1;//悬链线纵向长度
    scalar h = y1 - y0; //悬链线垂向高度
    vector horzVec = r - (r & g_)*g_;
    horzVec /= Foam::mag(horzVec);
    //求解系泊锚固点之间的标准悬链线长度 - l_p 
    //迭代求解，悬链线方程 y=a(coshx/a -1)
    scalar a1=999999; 
    scalar a2=0.0001; 
    scalar toleranceA=0.001;
    while (Foam::mag(((a1+a2)/2)*(Foam::cosh(span_/((a1+a2)/2))-1)-h) > toleranceA)
    {
        scalar y3 = ((a1+a2)/2) *(Foam::cosh ( span_/((a1+a2)/2)) -1);
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
    scalar maxCatenarycableLength = a * Foam::sinh (span_/a);

    // * * * * * * * * * * * 三种情况计算缆绳力 * * * * * * * * * * * * * * * * * //

    if (cableLength_ >= maxCatenarycableLength) //缆绳长度大于 l_p时，部分拖地
    {
        mooringstate=1;
        mooringstate1loopingZ:

            calculatingZX();
            while (calculatedZ < h) //增加seita继续迭代h
            {
                addSeita();
                calculatingZX();
            }
            while (calculatedX > span_)
            {
                calculatingFromNextPiece();
                goto mooringstate1loopingZ;
            }

        //restraintForce = Th * horzVec; //水平力
        //restraintForce = Tv * g_; //垂向力
        restraintForce = Tv * g_ + Th * horzVec;
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
            while (calculatedX < span_)
            {
                changeSupportForcePercent();
                goto mooringstate2loopingZ;
            }

        //restraintForce = Th * horzVec; //水平力
        //restraintForce = Tv * g_; //垂向力
        restraintForce = Tv * g_ + Th * horzVec;
    }

    else if(strightLineLength > cableLength_)
    {
        mooringstate=3; //The third state of the mooring chain, which is mentioned in section 4.1
        r/= (strightLineLength + VSMALL);
        scalar elongation = strightLineLength - cableLength_;
        scalar stiffness = EA_ / cableLength_;
        restraintForce = -stiffness * elongation * r;
    }
    
    // Data output
    if (motion.report())
    {
        Info<< " Force : "<<restraintForce
            << " horzVec " <<horzVec
            << " mooringstate: " << mooringstate
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
    

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::calculatingZX()const
{
    scalar l = catenaryCableLength / n_;  //每小段长度
    scalar verticalForce = (catenaryTotalMass/n_/2 * gMag_)*(1-supportForcePercent) ;//节点垂直力
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
        verticalForce = verticalForce + catenaryTotalMass * gMag_ / n_ ; //垂直力一段一段的加起来
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

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::addSeita()const
{
    firstPieceSeita += 0.0001;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::calculatingFromNextPiece()const
{
    catenaryCableLength = catenaryCableLength - catenaryCableLength/n_ ;
    catenaryTotalMass = catenaryTotalMass - catenaryTotalMass/n_ ;
}

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::changeSupportForcePercent()const
{
    supportForcePercent = supportForcePercent - 1;
}

bool Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::read
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

void Foam::sixDoFRigidBodyMotionRestraints::LumpedMassMethodMooring3D::write
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
