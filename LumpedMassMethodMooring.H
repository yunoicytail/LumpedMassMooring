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

#ifndef LumpedMassMethodMooring_H
#define LumpedMassMethodMooring_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace sixDoFRigidBodyMotionRestraints
{

/*-----------------------------------------------------------------------------------------------------------------------*\
			Class LumpedMassMethodMooring Declaration
\*-----------------------------------------------------------------------------------------------------------------------*/

class LumpedMassMethodMooring
:
        public sixDoFRigidBodyMotionRestraint
{
        // Private data

	//Anchor point, where the spring is attached to an immovable
	// object
	point anchor_;

	//- Reference point of attachment to the solid body
	point refAttachmentPt_;
	
	//- Total length of the mooring chain or cable
	scalar cableLength_;

	//- Total mass of the mooring chain or cable
	scalar totalMass_;

	//- The tensile stiffness of the mooring chain or cable
	scalar EA_;

	//- The number of segments
	scalar n_;

	//- The acceleration of gravity
	scalar g_;

	//- The cable length without lying on the sea bed, AKA the pendulous length
	mutable scalar catenaryCableLength;

	//- The cable mass without lying on the sea bed, AKA the pendulous mass
	mutable scalar catenaryTotalMass;

	//- The vertical mooring force of the sea bed acting on the mooring cable
	mutable scalar supportForcePercent;

	//- The intersection angle between the sea bed and mooring cable
	mutable scalar firstPieceSeita;

	//- The calculated value of the vertical position of the top side of the mooring cable
	mutable scalar calculatedZ;

	//- The calculated value of the horizontal position of the top side of the mooring cable
	mutable scalar calculatedX;

	//- The horizontal component of the mooring force
	mutable scalar Th;

	//- The vertical component of the mooring force
	mutable scalar Tv;

	//- The mooring force
	mutable scalar Ttotal;

	//- The state of the mooring cable
	mutable scalar mooringstate;
public:
        //- Runtime type information
        TypeName("LumpedMassMethodMooring");

        // Constructors

	//- Construct from components
	LumpedMassMethodMooring
	(
		const word& name,
		const dictionary& sDoFRBMRDict
	);

	//- Construct and return a clone
	virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const
	{
		return autoPtr<sixDoFRigidBodyMotionRestraint>
		(
		new LumpedMassMethodMooring(*this)
		);
	}

        //- Destructor
        virtual ~LumpedMassMethodMooring();

        // Member Functions
	//- Calculate the restraint position, force and moment.
	// Global reference frame vectors.
	virtual void restrain
	(
		const sixDoFRigidBodyMotion& motion,
		vector& restraintPosition,
		vector& restraintForce,
		vector& restraintMoment
	) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& sDoFRBMRCoeff);

        //- Write
        virtual void write(Ostream&) const;

        void calculatingZX()const;

        void addSeita()const;

        void calculatingFromNextPiece()const;

        void changeSupportForcePercent()const;

};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //