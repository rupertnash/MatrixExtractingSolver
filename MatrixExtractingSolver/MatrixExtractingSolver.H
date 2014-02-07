#ifndef MATRIXEXTRACTINGSOLVER_H
#define MATRIXEXTRACTINGSOLVER_H

// Create the header CrossPlatform.H in this directory.
// If you are using OpenFoam-extend, then you must define
// the preprocessor macro ON_EXTEND in that file.
// If you are not using -extend, then an empty file suffices.
#include "CrossPlatform.H"

#include "lduMatrix.H"

// Typedefs to be compatible between vanilla and extend
#ifdef ON_EXTEND
typedef Foam::lduMatrix::solverPerformance SolverPerformance;
#else
typedef Foam::solverPerformance SolverPerformance;
#endif

namespace Foam
{
  // Forward declare Foam::Time
  class Time;
  
  class MatrixExtractingSolver : public lduMatrix::solver
  {
    // Private Member Functions

    //- Disallow default bitwise copy construct
    MatrixExtractingSolver(const MatrixExtractingSolver&);

    //- Disallow default bitwise assignment
    void operator=(const MatrixExtractingSolver&);
        
  public:
    //- Runtime type information
    TypeName("MatrixExtractingSolver");
    
    //- Construct from matrix components and solver data stream
        MatrixExtractingSolver
	(
	 const word& fieldName,
	 const lduMatrix& matrix,
	 const FieldField<Field, scalar>& coupleBouCoeffs,
	 const FieldField<Field, scalar>& coupleIntCoeffs,
	 const lduInterfaceFieldPtrsList& interfaces,
	 const dictionary& dict
	 );
    // D'tor
    virtual ~MatrixExtractingSolver()
    {
      // Use of autoPtr means we don't need to do anything here.
    }

    // Member functions

    //- Solve the matrix with this solver
    virtual ::SolverPerformance solve
    (
     scalarField& x,
     const scalarField& b,
     const direction cmpt = 0
     ) const;

  protected:
    // Protected Member Functions                                                                                                                                                                 

    //- Read the control parameters from the controlDict_

    // NOTE: The OF implementors are clearly unaware that one cannot
    // call virtual functions from base class constructors. To work
    // around, call this from your derived class c'tor.
    virtual void readControls();
    
    virtual bool shouldWrite() const;
  private:
    // This is the actual solver to which we delegate the work.
    autoPtr<lduMatrix::solver> worker;
    // Reference to the application's Time object
    const Time& appTime;
  };
} // End namespace Foam

#endif // MATRIXEXTRACTINGSOLVER_H
