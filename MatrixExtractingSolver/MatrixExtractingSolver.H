#ifndef MATRIXEXTRACTINGSOLVER_H
#define MATRIXEXTRACTINGSOLVER_H

#include "lduMatrix.H"

namespace Foam
{
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
    virtual solverPerformance solve
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

    // Helper functions
    void writeVector(const scalarField& x, const std::string& fileName) const;
    void writeMatrix(const lduMatrix& A, const std::string& fileName) const;

  private:
    // This is the actual solver to which we delegate the work.
    autoPtr<lduMatrix::solver> worker;

  };
} // End namespace Foam

#endif // MATRIXEXTRACTINGSOLVER_H