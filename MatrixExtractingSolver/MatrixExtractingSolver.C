#include "MatrixExtractingSolver.H"
#include "fileStat.H"
#include <fstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(MatrixExtractingSolver, 0);

  lduMatrix::solver::addasymMatrixConstructorToTable<MatrixExtractingSolver>
  addMatrixExtractingSolverAsymMatrixConstructorToTable_;

  lduMatrix::solver::addsymMatrixConstructorToTable<MatrixExtractingSolver>
  addMatrixExtractingSolverSymMatrixConstructorToTable_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MatrixExtractingSolver::MatrixExtractingSolver
(
 const word& fieldName,
 const lduMatrix& matrix,
 const FieldField<Field, scalar>& coupleBouCoeffs,
 const FieldField<Field, scalar>& coupleIntCoeffs,
 const lduInterfaceFieldPtrsList& interfaces,
 const dictionary& dict
 )
  :
  lduMatrix::solver
  (
   fieldName,
   matrix,
   coupleBouCoeffs,
   coupleIntCoeffs,
   interfaces,
   dict
   )
{
  // Initialise the delegate
  const dictionary& workerDict = dict.subDict("worker");
  worker = lduMatrix::solver::New(fieldName, matrix, coupleBouCoeffs, coupleIntCoeffs, interfaces, workerDict);
  readControls();
}


// NOTE: The OF implementors are clearly unaware that one cannot
// call virtual functions from base class constructors. To work
// around, call this from your derived class c'tor.
void Foam::MatrixExtractingSolver::readControls()
{
  lduMatrix::solver::readControls();
  const dictionary& workerDict = controlDict_.subDict("worker");
  worker->read(workerDict);
}

//- Solve the matrix with this solver
Foam::solverPerformance Foam::MatrixExtractingSolver::solve
(
 scalarField& x,
 const scalarField& b,
 const direction cmpt
 ) const
{
  // Delegate actually solving to the worker.

  // Do this first so it sets up the solverPerformance object which we
  // can query for useful things.
  solverPerformance sPerf = worker->solve(x, b, cmpt);
  
  // Figure out a base file name
  word baseName = sPerf.fieldName() + sPerf.solverName();
  // Write the system.
  writeVector(b, baseName + ".vec");
  writeMatrix(matrix_, baseName + ".coo");

  return sPerf;
}

void Foam::MatrixExtractingSolver::writeVector(const scalarField& x, const std::string& fileName) const
{
  Info << "Write vector to file " << fileName << endl;
  // If the file exists, bail out.
  fileStat status(fileName);
  if (status.isValid())
    FatalErrorIn("Foam::MatrixExtractingSolver::writeVector(const scalarField&, const std::string&) const")
      << "File '" << fileName << "' already exists." << exit(FatalError);
  
  //Open the file.
  std::ofstream file(fileName.c_str());
  scalarField::const_iterator ptr = x.begin(), end = x.end();
  for (; ptr != end; ++ptr) {
    file << *ptr << std::endl;
  }
}

void Foam::MatrixExtractingSolver::writeMatrix(const lduMatrix& A, const std::string& fileName) const
{
  Info << "Write matrix to file " << fileName << endl;
  // If the file exists, bail out.
  fileStat status(fileName);
  if (status.isValid())
    FatalErrorIn("Foam::MatrixExtractingSolver::writeMatrix(const lduMatrix&, const std::string&) const")
      << "File '" << fileName << "' already exists." << exit(FatalError);
  
  //Open the file.
  std::ofstream file(fileName.c_str());
  file << "# row col val" << std::endl;

  // Diagonal first. If not present, then the matrix is in an error
  // state.
  if (!A.hasDiag())
    FatalErrorIn("Foam::MatrixExtractingSolver::writeMatrix(const lduMatrix&, const std::string&) const")
      << "Matrix lacks diagonal" << exit(FatalError);
  
  const scalarField& diag = A.diag();
  scalarField::const_iterator ptr = diag.begin(), end = diag.end();
  unsigned i = 0;
  for (; ptr != end; ++ptr, ++i) {
    file << i << " " << i << " " << *ptr << std::endl;
  }
  
  // Nothing more to be done for a diagonal matrix.
  if (A.diagonal())
    return;

  // Get the addressing helper
  const lduAddressing& addr = A.lduAddr();
  
  // Note that lduMatrix will return the same data for both upper and
  // lower if it's symmetric. This is helpful.
  
  // Get the indexing arrays
  labelList::const_iterator upAddrPtr = addr.upperAddr().begin();
  labelList::const_iterator loAddrPtr = addr.lowerAddr().begin();

  const scalarField& lo = A.lower();
  scalarField::const_iterator loPtr = lo.begin(), loEnd = lo.end();
  labelList::const_iterator& loRowPtr = upAddrPtr;
  labelList::const_iterator& loColPtr = loAddrPtr;

  const scalarField& up = A.upper();
  scalarField::const_iterator upPtr = up.begin(), upEnd = up.end();
  labelList::const_iterator& upRowPtr = loAddrPtr;
  labelList::const_iterator& upColPtr = upAddrPtr;
  
  for(; loPtr != loEnd; ++loPtr, ++upPtr, ++loRowPtr, ++loColPtr)
  {
    file << *loRowPtr << " " << *loColPtr << " " << *loPtr << std::endl;
    file << *upRowPtr << " " << *upColPtr << " " << *upPtr << std::endl;
  }
}
