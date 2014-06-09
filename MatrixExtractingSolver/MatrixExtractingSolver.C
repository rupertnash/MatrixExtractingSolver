#include "MatrixExtractingSolver.H"
#include "fileStat.H"
#include "Time.H"
#include "OFstream.H"
#include "processorFvPatchField.H"

#include <fstream>
#include <unistd.h>

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
   ),
  // Get the runTime object
  appTime(matrix.mesh().thisDb().time())
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

  // Work around different APIs on solver to get the control dictionary...
#ifdef ON_EXTEND
  const dictionary& solverDict = dict();
#else
  const dictionary& solverDict = controlDict_;
#endif

  const dictionary& workerDict = solverDict.subDict("worker");
  worker->read(workerDict);
}

// Helper function for getting the coefficients (because there are multiple possible types of patch field).
template <typename T>
void AddInterface(Foam::dictionary& interfacesDict, const Foam::processorFvPatchField<T>* ptr, const Foam::scalarField& coupleCoeffs)
{
  const int neigh = ptr->neighbProcNo();
  const Foam::word neighName(Foam::word("processor") + Foam::name(neigh));
  Foam::dictionary thisInterfaceDict;
  
  const Foam::fvPatch& patch = ptr->patch();
  const Foam::unallocLabelList& localCellIds = patch.faceCells();
  
  thisInterfaceDict.add("localCellIds", localCellIds);
  thisInterfaceDict.add("coeffs", coupleCoeffs);
  
  interfacesDict.add(neighName, thisInterfaceDict);
}

Foam::dictionary Foam::MatrixExtractingSolver::GetInterfaceCoeffs(const direction cmpt) const
{
  dictionary interfaceDict;
  forAll (interfaces_, interfaceI)
  {
    if (interfaces_.set(interfaceI))
    {
      const scalarField& coeffs = coupleBouCoeffs_[interfaceI];
      if (const processorFvPatchField<Vector<double> >* ptr = dynamic_cast<const processorFvPatchField<Vector<double> >*>(interfaces_(interfaceI)))
      {
	AddInterface(interfaceDict, ptr, coeffs);
      }
      else if (const processorFvPatchField<double>* ptr = dynamic_cast<const processorFvPatchField<double>*>(interfaces_(interfaceI)))
      {
	AddInterface(interfaceDict, ptr, coeffs);
      }
      else
      {
	// Error
	FatalErrorIn("Foam::MatrixExtractingSolver::solve(scalarField&, const scalarField&, const direction) const")
	  << "Can't cast interface to a concrete type" << exit(FatalError);
      }
    }
  }
  return interfaceDict;
}

//- Solve the matrix with this solver
SolverPerformance Foam::MatrixExtractingSolver::solve
(
 scalarField& x,
 const scalarField& b,
 const direction cmpt
 ) const
{
  // Delegate actually solving to the worker.

  // Do this first so it sets up the solverPerformance object which we
  // can query for useful things.
  ::SolverPerformance sPerf = worker->solve(x, b, cmpt);
  
  if (shouldWrite())
  {
    word dictMetaName = fieldName() + ".metadata";
    IOdictionary* metaDict = NULL;
    
    objectRegistry::const_iterator item = appTime.find(dictMetaName);
    if (item == appTime.objectRegistry::end())
      {
	metaDict = new IOdictionary
	  (
	   IOobject
	   (
	    dictMetaName,
	    appTime.constant(),
	    appTime,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	    )
	   );
	metaDict->store();
      }
    else
      {
	metaDict = dynamic_cast<IOdictionary*>(item());
      }
    
    label lastWriteTime = metaDict->lookupOrDefault("time", -1, false, false);
    label subCycle = 0;
    if (lastWriteTime == appTime.timeIndex())
      {
	// Iteration > 0
	subCycle = metaDict->lookupOrDefault("iteration", subCycle, false, false) + 1;
      }
    
    word dictName;
    {
      std::ostringstream ss;
      ss << fieldName() << "." << subCycle << ".system"; //makeDictName(subCycle);
      dictName = ss.str();
    }
    item = appTime.find(dictName);
    IOdictionary* systemDict = NULL;
    if (item == appTime.objectRegistry::end())
      {
	// Construct it
	systemDict = new IOdictionary
	  (
	   IOobject
	   (
	    dictName,
	    appTime.timeName(),
	    appTime,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	    )
	   );
	systemDict->store();
      }
    else
      {
	systemDict = dynamic_cast<IOdictionary*>(item());
      }
    systemDict->clear();
    
    // Add the source vector
    systemDict->add("source", b);

    // Construct a dictionary to hold the LDU matrix
    dictionary matrixDict;
    
    // Diagonal first. If not present, then the matrix is in an error
    // state.
    if (!matrix_.hasDiag())
      FatalErrorIn("Foam::MatrixExtractingSolver::solve(scalarField&, const scalarField&, const direction) const")
	<< "Matrix lacks diagonal" << exit(FatalError);
    
    const scalarField& diag = matrix_.diag();
    matrixDict.add("diag", diag);

    if (matrix_.hasLower()) 
    {
      const scalarField& lo = matrix_.lower();
      matrixDict.add("lower", lo);
    }

    if (matrix_.hasUpper()) 
    {
      const scalarField& up = matrix_.upper();
      matrixDict.add("upper", up);
    }
    
    // Only non-diagonal matrices need the indexing helpers.
    if (!matrix_.diagonal())
    {
      // Get the addressing helper
      const lduAddressing& addr = matrix_.lduAddr();
  
      const labelList& upAddr = addr.upperAddr();
      const labelList& loAddr = addr.lowerAddr();
      matrixDict.add("upperAddr", upAddr);
      matrixDict.add("lowerAddr", loAddr);
    }
    
    // Add the dict representation of the matrix to the system
    systemDict->add("matrix", matrixDict);
    
    dictionary interfaceDict = GetInterfaceCoeffs(cmpt);
    systemDict->add("interfaces", interfaceDict);
    
    // Write
    systemDict->regIOobject::write();
    
    // Update the cache metadata.
    metaDict->set("time", appTime.timeIndex());
    metaDict->set("iteration", subCycle);
  }
  
  return sPerf;
}

bool Foam::MatrixExtractingSolver::shouldWrite() const
{
  return true;
}

