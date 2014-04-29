void Isynap1(double *vm);
void IBackGrnd(double *vm);
void AuxRffTotal();
void RffTotal(double theta, double t);
void Gff(double theta, double t);
double XCordinate(int neuronIdx, double nA);
double YCordinate(int neuronIdx, double nA);
void GenConMat02();
double gaussrand( long int *seed );
void LinSpace(double startVal, double stopVal, double stepSize, double *outVector, int* nSteps);
void GenSparseConMat(sparseMat *sPtr[]);
void FreeSparseMat(sparseMat *sPtr[]);

