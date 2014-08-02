#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int conMat[2][2];
void ReadConMatFromFile(const char* filename, int* cmat, int nNeurons) {
  FILE* fp;
  int i, j;
  char buffer[nNeurons * nNeurons];
  size_t result;
  
  /*  long fileSize;*/
  fp = fopen(filename, "r");

  /*  fseek(fp, 0, SEEK_END);
  fileSize = ftell(fp);
  rewind(fp);
  buffer = (char *)malloc(sizeof(char) * fileSize);
  if(buffer == NULL) {
    fputs("buffer not allocated", stderr);
    exit(1);
    }
  */
  if(fp != NULL) {
    result = fread(buffer, sizeof(int), nNeurons * nNeurons, fp);
    for(i = 0; i < nNeurons; ++i) {
      for(j = 0; j < nNeurons; ++j) {
        cmat[i + j * nNeurons] = buffer[i + j * nNeurons] - '0'; // convert ascii decimal character representation to integer 
        conMat[i][j] = buffer[i + j * nNeurons] - '0'; // convert ascii decimal character representation to integer 
      }
    }
    fclose(fp);
  }
  else {
    printf("%s does not exist\n", filename);
  }
}

int main() {
  FILE* fp;
  int conmat[4] = {0, 1, 1, 0}; 
  int rconmat[4];
  int i, j, n = 2, cm;
  fp = fopen("conVec.csv", "w");
  for(i = 0; i < n*n; ++i){
    fprintf(fp, "%d", conmat[i]);
  }
  fclose(fp);
  /* read from file */
  ReadConMatFromFile("conVec.csv", rconmat, n);
  printf("\n");
  for(i = 0; i < n; ++i) {
    for(j = 0; j < n; ++j) {
      printf("%d ", conMat[i][j]);
    }
    printf("\n");
  }
}
