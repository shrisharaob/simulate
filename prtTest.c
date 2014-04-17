#include<stdio.h>

/* void test(int a, int *l){
  (*l)+= a;
}
*/
main()
{
  int k;
  int x[] = {1, 2, 3, 4};
  int y[2][4] = {{2,4,6,8}, {3,6,9,12}};
  FILE *fp;
  fp = fopen("tstFile", "w");
  printf("line this\n");
  for(k = 0; k < 4; ++k) 
    {
      printf("%d\n", k+1);
      fprintf(fp, "%d %d %d\n", x[k], y[0][k], y[1][k]);
    }
  //   printf("line this\n");
  fclose(fp);
}
