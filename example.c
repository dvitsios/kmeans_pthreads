/*************************************************************************/
/**   File:         example.c                                           **/
/**   Description:  Takes as input a file:                              **/
/**                 ascii  file: containing 1 data point per line       **/
/**                 binary file: first int is the number of objects     **/
/**                              2nd int is the no. of features of each **/
/**                              object                                 **/
/**                 This example performs a fuzzy c-means clustering    **/
/**                 on the data. Fuzzy clustering is performed using    **/
/**                 min to max clusters and the clustering that gets    **/
/**                 the best score according to a compactness and       **/
/**                 separation criterion are returned.                  **/
/**   Author:  Wei-keng Liao                                            **/
/**            ECE Department Northwestern University                   **/
/**            email: wkliao@ece.northwestern.edu                       **/
/**                                                                     **/
/**   Edited by: Jay Pisharath                                          **/
/**              Northwestern University.                               **/
/**                                                                     **/
/**   ================================================================  **/
/**																		**/
/**   Edited by: Sang-Ha  Lee											**/
/**				 University of Virginia									**/
/**																		**/
/**   Description:	No longer supports fuzzy c-means clustering;	 	**/
/**					only regular k-means clustering.					**/
/**					Simplified for main functionality: regular k-means	**/
/**					clustering.											**/
/**   ================================================================  **/
/**																		**/
/**   Edited by: Dimitrios Vitsios										**/
/**				 Aristotle University of Thessaloniki					**/
/**																		**/
/**   Description:	parallel implementation using pthreads              **/
/**																		**/
/*************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/time.h>
#include "getopt.h"
#include "kmeans.h"

extern double wtime(void);
extern float sum;


void usage(char *argv0) {
    char *help =
        "Usage: %s [switches] -i filename\n"
        "       -i filename     :  file containing data to be clustered\n"
        "       -b                 :input file is in binary format\n"
		"       -k                 : number of clusters (default is 8) \n"
        "       -t NUM_THREADS     : number of threads\n";
    fprintf(stderr, help, argv0);
    exit(-1);
}

int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     nclusters=5;
           char   *filename = 0;           
           float  *buf;
           float **attributes;
           float **cluster_centres=NULL;
           int     i, j;           
		   int     NUM_THREADS;
           int     numAttributes;
           int     numObjects;           
           char    line[1024];
           int     isBinaryFile = 0;
           int     nloops;
		   double  timing;



	while ( (opt=getopt(argc,argv,"i:k:t:b"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 'b': isBinaryFile = 1;
                      break;
            case 't': NUM_THREADS=atoi(optarg);
                      break;
            case 'k': nclusters = atoi(optarg);
                      break;
            case '?': usage(argv[0]);
                      break;
            default: usage(argv[0]);
                      break;
        }
    }

    if (filename == 0) usage(argv[0]);

    numAttributes = numObjects = 0;

    /* from the input file, get the numAttributes and numObjects ------------*/
   
    if (isBinaryFile) {
        int infile;
        if ((infile = open(filename, O_RDONLY, "0600")) == -1) {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            exit(1);
        }
        read(infile, &numObjects,    sizeof(int));
        read(infile, &numAttributes, sizeof(int));
   

        /* allocate space for attributes[] and read attributes of all objects */
        buf           = (float*) malloc(numObjects*numAttributes*sizeof(float));
        attributes    = (float**)malloc(numObjects*             sizeof(float*));
        attributes[0] = (float*) malloc(numObjects*numAttributes*sizeof(float));
        for (i=1; i<numObjects; i++)
            attributes[i] = attributes[i-1] + numAttributes;

        read(infile, buf, numObjects*numAttributes*sizeof(float));

        close(infile);
    }
    else {
        FILE *infile;
        if ((infile = fopen(filename, "r")) == NULL) {
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            exit(1);
        }
        while (fgets(line, 1024, infile) != NULL)
            if (strtok(line, " \t\n") != 0)
                numObjects++;
        rewind(infile);
        while (fgets(line, 1024, infile) != NULL) {
            if (strtok(line, " \t\n") != 0) {
                /* ignore the id (first attribute): numAttributes = 1; */
                while (strtok(NULL, " ,\t\n") != NULL) numAttributes++;
                break;
            }
        }
     

        /* allocate space for attributes[] and read attributes of all objects */
        buf           = (float*) malloc(numObjects*numAttributes*sizeof(float));
        attributes    = (float**)malloc(numObjects*             sizeof(float*));
        attributes[0] = (float*) malloc(numObjects*numAttributes*sizeof(float));
        for (i=1; i<numObjects; i++)
            attributes[i] = attributes[i-1] + numAttributes;
        rewind(infile);
        i = 0;
        while (fgets(line, 1024, infile) != NULL) {
            if (strtok(line, " \t\n") == NULL) continue; 
            for (j=0; j<numAttributes; j++) {
                buf[i] = atof(strtok(NULL, " ,\t\n"));
                i++;
            }
        }
        fclose(infile);
    }
  
    nloops = 1;	
	printf("I/O completed\n");

    
	memcpy(attributes[0], buf, numObjects*numAttributes*sizeof(float));

    struct  timeval  first,    second,    lapsed;
    struct  timezone  tzp;
    gettimeofday(&first,  &tzp);

    for (i=0; i<nloops; i++) {
        		
        cluster_centres = NULL;
        cluster(numObjects,
                numAttributes,
                attributes,           /* [numObjects][numAttributes] */
                nclusters,
                NUM_THREADS,
                &cluster_centres   
               );

     
    }

    
    gettimeofday  (&second,  &tzp);
    if  (first.tv_usec  >  second.tv_usec)  {
        second.tv_usec  +=  1000000;
        second.tv_sec--;
    }
    lapsed.tv_usec  =  second.tv_usec  -  first.tv_usec;
    lapsed.tv_sec  =  second.tv_sec  -  first.tv_sec;
 
 

    printf("number of Clusters %d\n",nclusters); 
    printf("number of Attributes %d\n\n",numAttributes); 
    printf("Cluster Centers Output\n"); 
    printf("The first number is cluster number and the following data is arribute value\n");
    printf("=============================================================================\n\n");
    
    for (i=0; i<nclusters; i++) {
      printf("%d: ", i);
      for (j=0; j<numAttributes; j++)
	printf("%f ", cluster_centres[i][j]);
      printf("\n\n");
    }
    
 
    printf("\nTotal sum of distances: %f\n\n",sum);
    
    printf("Time  elapsed  during  thread  execution:  %d,%d\n\n\n",lapsed.tv_sec,lapsed.tv_usec);
    
    free(attributes);
    free(cluster_centres[0]);
    free(cluster_centres);
    free(buf);
    pthread_exit(NULL);
}

