#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>
#include <fcntl.h>
#include "getopt.h"
#include <pthread.h>
#include "kmeans.h"

#define RANDOM_MAX 2147483647

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38
#endif


extern double wtime(void);


/*************************************************************************/
/**   File:         cluster.c, part: "cluster"                          **/
/**   Description:  Takes as input a file, containing 1 data point per  **/
/**                 per line, and performs a fuzzy c-means clustering   **/
/**                 on the data. Fuzzy clustering is performed using    **/
/**                 min to max clusters and the clustering that gets    **/
/**                 the best score according to a compactness and       **/
/**                 separation criterion are returned.                  **/
/**   Author:  Brendan McCane                                           **/
/**            James Cook University of North Queensland.               **/
/**            Australia. email: mccane@cs.jcu.edu.au                   **/
/**                                                                     **/
/**   Edited by: Jay Pisharath, Wei-keng Liao                           **/
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



int cluster(int      numObjects,      /* number of input objects */
            int      numAttributes,   /* size of attribute of each object */
            float  **attributes,      /* [numObjects][numAttributes] */
            int      num_nclusters,
            int      NUM_THREADS,       /* in: number of threads  */
            float ***cluster_centres /* out: [best_nclusters][numAttributes] */
            )
{
    int     nclusters;
    float **tmp_cluster_centres;
    

    nclusters=num_nclusters;

    srand(7);
	
    tmp_cluster_centres = kmeans_clustering(attributes,
                                            numAttributes,
                                            numObjects,
                                            nclusters,
                                            NUM_THREADS);

                                            

   
    if (*cluster_centres) {
		free((*cluster_centres)[0]);
        free(*cluster_centres);
	}
	*cluster_centres = tmp_cluster_centres;


    return 0;
}


/*************************************************************************/
/**   File:         cluster.c, part: "kmeans_clustering"                **/
/**   Description:  Implementation of regular k-means clustering        **/
/**                 algorithm                                           **/
/**   Author:  Wei-keng Liao                                            **/
/**            ECE Department, Northwestern University                  **/
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


typedef struct{
       float  **feature;  /*[npoints][nfeatures] */  
       int     npoints;
       int     nfeatures;
       float  **clusters;        /* [nclusters][nfeatures] */
       int     nclusters;
       int     rem;
       int     cluster_len;
       int     *membership;
       int     *new_centers_len; /* [nclusters]: no. of points in each cluster */
       float  **new_centers;     /* [nclusters][nfeatures] */
       float    total_sum; 
}POINTS_CENTERS_MATRICES;

POINTS_CENTERS_MATRICES data1;

float sum;
int NUM_THREADS;
pthread_mutex_t mut1;


void *calc_distance(void *arg)
{
    int  i, j, start, end, len;
    float dist;
    long offset;
    offset=(long)arg;
    len=data1.cluster_len;
    start = offset*len;
    if(offset==(NUM_THREADS - 1))
        end = start + len + data1.rem;
    else
        end = start + len;
    
    
    for (i=start; i<end; i++){ 
              dist = euclid_dist_2(data1.feature[i], data1.clusters[data1.membership[i]], data1.nfeatures);  // no need square root 
                  pthread_mutex_lock(&mut1);
                  data1.total_sum+=dist;
                  pthread_mutex_unlock(&mut1);
    }
   
    pthread_exit(NULL);
}

   
void *find_nearest_point(void *arg)
{ 
    int  i, j, l, start, end, len;
    long offset;
    offset=(long)arg;
    len=data1.cluster_len;
    start = offset*len;
    if(offset==(NUM_THREADS - 1))
        end = start + len + data1.rem;
    else
        end = start + len;
 

    /* find the cluster center id with min distance to a point */
    for (i=start; i<end; i++) {
        float max_dist=FLT_MAX;
        for(j=0; j<data1.nclusters; j++){
                 float dist;
                 dist = euclid_dist_2(data1.feature[i], data1.clusters[j], data1.nfeatures);  /* no need square root */
                    if (dist < max_dist) {
                          max_dist = dist;
                          data1.membership[i] = j;
                          data1.new_centers_len[j]++;   //the number of points contained by cluster with index 'j' is increased by '1'
	                      for (l=0; l<data1.nfeatures; l++)          
				          data1.new_centers[j][l] += data1.feature[i][l];
                    }
        }
    }
    pthread_exit(NULL);
}

/* multi-dimensional spatial Euclid distance square */
__inline
float euclid_dist_2(float *pt1,
                    float *pt2,
                    int    numdims)
{
    int i;
    float ans=0.0;

    for (i=0; i<numdims; i++)
        ans += (pt1[i]-pt2[i]) * (pt1[i]-pt2[i]);

    return(ans);
}



float** kmeans_clustering(float **feature,    /* [npoints][nfeatures] */
                          int     nfeatures,
                          int     npoints,
                          int     nclusters,
                          int     num_threads) 
{
    int      i, j, r, n=0, index;
    long     arg, t;
    void    *status, *status2;

    NUM_THREADS = num_threads;
    pthread_t *threads1;
    pthread_t *threads2;
    threads1 = (pthread_t *)malloc(NUM_THREADS * sizeof(pthread_t));
    threads2 = (pthread_t *)malloc(NUM_THREADS * sizeof(pthread_t));
    
    /* allocate space for returning variable data1.clusters[] */
    data1.clusters = NULL;
    data1.clusters = (float **) malloc(nclusters * sizeof(float *));
    data1.clusters[0] = (float *)  malloc(nclusters * nfeatures * sizeof(float));
    for (i=1; i<nclusters; i++)
        data1.clusters[i] = data1.clusters[i-1] + nfeatures;

    /* randomly pick cluster centers */
    for (i=0; i<nclusters; i++) {
        //n = (int)rand() % npoints;
        for (j=0; j<nfeatures; j++)
            data1.clusters[i][j] = feature[n][j];
		n++;
    }
    
    /* need to initialize new_centers_len and new_centers[0] all to 0 */
    data1.new_centers_len = (int *)calloc(nclusters, sizeof(int));

    data1.new_centers = NULL;
    data1.new_centers = (float **)malloc(nclusters * sizeof(float *));
    data1.new_centers[0] = (float *)  calloc(nclusters * nfeatures, sizeof(float));
    for (i=1; i<nclusters; i++)
        data1.new_centers[i] = data1.new_centers[i-1] + nfeatures;
    
    data1.membership = (int *) malloc(npoints * sizeof(int));
    for (i=0; i<npoints; i++)
		data1.membership[i] = -1;
    
    data1.feature=feature;
    data1.nfeatures=nfeatures;
    data1.nclusters=nclusters;
    data1.npoints=npoints;
    data1.total_sum=0;
   	data1.rem=data1.npoints % NUM_THREADS;
    data1.cluster_len=(data1.npoints - data1.rem) / NUM_THREADS;
        
  
    pthread_attr_t attr;
    pthread_attr_t attr2;
    pthread_mutex_init(&mut1, NULL);
    
    for(r=0;r<30;r++){
    
         
         pthread_attr_init(&attr);
         pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
             
         for (arg=0; arg<NUM_THREADS; arg++) {
            
            int rc;
            rc=pthread_create(&threads1[arg], &attr, find_nearest_point, (void *)arg);
            if (rc){
               printf("ERROR; return code from pthread_create() '-1-' is %d\n", rc);
               exit(-1);
            }
         }
        
        pthread_attr_destroy(&attr);
        
     
        /* Wait on the other threads */
        for(i=0; i<NUM_THREADS; i++)
        {
         pthread_join(threads1[i], &status);
        }
	    
	    
	/* replace old cluster centers with new_centers */
    	for (i=0; i<nclusters; i++) {
            for (j=0; j<nfeatures; j++) {
                if (data1.new_centers_len[i] > 0)
					data1.clusters[i][j] = data1.new_centers[i][j] / data1.new_centers_len[i];
				data1.new_centers[i][j] = 0.0;   /* set back to 0 */
			}
			data1.new_centers_len[i] = 0;   /* set back to 0 */
		}
		
    }
    
    pthread_attr_init(&attr2);
    pthread_attr_setdetachstate(&attr2, PTHREAD_CREATE_JOINABLE);
    
      for (t=0; t<NUM_THREADS; t++) {
            
            int rc;
            rc=pthread_create(&threads2[t], &attr2, calc_distance, (void *)t);
            if (rc){
               printf("ERROR; return code from pthread_create() '-2-' is %d\n", rc);
               exit(-1);
            }
      }
       
      pthread_attr_destroy(&attr2);
       
      for(i=0; i<NUM_THREADS; i++)
      {
       pthread_join(threads2[i], &status2);
      }

    sum=data1.total_sum;


    free(data1.new_centers[0]);
    free(data1.new_centers);
    free(data1.new_centers_len);
    return data1.clusters;
}

