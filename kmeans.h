#ifndef _H_KMEANS
#define _H_KMEANS

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38
#endif

/* cluster.c */
int     cluster(int, int, float**, int, int, float***);

/* kmeans_clustering.c */
float **kmeans_clustering(float**, int, int, int, int);
float   euclid_dist_2        (float*, float*, int);
void     *find_nearest_point   (void *args);

#endif
