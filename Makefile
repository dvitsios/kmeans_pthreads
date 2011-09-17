all:	cluster.c getopt.c kmeans.h unistd.h example.c getopt.h kmeans_clustering.c
	gcc -O3 example.c cluster.c getopt.c kmeans_clustering.c -o kmeans

clean:
	rm -f *.o *~
