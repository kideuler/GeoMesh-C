#include "GeoMesh.h"

#define NUM_POINTS 5
#define NUM_COMMANDS 2


int main(void){
    srand ( time(NULL) );
    int npoints = 100;
    struct DoubleMatrix xs = DoubleMatrix_create_Random(npoints,2,-10.0,10.0);
    struct Mesh msh = GeoMesh_Delaunay(&xs);
    //Mesh_print(&msh);

    FILE * gnuplotPipe = popen("gnuplot", "w");
    if (!gnuplotPipe) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output 'mesh.png'\n");
    fprintf(gnuplotPipe, "set xrange [-10:10] \n");
    fprintf(gnuplotPipe, "set yrange [-10:10] \n");
    fprintf(gnuplotPipe, "set pointsize 1.0\n" );
    fprintf(gnuplotPipe, "set style line 1 lc rgb 'red' pt 6 \n" );
    fprintf(gnuplotPipe, "plot '-' w p ls 1\n");
    int i;

    for (int i = 0; i < npoints; i++)
    {
        fprintf(gnuplotPipe, "%lf %lf\n", xs.data[i*xs.ncols], xs.data[i*xs.ncols+1]);
    }
    fprintf(gnuplotPipe, "e using 1:2 with lines\n");
    fprintf(gnuplotPipe, "replot\n");
    fprintf(gnuplotPipe, "exit\n");
    pclose(gnuplotPipe);
    return 0;
}