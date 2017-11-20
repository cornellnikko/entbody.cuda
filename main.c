#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>

#include "util.h"

#ifdef PLOT
#include "plot.h"
#endif

#include "io.h"

void simulate(int s);

#define NPOW        (7+7)
#define TPOW        (6 + NPOW%2)
#define N           (1 << NPOW)
#define NTHREADS    (1 << TPOW)
#define NBLOCKS     (1 << ((NPOW-TPOW)/2))

#define KICKFORCE 10.0f
#define FORCE_HERTZ_EPSILON 100.0f

#define ERROR_CHECK { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
    printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__);}}

#define CUDA_SAFE_CALL( call) {                                            \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#define CUT_DEVICE_INIT(dev) {                                               \
    int deviceCount;                                                         \
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));                        \
    if (deviceCount == 0) {                                                  \
        fprintf(stderr, "cutil error: no devices supporting CUDA.\n");       \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    if (dev < 0) dev = 0;                                                    \
    if (dev > deviceCount-1) dev = deviceCount - 1;                          \
    cudaDeviceProp deviceProp;                                               \
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));               \
    if (deviceProp.major < 1) {                                              \
        fprintf(stderr, "cutil error: device does not support CUDA.\n");     \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    CUDA_SAFE_CALL(cudaSetDevice(dev));                                      \
}


__device__ int mod_rvec(int a, int b, int p, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (p != 0){
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    *image = 0;
    return a;
}
__device__ float mymod(float a, float b){
  return a - b*(int)(a/b) + b*(a<0);
}
//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    int seed_in = 0;

    int device = 0;
    CUT_DEVICE_INIT(device);

    if (argc == 1)
        simulate(seed_in);
    else if (argc == 2)
        simulate(atoi(argv[1]));//seed_in);
    else {
        printf("usage:\n");
        printf("\t./entbody [seed]\n");
    }
    return 0;
}

__global__ void nbl_reset(int *cells, int *count, int size_total){
    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    int idy = blockDim.y*blockIdx.y + threadIdx.y;
    int i = idx + idy*NBLOCKS*NTHREADS;

    // reset the neighborlists
    if (i < size_total)
        count[i] = 0;
    if (i < size_total*NMAX)
        cells[i] = 0;
}

__global__ void nbl_build(float *x, int *cells, int *count, int *size,
        int size_total, float L){
    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    int idy = blockDim.y*blockIdx.y + threadIdx.y;
    int i = idx + idy*NBLOCKS*NTHREADS;
    volatile int index[2];

    index[0] = __float2int_rz(x[2*i+0]/L  * size[0]);
    index[1] = __float2int_rz(x[2*i+1]/L  * size[1]);
    volatile int t  = index[0] + index[1]*size[0];
    volatile unsigned int ct = atomicAdd(&count[t], 1);
    volatile unsigned int bt = NMAX*t + ct;
    cells[bt] = i;
}

__global__ void makecopy(float *x, float *copyx){
    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    int idy = blockDim.y*blockIdx.y + threadIdx.y;
    int i = idx + idy*NBLOCKS*NTHREADS;
    copyx[2*i+0] = x[2*i+0];
    copyx[2*i+1] = x[2*i+1];
}

//==================================================================
// the timestep - can be CPU or CUDA!
//==================================================================
__global__
void step(float *x, float *copyx, float *v, int *type, float *rad, float *col,
          int *cells, int *count, int *size, int size_total, int *key,
          float L, float R, int *pbc, float dt, float Tglobal, float colfact){

    int idx = blockDim.x*blockIdx.x + threadIdx.x;
    int idy = blockDim.y*blockIdx.y + threadIdx.y;
    int i = idx + idy*NBLOCKS*NTHREADS;
    int j;

    //=========================================
    // reset the neighborlists
    volatile int index[2];
    index[0] = __float2int_rz(copyx[2*i+0]/L  * size[0]);
    index[1] = __float2int_rz(copyx[2*i+1]/L  * size[1]);

    //==========================================
    // this is mainly for CUDA optimization
    int tt[2];
    int tix[2];
    int image[2];
    float dx[2];
    int goodcell, ind, tn;
    float dist;

    float px, py;
    float vx, vy;
    float fx, fy;
    float ox, oy;

    float trad;
    float tcol;
    float R2 = R*R;

    //==========================================
    // find forces on all particles
    tcol  = col[i]; trad = rad[i];
    px = copyx[2*i+0];
    py = copyx[2*i+1];
    vx = v[2*i+0];
    vy = v[2*i+1];

    fx = 0.0; fy = 0.0;
    ox = 0.0; oy = 0.0;

    #ifdef PLOT
    int ttype = type[i];
    if (key['w'] == 1){
        if (ttype == RED)
            oy = -KICKFORCE;
    }
    if (key['s'] == 1){
        if (ttype == RED)
            oy = KICKFORCE;
    }
    if (key['a'] == 1){
        if (ttype == RED)
            ox = -KICKFORCE;
    }
    if (key['d'] == 1){
        if (ttype == RED)
            ox = KICKFORCE;
    }
    #endif

    index[0] = (int)(px/L * size[0]);
    index[1] = (int)(py/L * size[1]);

    for (tt[0]=-1; tt[0]<=1; tt[0]++){
    for (tt[1]=-1; tt[1]<=1; tt[1]++){
        goodcell = 1;
        tix[0] = mod_rvec(index[0]+tt[0],size[0]-1,pbc[0],&image[0]);
        tix[1] = mod_rvec(index[1]+tt[1],size[1]-1,pbc[1],&image[1]);
        goodcell = !(pbc[0] < image[0] || pbc[1] < image[1]);

        if (goodcell){
            ind = tix[0] + tix[1]*size[0];
            for (j=0; j<count[ind]; j++){
                tn = cells[NMAX*ind+j];
                float px2 = copyx[2*tn+0];
                float py2 = copyx[2*tn+1];

                dist = 0.0;
                dx[0] = px2 - px;
                dx[0] += image[0]*L*tt[0];
                dist += dx[0]*dx[0];

                dx[1] = py2 - py;
                dx[1] += image[1]*L*tt[1];
                dist += dx[1]*dx[1];

                //===============================================
                // force calculation
                if (dist > EPSILON && dist < R2){
                    float r0 = trad+rad[tn];
                    float l = sqrt(dist);
                    float co = FORCE_HERTZ_EPSILON
                            * (1-l/r0)*(1-l/r0) * (l<r0);
                    fx += -co * dx[0];
                    fy += -co * dx[1];
                    tcol += co*co*(dx[0]*dx[0] + dx[1]*dx[1])/dist;
                }
             }
        }
    } }

    //=====================================
    // global forces
    fx += ox;
    fy += oy;

    float damp = 1.0;
    //float vlen = vx*vx + vy*vy;
    fx -= damp*vx;
    fy -= damp*vy;

    //=====================================
    // Newton-Stomer-Verlet
    vx += fx * dt;
    vy += fy * dt;

    px += vx * dt;
    py += vy * dt;

    //======================================
    // boundary conditions
    const float restoration = 1.0;
    if (pbc[0] == 1){
        if (px >= L-EPSILON || px < 0)
            px = mymod(px, L);
    }
    else {
        if (px >= L){px = 2*L-px; vx *= -restoration;}
        if (px < 0) {px = -px;    vx *= -restoration;}
        if (px >= L-EPSILON || px < 0){px = mymod(px, L);}
    }

    if (pbc[1] == 1){
        if (py >= L-EPSILON || py < 0)
            py = mymod(py, L);
    }
    else {
        if (py >= L){py = 2*L-py; vy *= -restoration;}
        if (py < 0) {py = -py;    vy *= -restoration;}
        if (py >= L-EPSILON || py < 0){py = mymod(py, L);}
    }

    tcol = tcol/colfact;

    col[i] = tcol;
    x[2*i+0] = px;  x[2*i+1] = py;
    v[2*i+0] = vx;  v[2*i+1] = vy;
}


//==================================================
// simulation
//==================================================
void simulate(int seed){
    ran_seed(seed);
#ifdef CUDA_NO_SM_11_ATOMIC_INTRINSICS
        printf("WARNING! Not using atomics!\n");
#endif
    int pbc[]     = {1,1};
    float L       = 0.0;
    float dt      = 1e-1;
    float t       = 0.0;
    float Tglobal = 0.0;

    float colfact = 16.0;

    int i;
    int mem_size_f = sizeof(float)*N;
    int mem_size_i = sizeof(int)*N;
    int mem_size_k = sizeof(int)*256;

    int *type    =   (int*)malloc(mem_size_i);
    float *rad   = (float*)malloc(mem_size_f);
    float *col   = (float*)malloc(mem_size_f);
    for (i=0; i<N; i++){ type[i] = 0; rad[i] = col[i] = 0.0;}

    float *x     = (float*)malloc(2*mem_size_f);
    float *v     = (float*)malloc(2*mem_size_f);
    float *copyx = (float*)malloc(2*mem_size_f);
    for (i=0; i<2*N; i++){x[i] = v[i] = copyx[i] = 0.0;}

    float time_end = 2e10;

	time_end = 10;

    #ifdef PLOT
        int *key;
        plot_init(900); // 2**18 - 450, 2**20 - 900, 2**21 - 1100
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #else
        int *key = (int*)malloc(mem_size_k);
        memset(key, 0, mem_size_k);
    #endif

    //==========================================
    // initialize
    float radius  = 1.0;
    L = 0.97*sqrt(pi*radius*radius*N);

    for (i=0; i<N; i++){
        rad[i] = radius;
        x[2*i+0] = L*ran_ran2();
        x[2*i+1] = L*ran_ran2();

        type[i] = BLACK;
        v[2*i+0] = 0.0;
        v[2*i+1] = 0.0;
        if (i==0) {
            type[i] = RED;
            rad[i] = 1*radius;
        }
     }


    // find out what happened in initialization
    float maxr = 0.0;
    for (i=0; i<N; i++)
        if (rad[i] > maxr) maxr = rad[i];
    float R = 2*maxr*1.0;

    // make boxes for the neighborlist
    int size[2];
    int size_total = 1;
    for (i=0; i<2; i++){
        size[i] = (int)(L / R);
        size_total *= size[i];
    }

    int *count = (int*)malloc(sizeof(int)*size_total);
    int *cells = (int*)malloc(sizeof(int)*size_total*NMAX);
    for (i=0; i<size_total; i++)
        count[i] = 0;
    for (i=0; i<size_total*NMAX; i++)
        cells[i] = 0;

    //==========================================================
    // where the magic happens
    //==========================================================
    int mem_size2 = sizeof(int)*2;
    int imem_size = sizeof(int)*N;
    int fmem_size = sizeof(float)*N;
    int fmem_siz2 = sizeof(float)*N*2;
    int mem_cell  = sizeof(int)*size_total;
    int mem_cell2 = sizeof(int)*size_total*NMAX;

    int *cu_count  = NULL;
    int *cu_cells  = NULL;
    int *cu_size   = NULL;
    int *cu_type   = NULL;
    int *cu_key    = NULL;
    float *cu_rad  = NULL;
    float *cu_col  = NULL;
    float *cu_x    = NULL;
    float *cu_copyx= NULL;
    float *cu_v    = NULL;
    int *cu_pbc    = NULL;

    cudaMalloc((void**) &cu_pbc,   2*sizeof(int));
    cudaMalloc((void**) &cu_count, mem_cell);
    cudaMalloc((void**) &cu_cells, mem_cell2);
    cudaMalloc((void**) &cu_size,  mem_size2);

    cudaMalloc((void**) &cu_key,   mem_size_k);
    cudaMalloc((void**) &cu_type,  imem_size);
    cudaMalloc((void**) &cu_rad,   fmem_size);
    cudaMalloc((void**) &cu_col,   fmem_size);
    cudaMalloc((void**) &cu_x,     fmem_siz2);
    cudaMalloc((void**) &cu_copyx, fmem_siz2);
    cudaMalloc((void**) &cu_v,     fmem_siz2);

    printf("Copying problem...\n");
    cudaMemcpy(cu_size,  size,  mem_size2, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_type,  type,  imem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_rad,   rad,   fmem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_col,   col,   fmem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_x,     x,     fmem_siz2, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_v,     v,     fmem_siz2, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_pbc, pbc,   2*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(cu_count, 0, mem_cell);
    cudaMemset(cu_cells, 0, mem_cell2);
    ERROR_CHECK

    int frames = 0;

    float rate;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    dim3 grid(NBLOCKS, NBLOCKS);
    dim3 block(NTHREADS, 1);

	FILE* fp = start_frames("particles.csv",N);
	printf("Ready to compute %f->%f\n",dt,time_end);


    for (t=0.0; t<time_end; t+=dt){
//	printf("Starting %f/%f\n",t,time_end);
        nbl_reset<<<grid, block>>>(cu_cells, cu_count, size_total);
        nbl_build<<<grid, block>>>(cu_x, cu_cells, cu_count, cu_size, size_total, L);
        makecopy<<<grid, block>>>(cu_x, cu_copyx);

        step<<<grid, block>>>(cu_x, cu_copyx, cu_v, cu_type, cu_rad, cu_col,
                    cu_cells, cu_count, cu_size, size_total, cu_key,
                    L, R, cu_pbc, dt, Tglobal, colfact);
        cudaThreadSynchronize();
        ERROR_CHECK

        if (frames % 10 == 0){
            clock_gettime(CLOCK_REALTIME, &end);
            rate = frames/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
            printf("N = %i\tframe = %i\trate = %f\r", N, frames, rate);
            fflush(stdout);
        }

        #ifdef PLOT
        if (frames % 1 == 0){
            cudaMemcpy(x, cu_x, fmem_siz2, cudaMemcpyDeviceToHost);
            cudaMemcpy(col, cu_col, fmem_size, cudaMemcpyDeviceToHost);

            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L, col);
            if (key['q'] == 1) break;

            cudaMemcpy(cu_key, key, mem_size_k, cudaMemcpyHostToDevice);
        }
        #endif

	//#ifdef IOFILE
	if (frames % 1 == 0)
	{
		cudaMemcpy(x, cu_x, fmem_siz2, cudaMemcpyDeviceToHost);
                cudaMemcpy(col, cu_col, fmem_size, cudaMemcpyDeviceToHost);
		write_frame(fp,x,rad,type,N,L,col);
	}
	//#endif
        frames++;
    }
    // end of the magic, cleanup
    //----------------------------------------------
	fflush(stdout);
	printf("\nEnd of magic, cleanup\n");

    free(cells);
    free(count);

    free(copyx);
    free(x);
    free(v);
    free(rad);
    free(type);
    free(col);


    cudaFree(cu_count);
    cudaFree(cu_cells);
    cudaFree(cu_key);
    cudaFree(cu_type);
    cudaFree(cu_rad);
    cudaFree(cu_col);
    cudaFree(cu_x);
    cudaFree(cu_v);
    cudaFree(cu_copyx);
    ERROR_CHECK

	printf("Passed end error check\n");
	end_frames(fp);
    #ifdef PLOT
    plot_clean();
    #else
    free(key);
    #endif
}





